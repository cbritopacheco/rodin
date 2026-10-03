/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Taylor--Hood velocity, pressure and divergence on exact P2 maps. */

#include "../../CurvedGeometry.h"
#ifdef RODIN_CURVED_STOKES_PETSC
#include "../../PETScStokesProblem.h"
#else
#include "../../StokesProblem.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::Isoparametric::Stokes
{
  // Mapped nonpolynomial integrands: higher-order quadrature is checked separately.
  constexpr size_t AssemblyOrder = 12;
  constexpr size_t RefinedOrder = 16;
  constexpr Real PatchTolerance = 1e-9; // Dimensionless fields; roundoff allowance.
  constexpr Real RateMargin = 0.55; // Finite-resolution policy, not an error theorem.
  constexpr Real DerivativeMargin = 0.45;
  constexpr Real SensitivityTolerance = 1e-6; // Relative norm contamination budget.
  constexpr Real WrongViscosity = 2; // Exact manufactured sources use unit viscosity.
  constexpr Real WrongPressureL2 = 0.1;
  constexpr Real WrongPressureH1 = 1;
  constexpr Real GaugeTolerance = 1e-12; // Analytic unit volume and zero pressure mean.

#if defined(RODIN_CURVED_STOKES_PETSC) && defined(RODIN_USE_MPI)
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  /**
   * @brief Curved mixed workload reusing the established Stokes solver contracts.
   * @par Architecture
   * Exact P2 geometry is installed on a fresh grid and all its traces. The
   * existing native/PETSc StokesProblem owns space and system construction,
   * residual checks, pressure gauge, and physical error integration. The
   * map preserves volume and x0, so the analytic pressure remains mean zero.
   * MPI ownership and halo metadata are unchanged by geometry installation.
   */
  template <class ContextType>
  class Workload
  {
    public:
      Workload(Polytope::Type geometry, size_t n)
        : m_mesh(makeMesh(geometry, n)),
          m_geometry(m_mesh)
      {
        m_geometry.template install<2>();
      }

      const auto& getMesh() const
      {
        return m_mesh;
      }

      StokesErrors solve(
        StokesData::Field field, Real viscosity = 1, size_t order = AssemblyOrder) const
      {
        const StokesData data(m_mesh.getDimension(), field);
#ifdef RODIN_CURVED_STOKES_PETSC
        return PETScStokesProblem(m_mesh, data, order).template solve<2>(viscosity);
#else
        return StokesProblem(m_mesh, data, order).solve<2>(viscosity);
#endif
      }

    private:
      static Mesh<ContextType> makeMesh(Polytope::Type geometry, size_t n)
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(geometry).makeMesh(n);
#if defined(RODIN_CURVED_STOKES_PETSC) && defined(RODIN_USE_MPI)
        else
          return DistributedUniformGrid(Context::MPI(*environment, *world), geometry)
            .makeMesh(n);
#endif
      }
      Mesh<ContextType> m_mesh;
      CurvedGeometry<Mesh<ContextType>> m_geometry;
  };

  template <class ContextType>
  class CurvedStokesTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      void checkRates() const
      {
        const auto levels = UniformGrid(this->GetParam()).getDimension() == 3
          ? std::initializer_list<size_t>{3, 4, 5}
          : std::initializer_list<size_t>{3, 5, 9};
        ErrorHistory velocity, pressure;
        NormHistory divergence;
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "n=" << n);
          Workload<ContextType> problem(this->GetParam(), n);
          const auto error = problem.solve(StokesData::Field::Cubic);
          ASSERT_TRUE(std::isfinite(error.divergence));
          ASSERT_GE(error.divergence, 0);
          EXPECT_LE(error.divergence,
            std::sqrt(Real(problem.getMesh().getDimension())) *
              error.velocity.getH1Seminorm());
          const Real h = Real(1) / Real(n - 1);
          velocity.append(h, error.velocity);
          pressure.append(h, error.pressure);
          divergence.append(h, error.divergence);
        }
        ASSERT_EQ(velocity.getSize(), 3u);
        for (size_t i = 1; i < velocity.getSize(); ++i)
        {
          for (bool isVelocity : {true, false})
          {
            const auto& history = isVelocity ? velocity : pressure;
            const auto& coarse = history.getSample(i - 1).error;
            const auto& fine = history.getSample(i).error;
            ASSERT_TRUE(coarse.isFinite());
            ASSERT_TRUE(fine.isFinite());
            ASSERT_GT(fine.getL2(), 0);
            ASSERT_GT(fine.getH1Seminorm(), 0);
            ASSERT_GT(coarse.getL2(), fine.getL2());
            ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
            const auto rate = history.getAlgebraicRates(i);
            SCOPED_TRACE(::testing::Message()
              << "velocity=" << isVelocity << " interval=" << i << " L2 "
              << coarse.getL2() << " -> " << fine.getL2() << " H1 "
              << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm()
              << " rates=" << rate.getL2() << ", " << rate.getH1Seminorm());
            const Real l2Order = isVelocity ? 3 : 2;
            const Real derivativeOrder = l2Order - 1;
            EXPECT_GT(rate.getL2(), l2Order - RateMargin);
            EXPECT_LT(rate.getL2(), l2Order + RateMargin);
            EXPECT_GT(rate.getH1Seminorm(), derivativeOrder - DerivativeMargin);
            EXPECT_LT(rate.getH1Seminorm(), derivativeOrder + DerivativeMargin);
          }
          const Real coarse = divergence.getSample(i - 1).error;
          const Real fine = divergence.getSample(i).error;
          ASSERT_TRUE(std::isfinite(coarse));
          ASSERT_TRUE(std::isfinite(fine));
          ASSERT_GE(coarse, 0);
          ASSERT_GE(fine, 0);
          // Divergence need not have an exact asymptotic order: cancellation
          // may make it smaller than the velocity H1 error. Its bound follows
          // from |tr(D(u-uh))| <= sqrt(d) |D(u-uh)|_F, since div(u)=0.
          const Real bound =
            std::sqrt(Real(UniformGrid(this->GetParam()).getDimension())) *
            velocity.getSample(i).error.getH1Seminorm();
          EXPECT_LE(fine, bound);
        }
      }

      void checkPatch() const
      {
        Workload<ContextType> problem(this->GetParam(), 3);
        const auto error = problem.solve(StokesData::Field::Affine);
        EXPECT_LT(error.velocity.getL2(), PatchTolerance);
        EXPECT_LT(error.velocity.getH1Seminorm(), PatchTolerance);
        EXPECT_LT(error.pressure.getL2(), PatchTolerance);
        EXPECT_LT(error.pressure.getH1Seminorm(), PatchTolerance);
        EXPECT_LT(error.divergence, PatchTolerance);
      }

      void checkWrongViscosity() const
      {
        Workload<ContextType> problem(this->GetParam(), 3);
        const auto error = problem.solve(StokesData::Field::Quadratic, WrongViscosity);
        ASSERT_TRUE(error.pressure.isFinite());
        EXPECT_GT(error.pressure.getL2(), WrongPressureL2);
        EXPECT_GT(error.pressure.getH1Seminorm(), WrongPressureH1);
      }

      void checkSensitivity() const
      {
        Workload<ContextType> problem(this->GetParam(), 3);
        const auto baseline = problem.solve(StokesData::Field::Cubic);
        const auto refined = problem.solve(StokesData::Field::Cubic, 1, RefinedOrder);
        for (const auto pair : {std::pair{baseline.velocity, refined.velocity},
               std::pair{baseline.pressure, refined.pressure}})
          for (const auto norms : {std::pair{pair.first.getL2(), pair.second.getL2()},
                 std::pair{pair.first.getH1Seminorm(), pair.second.getH1Seminorm()}})
          {
            ASSERT_TRUE(std::isfinite(norms.first));
            ASSERT_TRUE(std::isfinite(norms.second));
            ASSERT_GT(norms.first, 0);
            EXPECT_LT(std::abs(norms.second / norms.first - 1), SensitivityTolerance);
          }
      }

      void checkGauge() const
      {
        Workload<ContextType> problem(this->GetParam(), 3);
        const StokesData data(problem.getMesh().getDimension(), StokesData::Field::Cubic);
        const auto& mesh = problem.getMesh();
        const auto pressure = data.getPressure();
        Real volume = 0, mean = 0;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          if constexpr (requires { mesh.getShard(); })
            if (!mesh.getShard().isOwned(mesh.getDimension(), cell->getIndex()))
              continue;
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(AssemblyOrder, cell->getGeometry());
          const auto& quadrature = cell->getQuadrature(qf);
          for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
          {
            const auto& point = quadrature.getPoint(qp);
            const Real weight = qf.getWeight(qp) * point.getDistortion();
            volume += weight;
            mean += weight * pressure(point);
          }
        }
#if defined(RODIN_CURVED_STOKES_PETSC) && defined(RODIN_USE_MPI)
        if constexpr (std::is_same_v<ContextType, Context::MPI>)
        {
          const auto& comm = mesh.getContext().getCommunicator();
          volume = boost::mpi::all_reduce(comm, volume, std::plus<Real>());
          mean = boost::mpi::all_reduce(comm, mean, std::plus<Real>());
        }
#endif
        EXPECT_NEAR(volume, 1, GaugeTolerance);
        EXPECT_NEAR(mean, 0, GaugeTolerance);
      }
  };

  using LocalTest = CurvedStokesTest<Context::Local>;
  TEST_P(LocalTest, AffinePhysicalPatch)
  {
    checkPatch();
  }
  TEST_P(LocalTest, OptimalVelocityPressureRates)
  {
    checkRates();
  }
  TEST_P(LocalTest, RejectsWrongViscosity)
  {
    checkWrongViscosity();
  }
  TEST_P(LocalTest, QuadratureSensitivity)
  {
    checkSensitivity();
  }
  TEST_P(LocalTest, PhysicalPressureGauge)
  {
    checkGauge();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, LocalTest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

#if defined(RODIN_CURVED_STOKES_PETSC) && defined(RODIN_USE_MPI)
  using MPITest = CurvedStokesTest<Context::MPI>;
  TEST_P(MPITest, AffinePhysicalPatch)
  {
    checkPatch();
  }
  TEST_P(MPITest, OptimalVelocityPressureRates)
  {
    checkRates();
  }
  TEST_P(MPITest, RejectsWrongViscosity)
  {
    checkWrongViscosity();
  }
  TEST_P(MPITest, QuadratureSensitivity)
  {
    checkSensitivity();
  }
  TEST_P(MPITest, PhysicalPressureGauge)
  {
    checkGauge();
  }
  TEST_P(MPITest, DivergenceNormCountsOwnedCurvedCellsOnce)
  {
    Workload<Context::MPI> problem(GetParam(), 3);
    const auto& mesh = problem.getMesh();
    H1<2, Math::SpatialVector<Real>, Mesh<Context::MPI>> space(
      std::integral_constant<size_t, 2>{}, mesh, mesh.getDimension());
    PETSc::Variational::GridFunction u(space);
    u = VectorFunction(mesh.getDimension(), [dim = mesh.getDimension()](const Point& x) {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      value.setZero();
      value(0) = x(0);
      return value;
    });
    EXPECT_NEAR(
      ErrorNorm::computeDivergenceL2(mesh, u, AssemblyOrder), 1, GaugeTolerance);
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, MPITest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
#endif
}

#ifdef RODIN_CURVED_STOKES_PETSC
int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::Isoparametric::Stokes::environment = &env;
    Rodin::Tests::Convergence::Isoparametric::Stokes::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
#endif
