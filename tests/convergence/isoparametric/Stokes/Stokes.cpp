/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Taylor--Hood fields on exact and approximated curved domains. */

#include "../../CurvedGeometry.h"
#include "../../LiftedConvergence.h"
#include "../../SineMap.h"
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
  constexpr size_t NormOrder = 14, RefinedNormOrder = 18;

  struct LiftedErrors
  {
      LiftedErrors()
        : velocity(),
          pressure(),
          divergence()
      {}
      LiftedErrorNorm::Result velocity, pressure;
      std::array<Real, 3> divergence;
  };

#if defined(RODIN_CURVED_STOKES_PETSC) && defined(RODIN_USE_MPI)
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  /**
   * @brief Curved mixed workload reusing the established Stokes solver contracts.
   * @par Architecture
   * Prescribed-degree geometry is installed on a fresh grid and its traces. The
   * existing native/PETSc StokesProblem owns space and system construction,
   * residual checks, pressure gauge, and physical error integration. The
   * Volume and pressure mean are independently checked; x0 is preserved.
   * Solve-scoped const observers reuse the common exact-domain lift integrator.
   * MPI ownership and halo metadata are unchanged by geometry installation.
   */
  template <class ContextType, size_t Q = 2>
  class Workload
  {
    public:
      using Map = typename CurvedGeometry<Mesh<ContextType>>::Map;
      Workload(
        Polytope::Type geometry, size_t n, Map map = Map::Quadratic, bool lifted = false)
        : m_mesh(makeMesh(geometry, n)),
          m_reference(),
          m_geometry(m_mesh, map),
          m_sine(map == Map::Sine)
      {
        if (lifted)
          m_reference.emplace(m_mesh);
        m_geometry.template install<Q>();
      }

      const auto& getMesh() const
      {
        return m_mesh;
      }

      template <size_t K = 2>
      StokesErrors solve(StokesData::Field field, Real viscosity = 1,
        size_t order = AssemblyOrder, size_t normOrder = 0,
        LiftedErrors* lifted = nullptr) const
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const StokesData data(dim, field, m_sine ? dim - 1 : 1);
        const auto observe = [&](const auto& velocity, const auto& pressure,
                               const StokesData& exact) {
          if (!lifted)
            return;
          assert(m_reference);
          struct VelocityData
          {
              const StokesData& data;
              auto getSolution(const Math::SpatialPoint& x) const
              {
                return data.getVelocity(x);
              }
              auto getJacobian(const Math::SpatialPoint& x) const
              {
                return data.getVelocityJacobian(x);
              }
          };
          struct PressureData
          {
              const StokesData& data;
              Real getSolution(const Math::SpatialPoint& x) const
              {
                return data.getPressure(x);
              }
              auto getGradient(const Math::SpatialPoint& x) const
              {
                return data.getPressureGradient(x);
              }
          };
          std::array<Real, 3> squared{};
          const auto divergence = [&](const auto& derivatives, Real weight) {
            for (size_t component = 0; component < derivatives.size(); ++component)
            {
              Real trace = 0;
              for (size_t j = 0; j < dim; ++j)
                trace += derivatives[component](j, j);
              squared[component] += weight * trace * trace;
            }
          };
          const size_t integrationOrder = normOrder == 0 ? order : normOrder;
          lifted->velocity = LiftedErrorNorm::compute(*m_reference, m_mesh, velocity,
            VelocityData{exact}, SineMap(), integrationOrder, divergence);
          lifted->pressure = LiftedErrorNorm::compute(*m_reference, m_mesh, pressure,
            PressureData{exact}, SineMap(), integrationOrder);
#ifdef RODIN_USE_MPI
          if constexpr (requires { m_mesh.getShard(); })
            for (Real& value : squared)
              value = boost::mpi::all_reduce(
                m_mesh.getContext().getCommunicator(), value, std::plus<Real>());
#endif
          for (size_t component = 0; component < squared.size(); ++component)
            lifted->divergence[component] = std::sqrt(squared[component]);
        };
#ifdef RODIN_CURVED_STOKES_PETSC
        return PETScStokesProblem(m_mesh, data, order)
          .template solve<K>(viscosity, normOrder, observe);
#else
        return StokesProblem(m_mesh, data, order).solve<K>(viscosity, normOrder, observe);
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
      Optional<Mesh<ContextType>> m_reference;
      CurvedGeometry<Mesh<ContextType>> m_geometry;
      bool m_sine;
  };

  template <class ContextType, size_t Q = 2>
  class CurvedStokesTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      using Map = typename Workload<ContextType>::Map;
      static constexpr Real PressureGeometryTolerance = 1e-10;

      /** @brief Affine shear isolates geometry without fitting zero pressure errors.
       * The map changes the last coordinate but preserves x0 and volume.
       * Velocity is represented at K=max(2,Q); affine pressure is represented
       * at K-1. Field/pressure/divergence budgets are absolute and dimensionless.
       */
      void checkRepresentable(
        const StokesErrors& represented, const LiftedErrors& lifted) const
      {
        checkLifted(lifted);
        LiftedConvergence::expectRepresentable(
          represented.velocity, lifted.velocity, PatchTolerance);
        for (const auto& error : {represented.pressure, lifted.pressure.field,
               lifted.pressure.total})
        {
          EXPECT_TRUE(error.isFinite());
          EXPECT_LT(error.getL2(), PatchTolerance);
          EXPECT_LT(error.getH1Seminorm(), PatchTolerance);
        }
        EXPECT_TRUE(std::isfinite(represented.divergence));
        EXPECT_GE(represented.divergence, 0);
        EXPECT_LT(represented.divergence, PatchTolerance);
        EXPECT_LT(lifted.divergence[0], PatchTolerance);
        const Real divergenceBudget =
          std::sqrt(Real(UniformGrid(this->GetParam()).getDimension())) * PatchTolerance;
        EXPECT_NEAR(lifted.divergence[2], lifted.divergence[1], divergenceBudget);
      }

      void liftedAffineVelocityRates() const
      {
        constexpr size_t K = std::max(size_t(2), Q);
        LiftedConvergence velocity;
        for (size_t n : {3u, 5u, 9u})
        {
          SCOPED_TRACE(::testing::Message() << "geometry degree=" << Q
            << " velocity degree=" << K << " pressure degree=" << K - 1 << " n=" << n);
          Workload<ContextType, Q> problem(this->GetParam(), n, Map::Sine, true);
          LiftedErrors lifted;
          const auto represented = problem.template solve<K>(
            StokesData::Field::Affine, 1, AssemblyOrder, NormOrder, &lifted);
          checkRepresentable(represented, lifted);
          velocity.appendRepresentable(
            Real(1) / Real(n - 1), represented.velocity, lifted.velocity, PatchTolerance);
        }
        velocity.expectGeometryRates(Q);
      }

      void liftedQuadratureSensitivity() const
      {
        constexpr size_t K = std::max(size_t(2), Q);
        Workload<ContextType, Q> problem(this->GetParam(), 5, Map::Sine, true);
        std::array<LiftedErrors, 3> errors;
        for (size_t i = 0; i < errors.size(); ++i)
        {
          SCOPED_TRACE(::testing::Message() << "geometry degree=" << Q << " control=" << i);
          const auto represented = problem.template solve<K>(StokesData::Field::Affine,
            1, i == 1 ? RefinedOrder : AssemblyOrder,
            i == 2 ? RefinedNormOrder : NormOrder, &errors[i]);
          checkRepresentable(represented, errors[i]);
        }
        for (size_t i = 1; i < errors.size(); ++i)
        {
          LiftedConvergence::expectGeometrySensitivity(errors[0].velocity, errors[i].velocity);
          for (size_t component : {1u, 2u})
            EXPECT_NEAR(errors[0].divergence[component],
              errors[i].divergence[component], PatchTolerance);
        }
      }

      void checkLifted(const LiftedErrors& error) const
      {
        LiftedConvergence::expectDecomposition(error.velocity);
        LiftedConvergence::expectDecomposition(error.pressure);
        // The sine shear preserves x0 and unit volume. Pressure is a function
        // of x0 only, so its geometry defect is zero, not a positive rate study.
        EXPECT_LT(error.pressure.geometry.getL2(), PressureGeometryTolerance);
        EXPECT_LT(error.pressure.geometry.getH1Seminorm(), PressureGeometryTolerance);
        const auto velocity =
          std::array{error.velocity.field, error.velocity.geometry, error.velocity.total};
        const Real dimension = UniformGrid(this->GetParam()).getDimension();
        for (size_t component = 0; component < velocity.size(); ++component)
        {
          EXPECT_TRUE(std::isfinite(error.divergence[component]));
          EXPECT_GE(error.divergence[component], 0);
          EXPECT_LE(error.divergence[component],
            std::sqrt(dimension) * velocity[component].getH1Seminorm() +
              LiftedConvergence::RoundoffTolerance);
        }
      }

      void checkApproximatedRates() const
      {
        const auto levels = this->GetParam() == Polytope::Type::Tetrahedron
          ? std::initializer_list<size_t>{9, 11, 13}
          : this->GetParam() == Polytope::Type::Wedge
          ? std::initializer_list<size_t>{5, 7, 9}
          : UniformGrid(this->GetParam()).getDimension() == 3
          ? std::initializer_list<size_t>{3, 4, 5}
          : this->GetParam() == Polytope::Type::Triangle
          ? std::initializer_list<size_t>{9, 17, 33}
          : std::initializer_list<size_t>{3, 5, 9};
        LiftedConvergence velocity;
        std::array<ErrorHistory, 3> pressure;
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "n=" << n);
          Workload<ContextType> problem(this->GetParam(), n, Map::Sine, true);
          LiftedErrors lifted;
          const auto represented =
            problem.solve(StokesData::Field::Cubic, 1, AssemblyOrder, NormOrder, &lifted);
          checkLifted(lifted);
          const Real h = Real(1) / Real(n - 1);
          velocity.append(h, represented.velocity, lifted.velocity);
          const auto errors = std::array{
            represented.pressure, lifted.pressure.field, lifted.pressure.total};
          for (size_t component = 0; component < errors.size(); ++component)
          {
            ASSERT_TRUE(errors[component].isFinite());
            ASSERT_GT(errors[component].getL2(), 0);
            ASSERT_GT(errors[component].getH1Seminorm(), 0);
            pressure[component].append(h, errors[component]);
          }
          EXPECT_LE(represented.divergence,
            std::sqrt(Real(problem.getMesh().getSpaceDimension())) *
                represented.velocity.getH1Seminorm() +
              LiftedConvergence::RoundoffTolerance);
        }
        velocity.expectRates(2, 2);
        for (size_t component = 0; component < pressure.size(); ++component)
          for (size_t i = 1; i < pressure[component].getSize(); ++i)
          {
            SCOPED_TRACE(::testing::Message()
              << "pressure component=" << component << " interval=" << i);
            const auto& coarse = pressure[component].getSample(i - 1).error;
            const auto& fine = pressure[component].getSample(i).error;
            const auto rate = pressure[component].getAlgebraicRates(i);
            SCOPED_TRACE(::testing::Message()
              << "L2=" << coarse.getL2() << " -> " << fine.getL2()
              << " H1=" << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm()
              << " rates=" << rate.getL2() << "," << rate.getH1Seminorm());
            EXPECT_GT(coarse.getL2(), fine.getL2());
            EXPECT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
            EXPECT_GT(rate.getL2(), 2 - RateMargin);
            EXPECT_LT(rate.getL2(), 2 + RateMargin);
            EXPECT_GT(rate.getH1Seminorm(), 1 - DerivativeMargin);
            EXPECT_LT(rate.getH1Seminorm(), 1 + DerivativeMargin);
          }
      }

      void checkApproximatedSensitivity() const
      {
        Workload<ContextType> problem(this->GetParam(), 5, Map::Sine, true);
        std::vector<LiftedConvergence::Components> velocity, pressure;
        for (size_t i = 0; i < 3; ++i)
        {
          LiftedErrors lifted;
          const auto represented = problem.solve(StokesData::Field::Cubic, 1,
            i == 1 ? RefinedOrder : AssemblyOrder, i == 2 ? RefinedNormOrder : NormOrder,
            &lifted);
          checkLifted(lifted);
          velocity.push_back(
            LiftedConvergence::components(represented.velocity, lifted.velocity));
          pressure.push_back(
            LiftedConvergence::components(represented.pressure, lifted.pressure));
        }
        for (size_t i = 1; i < velocity.size(); ++i)
        {
          SCOPED_TRACE(::testing::Message() << "control=" << i);
          LiftedConvergence::expectSensitivity(velocity[0], velocity[i]);
          // Geometry-pressure norms are analytically zero: no relative comparison.
          for (size_t component : {0u, 1u, 3u})
            for (const auto& pair :
              {std::pair{pressure[0][component].getL2(), pressure[i][component].getL2()},
                std::pair{pressure[0][component].getH1Seminorm(),
                  pressure[i][component].getH1Seminorm()}})
            {
              ASSERT_GT(pair.first, 0);
              ASSERT_TRUE(std::isfinite(pair.second));
              EXPECT_LT(std::abs(pair.second / pair.first - 1), SensitivityTolerance);
            }
        }
      }

      void checkApproximatedPatch() const
      {
        Workload<ContextType> problem(this->GetParam(), 3, Map::Sine, true);
        LiftedErrors lifted;
        const auto represented =
          problem.solve(StokesData::Field::Affine, 1, AssemblyOrder, NormOrder, &lifted);
        checkLifted(lifted);
        for (const auto& error : {represented.velocity, represented.pressure,
               lifted.velocity.field, lifted.pressure.field, lifted.pressure.total})
        {
          EXPECT_LT(error.getL2(), PatchTolerance);
          EXPECT_LT(error.getH1Seminorm(), PatchTolerance);
        }
        EXPECT_LT(represented.divergence, PatchTolerance);
        EXPECT_LT(lifted.divergence[0], PatchTolerance);
        // This shear lift need not remain divergence-free on the exact domain.
        EXPECT_GT(lifted.velocity.geometry.getH1Seminorm(), 0);
        EXPECT_GT(lifted.divergence[1], 0);
      }

      void checkApproximatedControl() const
      {
        Workload<ContextType> problem(this->GetParam(), 5, Map::Sine, true);
        LiftedErrors base, wrong;
        problem.solve(StokesData::Field::Quadratic, 1, AssemblyOrder, NormOrder, &base);
        const auto incorrect = problem.solve(
          StokesData::Field::Quadratic, WrongViscosity, AssemblyOrder, NormOrder, &wrong);
        checkLifted(base);
        checkLifted(wrong);
        // Changing viscosity can be absorbed in pressure for this shear profile.
        // It need not increase velocity error, so pressure is the operator oracle.
        for (const auto& error :
          {incorrect.pressure, wrong.pressure.field, wrong.pressure.total})
        {
          EXPECT_GT(error.getL2(), WrongPressureL2);
          EXPECT_GT(error.getH1Seminorm(), WrongPressureH1);
        }
        EXPECT_EQ(base.velocity.geometry.getL2(), wrong.velocity.geometry.getL2());
        EXPECT_EQ(base.velocity.geometry.getH1Seminorm(),
          wrong.velocity.geometry.getH1Seminorm());
      }

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
        for (const auto& pair : {std::pair{baseline.velocity, refined.velocity},
               std::pair{baseline.pressure, refined.pressure}})
          for (const auto& norms : {std::pair{pair.first.getL2(), pair.second.getL2()},
                 std::pair{pair.first.getH1Seminorm(), pair.second.getH1Seminorm()}})
          {
            ASSERT_TRUE(std::isfinite(norms.first));
            ASSERT_TRUE(std::isfinite(norms.second));
            ASSERT_GT(norms.first, 0);
            EXPECT_LT(std::abs(norms.second / norms.first - 1), SensitivityTolerance);
          }
      }

      void checkGauge(Map map = Map::Quadratic) const
      {
        Workload<ContextType, Q> problem(this->GetParam(), 3, map);
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
          const Geometry::PolytopeQuadrature quadrature(*cell, qf);
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

  using LocalQ1Test = CurvedStokesTest<Context::Local, 1>;
  TEST_P(LocalQ1Test, LiftedAffineVelocityRates)
  {
    liftedAffineVelocityRates();
  }
  TEST_P(LocalQ1Test, LiftedQuadratureSensitivity)
  {
    liftedQuadratureSensitivity();
  }
  TEST_P(LocalQ1Test, LiftedPhysicalPressureGauge)
  {
    checkGauge(Map::Sine);
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, LocalQ1Test,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
  using LocalQ3Test = CurvedStokesTest<Context::Local, 3>;
  TEST_P(LocalQ3Test, LiftedAffineVelocityRates)
  {
    liftedAffineVelocityRates();
  }
  TEST_P(LocalQ3Test, LiftedQuadratureSensitivity)
  {
    liftedQuadratureSensitivity();
  }
  TEST_P(LocalQ3Test, LiftedPhysicalPressureGauge)
  {
    checkGauge(Map::Sine);
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, LocalQ3Test,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
#ifndef RODIN_CURVED_STOKES_PETSC
  /** @brief Fixed-mesh pressure forward-error regression, not a rate study.
   * The affine pressure is represented exactly on the q3 n5 wedge mesh.
   * A 1e-11 dimensionless absolute H1 budget separates the corrected solve
   * from the observed 2.47e-10 uncorrected error. The three-level rate study
   * retains its independent 1e-9 patch budget and n9 endpoint.
   */
  TEST(Rodin_Convergence_CurvedStokes, NativeCubicWedgePressureForwardAccuracy)
  {
    using Problem = Workload<Context::Local,3>;
    Problem problem(Polytope::Type::Wedge,5,Problem::Map::Sine,true);
    LiftedErrors lifted;
    const auto errors = problem.solve<3>(StokesData::Field::Affine,1,
      AssemblyOrder,NormOrder,&lifted);
    constexpr Real PressureForwardTolerance = 1e-11;
    for (const auto& pressure : {errors.pressure,lifted.pressure.field,lifted.pressure.total})
    {
      ASSERT_TRUE(pressure.isFinite());
      EXPECT_LT(pressure.getL2(),PressureForwardTolerance);
      EXPECT_LT(pressure.getH1Seminorm(),PressureForwardTolerance);
    }
  }
#endif
  using LocalTest = CurvedStokesTest<Context::Local>;
  TEST_P(LocalTest, ApproximatedPhysicalPressureGauge)
  {
    checkGauge(Map::Sine);
  }
  TEST_P(LocalTest, ApproximatedVelocityPressureRates)
  {
    checkApproximatedRates();
  }
  TEST_P(LocalTest, ApproximatedSensitivity)
  {
    checkApproximatedSensitivity();
  }
  TEST_P(LocalTest, ApproximatedAffinePatch)
  {
    checkApproximatedPatch();
  }
  TEST_P(LocalTest, ApproximatedRejectsWrongViscosity)
  {
    checkApproximatedControl();
  }
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
  using MPIQ1Test = CurvedStokesTest<Context::MPI, 1>;
  TEST_P(MPIQ1Test, LiftedAffineVelocityRates)
  {
    liftedAffineVelocityRates();
  }
  TEST_P(MPIQ1Test, LiftedQuadratureSensitivity)
  {
    liftedQuadratureSensitivity();
  }
  TEST_P(MPIQ1Test, LiftedPhysicalPressureGauge)
  {
    checkGauge(Map::Sine);
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, MPIQ1Test,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
  using MPIQ3Test = CurvedStokesTest<Context::MPI, 3>;
  TEST_P(MPIQ3Test, LiftedAffineVelocityRates)
  {
    liftedAffineVelocityRates();
  }
  TEST_P(MPIQ3Test, LiftedQuadratureSensitivity)
  {
    liftedQuadratureSensitivity();
  }
  TEST_P(MPIQ3Test, LiftedPhysicalPressureGauge)
  {
    checkGauge(Map::Sine);
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, MPIQ3Test,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
  using MPITest = CurvedStokesTest<Context::MPI>;
  TEST_P(MPITest, ApproximatedPhysicalPressureGauge)
  {
    checkGauge(Map::Sine);
  }
  TEST_P(MPITest, ApproximatedVelocityPressureRates)
  {
    checkApproximatedRates();
  }
  TEST_P(MPITest, ApproximatedSensitivity)
  {
    checkApproximatedSensitivity();
  }
  TEST_P(MPITest, ApproximatedAffinePatch)
  {
    checkApproximatedPatch();
  }
  TEST_P(MPITest, ApproximatedRejectsWrongViscosity)
  {
    checkApproximatedControl();
  }
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
