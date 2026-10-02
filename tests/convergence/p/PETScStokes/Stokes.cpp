/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local/MPI degree refinement of both Stokes fields. */

#include "../../PETScStokesProblem.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::P::PETScStokes
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  template <class ContextType>
  class PETScStokesPTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      auto makeMesh(size_t n) const
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(GetParam()).makeMesh(n);
#ifdef RODIN_USE_MPI
        else
          return DistributedUniformGrid(Context::MPI(*environment, *world), GetParam())
            .makeMesh(n);
#endif
      }

      void checkDegrees() const
      {
        auto mesh = makeMesh(3);
        const PETScStokesProblem problem(
          mesh, StokesData(mesh.getDimension(), StokesData::Field::Smooth));
        const auto e2 = problem.template solve<2>();
        const auto e3 = problem.template solve<3>();
        const auto e4 = problem.template solve<4>();
        ErrorHistory velocity, pressure;
        velocity.append(2, e2.velocity).append(3, e3.velocity).append(4, e4.velocity);
        pressure.append(1, e2.pressure).append(2, e3.pressure).append(3, e4.pressure);
        for (bool isVelocity : {true, false})
        {
          const auto& history = isVelocity ? velocity : pressure;
          ASSERT_EQ(history.getSize(), 3u);
          for (size_t i = 1; i < history.getSize(); ++i)
          {
            const auto& coarse = history.getSample(i - 1).error;
            const auto& fine = history.getSample(i).error;
            SCOPED_TRACE(::testing::Message()
              << "velocity=" << isVelocity << " interval=" << i << " L2 "
              << coarse.getL2() << " -> " << fine.getL2() << " H1 "
              << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm());
            ASSERT_TRUE(coarse.isFinite());
            ASSERT_TRUE(fine.isFinite());
            ASSERT_GT(fine.getL2(), 0);
            ASSERT_GT(fine.getH1Seminorm(), 0);
            ASSERT_GT(coarse.getL2(), fine.getL2());
            ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
            const auto decay = history.getExponentialRates(i);
            EXPECT_GT(decay.getL2(), 0.1);
            EXPECT_GT(decay.getH1Seminorm(), 0.1);
          }
        }
      }

      void checkPatches() const
      {
        auto mesh = makeMesh(3);
        const PETScStokesProblem quadratic(
          mesh, StokesData(mesh.getDimension(), StokesData::Field::Quadratic));
        const PETScStokesProblem cubic(
          mesh, StokesData(mesh.getDimension(), StokesData::Field::Cubic));
        const PETScStokesProblem quartic(
          mesh, StokesData(mesh.getDimension(), StokesData::Field::Quartic));
        for (const auto& errors : {quadratic.template solve<2>(),
               cubic.template solve<3>(), quartic.template solve<4>()})
        {
          EXPECT_LT(errors.velocity.getL2(), 1e-9);
          EXPECT_LT(errors.velocity.getH1Seminorm(), 1e-9);
          EXPECT_LT(errors.pressure.getL2(), 1e-9);
          EXPECT_LT(errors.pressure.getH1Seminorm(), 1e-9);
          EXPECT_LT(errors.divergence, 1e-9);
        }
      }

      void checkWrongViscosity() const
      {
        auto mesh = makeMesh(3);
        const PETScStokesProblem problem(
          mesh, StokesData(mesh.getDimension(), StokesData::Field::Quadratic));
        for (const auto& errors : {problem.template solve<2>(2),
               problem.template solve<3>(2), problem.template solve<4>(2)})
        {
          EXPECT_LT(errors.velocity.getL2(), 1e-9);
          EXPECT_LT(errors.velocity.getH1Seminorm(), 1e-9);
          EXPECT_GT(errors.pressure.getL2(), 0.1);
          EXPECT_GT(errors.pressure.getH1Seminorm(), 1);
        }
      }

      void checkQuadrature() const
      {
        auto mesh = makeMesh(3);
        const StokesData data(mesh.getDimension(), StokesData::Field::Smooth);
        const auto a = PETScStokesProblem(mesh, data).template solve<4>();
        const auto b = PETScStokesProblem(mesh, data, 18).template solve<4>();
        for (const auto& pair :
          {std::pair{a.velocity, b.velocity}, std::pair{a.pressure, b.pressure}})
        {
          ASSERT_TRUE(pair.first.isFinite());
          ASSERT_TRUE(pair.second.isFinite());
          ASSERT_GT(pair.first.getL2(), 0);
          ASSERT_GT(pair.first.getH1Seminorm(), 0);
          EXPECT_LT(std::abs(pair.second.getL2() / pair.first.getL2() - 1), 1e-6);
          EXPECT_LT(
            std::abs(pair.second.getH1Seminorm() / pair.first.getH1Seminorm() - 1), 1e-6);
        }
      }
  };

  using PETScStokesPLocalTest = PETScStokesPTest<Context::Local>;
  TEST_P(PETScStokesPLocalTest, BothFieldsImproveAtEveryDegree)
  {
    checkDegrees();
  }
  TEST_P(PETScStokesPLocalTest, EveryPairReproducesPolynomialPatches)
  {
    checkPatches();
  }
  TEST_P(PETScStokesPLocalTest, EveryPairRejectsIncorrectViscosity)
  {
    checkWrongViscosity();
  }
  TEST_P(PETScStokesPLocalTest, P4QuadratureSensitivity)
  {
    checkQuadrature();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScStokesPLocalTest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

#ifdef RODIN_USE_MPI
  using PETScStokesPMPITest = PETScStokesPTest<Context::MPI>;
  TEST_P(PETScStokesPMPITest, BothFieldsImproveAtEveryDegree)
  {
    checkDegrees();
  }
  TEST_P(PETScStokesPMPITest, EveryPairReproducesPolynomialPatches)
  {
    checkPatches();
  }
  TEST_P(PETScStokesPMPITest, EveryPairRejectsIncorrectViscosity)
  {
    checkWrongViscosity();
  }
  TEST_P(PETScStokesPMPITest, P4QuadratureSensitivity)
  {
    checkQuadrature();
  }
  TEST_P(PETScStokesPMPITest, P4NormsCountOwnedCellsOnce)
  {
    auto mesh = makeMesh(2);
    const auto dim = mesh.getDimension();
    H1<4, Math::SpatialVector<Real>, decltype(mesh)> velocitySpace(
      std::integral_constant<size_t, 4>{}, mesh, dim);
    H1<3, Real, decltype(mesh)> pressureSpace(std::integral_constant<size_t, 3>{}, mesh);
    PETSc::Variational::GridFunction u(velocitySpace);
    PETSc::Variational::GridFunction p(pressureSpace);
    const VectorFunction reference(dim, [dim](const Point& x) {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      value.setZero();
      value(0) = 2 * x(0);
      return value;
    });
    u = 0.5 * reference;
    p = Zero();
    const auto velocity = ErrorNorm::computeVector(
      mesh, u, reference,
      [dim](const Point&) {
        Math::SpatialMatrix<Real> j(
          static_cast<std::uint8_t>(dim), static_cast<std::uint8_t>(dim));
        j.setZero();
        j(0, 0) = 2;
        return j;
      },
      16);
    const auto pressure = ErrorNorm::compute(mesh, p, RealFunction(1), Zero(dim), 16);
    EXPECT_NEAR(velocity.getL2(), 1 / std::sqrt(Real(3)), 1e-12);
    EXPECT_NEAR(velocity.getH1Seminorm(), 1, 1e-12);
    EXPECT_NEAR(pressure.getL2(), 1, 1e-12);
    EXPECT_NEAR(pressure.getH1Seminorm(), 0, 1e-12);
    EXPECT_NEAR(ErrorNorm::computeDivergenceL2(mesh, u, 16), 1, 1e-12);
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScStokesPMPITest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
#endif
}

int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  // Centralized RHS/solution matches the existing PETSc Stokes h workload.
  for (const char* name : {"-mat_mumps_icntl_20", "-mat_mumps_icntl_21"})
  {
    PetscBool set = PETSC_FALSE;
    if (PetscOptionsHasName(nullptr, nullptr, name, &set) != PETSC_SUCCESS)
      return 1;
    if (!set && PetscOptionsSetValue(nullptr, name, "0") != PETSC_SUCCESS)
      return 1;
  }
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::P::PETScStokes::environment = &env;
    Rodin::Tests::Convergence::P::PETScStokes::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
