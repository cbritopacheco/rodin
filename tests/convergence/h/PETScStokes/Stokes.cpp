/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local/MPI Taylor--Hood patches and h-convergence. */

#include <gtest/gtest.h>
#include "../../PETScStokesProblem.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::PETScStokes
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  template <class ContextType>
  auto makeMesh(Polytope::Type geometry, size_t n)
  {
    if constexpr (std::is_same_v<ContextType, Context::Local>)
      return UniformGrid(geometry).makeMesh(n);
#ifdef RODIN_USE_MPI
    else
      return DistributedUniformGrid(Context::MPI(*environment, *world), geometry)
        .makeMesh(n);
#endif
  }

  template <class MeshType>
  StokesErrors solve(
    const MeshType& mesh, StokesData::Field field, Real viscosity = 1, size_t order = 12)
  {
    return PETScStokesProblem(mesh, StokesData(mesh.getDimension(), field), order)
      .template solve<2>(viscosity);
  }

  template <class ContextType>
  void checkPatch(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 3);
    for (auto field : {StokesData::Field::Affine, StokesData::Field::Quadratic})
    {
      const auto errors = solve(mesh, field);
      EXPECT_LT(errors.velocity.getL2(), 1e-9);
      EXPECT_LT(errors.velocity.getH1Seminorm(), 1e-9);
      EXPECT_LT(errors.pressure.getL2(), 1e-9);
      EXPECT_LT(errors.pressure.getH1Seminorm(), 1e-9);
      EXPECT_LT(errors.divergence, 1e-9);
    }
  }

  template <class ContextType>
  void checkRates(Polytope::Type geometry)
  {
    ErrorHistory velocity, pressure;
    for (size_t n : {3, 5, 9})
    {
      SCOPED_TRACE(::testing::Message() << "n=" << n);
      auto mesh = makeMesh<ContextType>(geometry, n);
      const auto errors = solve(mesh, StokesData::Field::Cubic);
      velocity.append(1 / Real(n - 1), errors.velocity);
      pressure.append(1 / Real(n - 1), errors.pressure);
    }
    for (size_t i = 1; i < 3; ++i)
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
          << "velocity=" << isVelocity << " interval=" << i << " L2 " << coarse.getL2()
          << " -> " << fine.getL2() << " H1 " << coarse.getH1Seminorm() << " -> "
          << fine.getH1Seminorm() << " rates " << rate.getL2() << ", "
          << rate.getH1Seminorm());
        EXPECT_GT(rate.getL2(), isVelocity ? 2.45 : 1.45);
        EXPECT_LT(rate.getL2(), isVelocity ? 3.55 : 2.55);
        EXPECT_GT(rate.getH1Seminorm(), isVelocity ? 1.55 : 0.55);
        EXPECT_LT(rate.getH1Seminorm(), isVelocity ? 2.45 : 1.45);
      }
  }

  template <class ContextType>
  void checkNegative(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 3);
    const auto error = solve(mesh, StokesData::Field::Quadratic, 2);
    EXPECT_GT(error.pressure.getL2(), 0.1);
    EXPECT_GT(error.pressure.getH1Seminorm(), 1);
  }

  template <class ContextType>
  void checkSensitivity(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 3);
    const auto a = solve(mesh, StokesData::Field::Cubic);
    const auto b = solve(mesh, StokesData::Field::Cubic, 1, 14);
    for (const auto& pair :
      {std::pair{a.velocity, b.velocity}, std::pair{a.pressure, b.pressure}})
    {
      ASSERT_GT(pair.first.getL2(), 0);
      ASSERT_GT(pair.first.getH1Seminorm(), 0);
      EXPECT_LT(std::abs(pair.second.getL2() / pair.first.getL2() - 1), 1e-8);
      EXPECT_LT(
        std::abs(pair.second.getH1Seminorm() / pair.first.getH1Seminorm() - 1), 1e-8);
    }
  }

  class PETScStokesLocalTest : public ::testing::TestWithParam<Polytope::Type>
  {};
  TEST_P(PETScStokesLocalTest, Patches)
  {
    checkPatch<Context::Local>(GetParam());
  }
  TEST_P(PETScStokesLocalTest, OptimalRates)
  {
    checkRates<Context::Local>(GetParam());
  }
  TEST_P(PETScStokesLocalTest, RejectsWrongViscosity)
  {
    checkNegative<Context::Local>(GetParam());
  }
  TEST_P(PETScStokesLocalTest, QuadratureSensitivity)
  {
    checkSensitivity<Context::Local>(GetParam());
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScStokesLocalTest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

#ifdef RODIN_USE_MPI
  class PETScStokesMPITest : public ::testing::TestWithParam<Polytope::Type>
  {};
  TEST_P(PETScStokesMPITest, Patches)
  {
    checkPatch<Context::MPI>(GetParam());
  }
  TEST_P(PETScStokesMPITest, OptimalRates)
  {
    checkRates<Context::MPI>(GetParam());
  }
  TEST_P(PETScStokesMPITest, RejectsWrongViscosity)
  {
    checkNegative<Context::MPI>(GetParam());
  }
  TEST_P(PETScStokesMPITest, QuadratureSensitivity)
  {
    checkSensitivity<Context::MPI>(GetParam());
  }
  TEST_P(PETScStokesMPITest, DivergenceNormCountsOwnedCellsOnce)
  {
    auto mesh = makeMesh<Context::MPI>(GetParam(), 2);
    H1<2, Math::SpatialVector<Real>, decltype(mesh)> space(
      std::integral_constant<size_t, 2>{}, mesh, mesh.getDimension());
    PETSc::Variational::GridFunction u(space);
    u = VectorFunction(mesh.getDimension(), [dim = mesh.getDimension()](const Point& x) {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      value.setZero();
      value(0) = x(0);
      return value;
    });
    EXPECT_NEAR(ErrorNorm::computeDivergenceL2(mesh, u, 12), 1, 1e-12);
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScStokesMPITest,
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
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::H::PETScStokes::environment = &env;
    Rodin::Tests::Convergence::H::PETScStokes::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
