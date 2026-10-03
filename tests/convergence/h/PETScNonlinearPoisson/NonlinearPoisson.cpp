/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc/SNES local and distributed semilinear Poisson h-rates. */

#include <gtest/gtest.h>
#include "Convergence.h"
#include "../../PETScNonlinearPoisson.h"
#ifdef RODIN_USE_MPI
#include <boost/mpi/environment.hpp>
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::H::PETScNonlinearPoisson
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  template <class ContextType>
  auto makeMesh(Polytope::Type geometry, size_t level)
  {
    if constexpr (std::is_same_v<ContextType, Context::Local>)
      return UniformGrid(geometry).makeMesh(level);
#ifdef RODIN_USE_MPI
    else
    {
      Context::MPI context(*environment, *world);
      return DistributedUniformGrid(context, geometry).makeMesh(level);
    }
#endif
  }

  template <class ContextType, size_t K>
  void checkRates(Polytope::Type geometry)
  {
    UniformGridHierarchy hierarchy(geometry,
      K == 1 ? std::initializer_list<size_t>{5, 9, 17}
             : std::initializer_list<size_t>{3, 5, 9});
    ErrorHistory history;
    for (size_t level : hierarchy.getLevels())
    {
      SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << level);
      auto mesh = makeMesh<ContextType>(geometry, level);
      history.append(hierarchy.getMeshSize(level),
        PETScNonlinearPoissonProblem<K, decltype(mesh)>(mesh).solve());
    }
    ASSERT_EQ(history.getSize(), 3u);
    for (size_t i = 1; i < history.getSize(); ++i)
    {
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
        << "L2 " << coarse.getL2() << " -> " << fine.getL2() << ", H1 "
        << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm() << ", rates "
        << rate.getL2() << ", " << rate.getH1Seminorm());
      EXPECT_GT(rate.getL2(), Real(K) + 0.5);
      EXPECT_LT(rate.getL2(), Real(K) + 1.5);
      EXPECT_GT(rate.getH1Seminorm(), Real(K) - 0.3);
      EXPECT_LT(rate.getH1Seminorm(), Real(K) + 0.5);
    }
  }

  template <class ContextType>
  void checkTangent(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 3);
    EXPECT_LT(
      (PETScNonlinearPoissonProblem<1, decltype(mesh)>(mesh).tangentDefect()), 1e-6);
    EXPECT_LT(
      (PETScNonlinearPoissonProblem<2, decltype(mesh)>(mesh).tangentDefect()), 1e-6);
    EXPECT_GT((PETScNonlinearPoissonProblem<2, decltype(mesh)>(mesh, 12, 1, false, true)
                  .tangentDefect()),
      1e-3);
  }

  template <class ContextType>
  void checkSensitivity(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 3);
    const auto baseline = PETScNonlinearPoissonProblem<2, decltype(mesh)>(mesh).solve();
    const auto refined =
      PETScNonlinearPoissonProblem<2, decltype(mesh)>(mesh, 14).solve(1e-12);
    ASSERT_GT(baseline.getL2(), 0);
    ASSERT_GT(baseline.getH1Seminorm(), 0);
    EXPECT_LT(std::abs(refined.getL2() / baseline.getL2() - 1), 1e-6);
    EXPECT_LT(std::abs(refined.getH1Seminorm() / baseline.getH1Seminorm() - 1), 1e-6);
  }

  template <class ContextType>
  void checkNegativeControl(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 9);
    const auto correct =
      PETScNonlinearPoissonProblem<2, decltype(mesh)>(mesh, 12, 4).solve();
    const auto wrong =
      PETScNonlinearPoissonProblem<2, decltype(mesh)>(mesh, 12, 4, true).solve();
    EXPECT_LT(correct.getL2(), 0.05);
    EXPECT_LT(correct.getH1Seminorm(), 0.2);
    EXPECT_GT(wrong.getL2(), 0.05);
    EXPECT_GT(wrong.getH1Seminorm(), 0.2);
  }

  class PETScNonlinearPoissonLocalTest : public ::testing::TestWithParam<Polytope::Type>
  {};
  TEST_P(PETScNonlinearPoissonLocalTest, P1OptimalRates)
  {
    checkRates<Context::Local, 1>(GetParam());
  }
  TEST_P(PETScNonlinearPoissonLocalTest, P2OptimalRates)
  {
    checkRates<Context::Local, 2>(GetParam());
  }
  TEST_P(PETScNonlinearPoissonLocalTest, SNESResidualTangentConsistency)
  {
    checkTangent<Context::Local>(GetParam());
  }
  TEST_P(PETScNonlinearPoissonLocalTest, QuadratureAndSolverSensitivity)
  {
    checkSensitivity<Context::Local>(GetParam());
  }
  TEST_P(PETScNonlinearPoissonLocalTest, RejectsOmittedCubicReaction)
  {
    checkNegativeControl<Context::Local>(GetParam());
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScNonlinearPoissonLocalTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

#ifdef RODIN_USE_MPI
  class PETScNonlinearPoissonMPITest : public ::testing::TestWithParam<Polytope::Type>
  {};
  TEST_P(PETScNonlinearPoissonMPITest, P1OptimalRates)
  {
    checkRates<Context::MPI, 1>(GetParam());
  }
  TEST_P(PETScNonlinearPoissonMPITest, P2OptimalRates)
  {
    checkRates<Context::MPI, 2>(GetParam());
  }
  TEST_P(PETScNonlinearPoissonMPITest, SNESResidualTangentConsistency)
  {
    checkTangent<Context::MPI>(GetParam());
  }
  TEST_P(PETScNonlinearPoissonMPITest, QuadratureAndSolverSensitivity)
  {
    checkSensitivity<Context::MPI>(GetParam());
  }
  TEST_P(PETScNonlinearPoissonMPITest, RejectsOmittedCubicReaction)
  {
    checkNegativeControl<Context::MPI>(GetParam());
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScNonlinearPoissonMPITest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
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
    Rodin::Tests::Convergence::H::PETScNonlinearPoisson::environment = &env;
    Rodin::Tests::Convergence::H::PETScNonlinearPoisson::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  return PetscFinalize() == PETSC_SUCCESS ? result : 1;
}
