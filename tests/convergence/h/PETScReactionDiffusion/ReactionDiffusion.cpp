/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local/MPI coupled reaction-diffusion patches and h-rates. */

#include <array>
#include <gtest/gtest.h>
#include <petsc.h>

#include "Convergence.h"
#include "../../PETScReactionDiffusionProblem.h"

#ifdef RODIN_USE_MPI
#include <boost/mpi/environment.hpp>
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::PETScReactionDiffusion
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

  template <size_t K, class MeshType>
  std::array<ErrorNorms, 2> solve(const MeshType& mesh,
    ReactionDiffusionData::Field field, bool omitCoupling = false,
    size_t quadratureOrder = 12, Real tolerance = 1e-13)
  {
    return PETScReactionDiffusionProblem<K, MeshType>(mesh, field, quadratureOrder)
      .solve(omitCoupling, tolerance);
  }

  template <class ContextType, size_t K>
  void checkPatch(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 5);
    const auto field = K == 1 ? ReactionDiffusionData::Field::Affine
                              : ReactionDiffusionData::Field::Quadratic;
    for (const auto& error : solve<K>(mesh, field))
    {
      EXPECT_LT(error.getL2(), 1e-9);
      EXPECT_LT(error.getH1Seminorm(), 1e-9);
    }
  }

  template <class ContextType, size_t K>
  void checkRates(Polytope::Type geometry)
  {
    UniformGridHierarchy hierarchy(geometry,
      K == 1 ? std::initializer_list<size_t>{5, 9, 17}
             : std::initializer_list<size_t>{3, 5, 9});
    std::array<ErrorHistory, 2> histories;
    for (size_t level : hierarchy.getLevels())
    {
      SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << level);
      auto mesh = makeMesh<ContextType>(geometry, level);
      const auto errors = solve<K>(mesh, ReactionDiffusionData::Field::Smooth);
      for (size_t component = 0; component < 2; ++component)
        histories[component].append(hierarchy.getMeshSize(level), errors[component]);
    }
    for (size_t component = 0; component < 2; ++component)
    {
      SCOPED_TRACE(::testing::Message() << "component=" << component);
      const auto& history = histories[component];
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
        EXPECT_GT(rate.getL2(), Real(K) + 0.6);
        EXPECT_LT(rate.getL2(), Real(K) + 1.4);
        EXPECT_GT(rate.getH1Seminorm(), Real(K) - 0.25);
        EXPECT_LT(rate.getH1Seminorm(), Real(K) + 0.4);
      }
    }
  }

  template <class ContextType>
  void checkNegativeControl(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 5);
    for (const auto& error : solve<1>(mesh, ReactionDiffusionData::Field::Affine, true))
    {
      EXPECT_GT(error.getL2(), 1e-3);
      EXPECT_GT(error.getH1Seminorm(), 1e-2);
    }
  }

  template <class ContextType>
  void checkSensitivity(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 5);
    const auto baseline = solve<2>(mesh, ReactionDiffusionData::Field::Smooth);
    const auto refined =
      solve<2>(mesh, ReactionDiffusionData::Field::Smooth, false, 14, 1e-14);
    for (size_t component = 0; component < 2; ++component)
    {
      ASSERT_GT(baseline[component].getL2(), 0);
      ASSERT_GT(baseline[component].getH1Seminorm(), 0);
      EXPECT_LT(
        std::abs(refined[component].getL2() / baseline[component].getL2() - 1), 1e-6);
      EXPECT_LT(
        std::abs(
          refined[component].getH1Seminorm() / baseline[component].getH1Seminorm() - 1),
        1e-6);
    }
  }

  class PETScReactionDiffusionLocalTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(PETScReactionDiffusionLocalTest, AffineP1Patch)
  {
    checkPatch<Context::Local, 1>(GetParam());
  }
  TEST_P(PETScReactionDiffusionLocalTest, QuadraticP2Patch)
  {
    checkPatch<Context::Local, 2>(GetParam());
  }
  TEST_P(PETScReactionDiffusionLocalTest, P1OptimalRates)
  {
    checkRates<Context::Local, 1>(GetParam());
  }
  TEST_P(PETScReactionDiffusionLocalTest, P2OptimalRates)
  {
    checkRates<Context::Local, 2>(GetParam());
  }
  TEST_P(PETScReactionDiffusionLocalTest, AffinePatchRejectsOmittedCoupling)
  {
    checkNegativeControl<Context::Local>(GetParam());
  }
  TEST_P(PETScReactionDiffusionLocalTest, QuadratureAndSolverSensitivity)
  {
    checkSensitivity<Context::Local>(GetParam());
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScReactionDiffusionLocalTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

#ifdef RODIN_USE_MPI
  class PETScReactionDiffusionMPITest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(PETScReactionDiffusionMPITest, AffineP1Patch)
  {
    checkPatch<Context::MPI, 1>(GetParam());
  }
  TEST_P(PETScReactionDiffusionMPITest, QuadraticP2Patch)
  {
    checkPatch<Context::MPI, 2>(GetParam());
  }
  TEST_P(PETScReactionDiffusionMPITest, P1OptimalRates)
  {
    checkRates<Context::MPI, 1>(GetParam());
  }
  TEST_P(PETScReactionDiffusionMPITest, P2OptimalRates)
  {
    checkRates<Context::MPI, 2>(GetParam());
  }
  TEST_P(PETScReactionDiffusionMPITest, AffinePatchRejectsOmittedCoupling)
  {
    checkNegativeControl<Context::MPI>(GetParam());
  }
  TEST_P(PETScReactionDiffusionMPITest, QuadratureAndSolverSensitivity)
  {
    checkSensitivity<Context::MPI>(GetParam());
  }

  TEST_P(PETScReactionDiffusionMPITest, GlobalNormCountsOwnedCellsOnce)
  {
    auto mesh = makeMesh<Context::MPI>(GetParam(), 5);
    H1<1, Real, Mesh<Context::MPI>> space(std::integral_constant<size_t, 1>{}, mesh);
    PETSc::Variational::TrialFunction u(space), w(space);
    u.getSolution() = RealFunction(0);
    w.getSolution() = RealFunction(0);
    const VectorFunction gradient(mesh.getDimension(), [](const Point& p) {
      return Math::SpatialVector<Real>::Zero(p.getCoordinates().size());
    });
    for (const auto* solution : {&u.getSolution(), &w.getSolution()})
    {
      const auto error =
        ErrorNorm::compute(mesh, *solution, RealFunction(1), gradient, 12);
      EXPECT_NEAR(error.getL2(), 1, 1e-12);
      EXPECT_NEAR(ErrorNorm::computeL2(mesh, *solution, RealFunction(1), 12), 1, 1e-12);
      EXPECT_EQ(error.getH1Seminorm(), 0);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScReactionDiffusionMPITest,
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
    Rodin::Tests::Convergence::H::PETScReactionDiffusion::environment = &env;
    Rodin::Tests::Convergence::H::PETScReactionDiffusion::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
