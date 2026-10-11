/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local/MPI variable-conductivity patches and h-rates. */

#include <gtest/gtest.h>
#include <petsc.h>

#include "Convergence.h"
#include "../../Conductivity.h"
#include "Rodin/Assembly.h"
#include "Rodin/PETSc.h"

#ifdef RODIN_USE_MPI
#include <boost/mpi/environment.hpp>
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::PETScConductivity
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
  ErrorNorms solve(const MeshType& mesh, ConductivityData::Field field,
    bool omitCoefficient = false, size_t quadratureOrder = 12, Real tolerance = 1e-13)
  {
    const ConductivityData data(mesh.getDimension(), field);
    const auto exact = data.getSolution();
    const auto source = data.getSource();
    const auto gradient = data.getGradient();
    const auto gamma = data.getCoefficient(omitCoefficient);
    H1<K, Real, MeshType> space(std::integral_constant<size_t, K>{}, mesh);
    PETSc::Variational::TrialFunction u(space);
    PETSc::Variational::TestFunction v(space);
    auto stiffness = Integral(gamma * Grad(u), Grad(v));
    auto load = Integral(source, v);
    stiffness.setOrder(quadratureOrder);
    load.setOrder(quadratureOrder);
    Problem problem(u, v);
    problem = stiffness - load + DirichletBC(u, exact);
    PETSc::Solver::CG solver(problem);
    solver.setTolerances(tolerance, 1e-14, 1e5, 20000);
    solver.solve();
    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    EXPECT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
    EXPECT_GT(reason, 0);
    EXPECT_TRUE(std::isfinite(solver.getError()));
    EXPECT_LT(solver.getError(), 1e-9);
    return ErrorNorm::compute(mesh, u.getSolution(), exact, gradient, quadratureOrder);
  }

  template <class ContextType, size_t K>
  void checkPatch(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 5);
    const auto field =
      K == 1 ? ConductivityData::Field::Affine : ConductivityData::Field::Quadratic;
    const auto error = solve<K>(mesh, field);
    EXPECT_LT(error.getL2(), 1e-9);
    EXPECT_LT(error.getH1Seminorm(), 1e-9);
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
      history.append(
        hierarchy.getMeshSize(level), solve<K>(mesh, ConductivityData::Field::Smooth));
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
      EXPECT_GT(rate.getL2(), Real(K) + 0.6);
      EXPECT_LT(rate.getL2(), Real(K) + 1.4);
      EXPECT_GT(rate.getH1Seminorm(), Real(K) - 0.25);
      EXPECT_LT(rate.getH1Seminorm(), Real(K) + 0.4);
    }
  }

  template <class ContextType>
  void checkNegativeControl(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 5);
    const auto error = solve<1>(mesh, ConductivityData::Field::Affine, true);
    EXPECT_GT(error.getL2(), 1e-3);
    EXPECT_GT(error.getH1Seminorm(), 1e-2);
  }

  template <class ContextType>
  void checkSensitivity(Polytope::Type geometry)
  {
    auto mesh = makeMesh<ContextType>(geometry, 5);
    const auto baseline = solve<2>(mesh, ConductivityData::Field::Smooth);
    const auto refined =
      solve<2>(mesh, ConductivityData::Field::Smooth, false, 14, 1e-14);
    ASSERT_GT(baseline.getL2(), 0);
    ASSERT_GT(baseline.getH1Seminorm(), 0);
    EXPECT_LT(std::abs(refined.getL2() / baseline.getL2() - 1), 1e-6);
    EXPECT_LT(std::abs(refined.getH1Seminorm() / baseline.getH1Seminorm() - 1), 1e-6);
  }

  class PETScConductivityLocalTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(PETScConductivityLocalTest, AffineP1Patch)
  {
    checkPatch<Context::Local, 1>(GetParam());
  }
  TEST_P(PETScConductivityLocalTest, QuadraticP2Patch)
  {
    checkPatch<Context::Local, 2>(GetParam());
  }
  TEST_P(PETScConductivityLocalTest, P1OptimalRates)
  {
    checkRates<Context::Local, 1>(GetParam());
  }
  TEST_P(PETScConductivityLocalTest, P2OptimalRates)
  {
    checkRates<Context::Local, 2>(GetParam());
  }
  TEST_P(PETScConductivityLocalTest, AffinePatchRejectsOmittedCoefficient)
  {
    checkNegativeControl<Context::Local>(GetParam());
  }
  TEST_P(PETScConductivityLocalTest, QuadratureAndSolverSensitivity)
  {
    checkSensitivity<Context::Local>(GetParam());
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScConductivityLocalTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

#ifdef RODIN_USE_MPI
  class PETScConductivityMPITest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(PETScConductivityMPITest, AffineP1Patch)
  {
    checkPatch<Context::MPI, 1>(GetParam());
  }
  TEST_P(PETScConductivityMPITest, QuadraticP2Patch)
  {
    checkPatch<Context::MPI, 2>(GetParam());
  }
  TEST_P(PETScConductivityMPITest, P1OptimalRates)
  {
    checkRates<Context::MPI, 1>(GetParam());
  }
  TEST_P(PETScConductivityMPITest, P2OptimalRates)
  {
    checkRates<Context::MPI, 2>(GetParam());
  }
  TEST_P(PETScConductivityMPITest, AffinePatchRejectsOmittedCoefficient)
  {
    checkNegativeControl<Context::MPI>(GetParam());
  }
  TEST_P(PETScConductivityMPITest, QuadratureAndSolverSensitivity)
  {
    checkSensitivity<Context::MPI>(GetParam());
  }

  TEST_P(PETScConductivityMPITest, GlobalNormCountsOwnedCellsOnce)
  {
    auto mesh = makeMesh<Context::MPI>(GetParam(), 5);
    H1<1, Real, Mesh<Context::MPI>> space(std::integral_constant<size_t, 1>{}, mesh);
    PETSc::Variational::TrialFunction u(space);
    u.getSolution() = RealFunction(0);
    const RealFunction exact(1);
    const VectorFunction gradient(mesh.getDimension(), [](const Point& p) {
      return Math::SpatialVector<Real>::Zero(p.getCoordinates().size());
    });
    const auto error = ErrorNorm::compute(mesh, u.getSolution(), exact, gradient, 12);
    EXPECT_NEAR(error.getL2(), 1, 1e-12);
    EXPECT_NEAR(ErrorNorm::computeL2(mesh, u.getSolution(), exact, 12), 1, 1e-12);
    EXPECT_EQ(error.getH1Seminorm(), 0);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScConductivityMPITest,
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
    Rodin::Tests::Convergence::H::PETScConductivity::environment = &env;
    Rodin::Tests::Convergence::H::PETScConductivity::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
