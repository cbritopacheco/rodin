/*
 * @file ProblemTest.cpp
 * @brief Unit tests for Problem, SparseProblem, and DenseProblem classes.
 */
#include <gtest/gtest.h>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Variational_Problem, ConstructionFromTrialTest)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  P1 Vh(mesh);
  TrialFunction u(Vh);
  TestFunction v(Vh);

  Problem problem(u, v);
  // Construction succeeds
  SUCCEED();
}

TEST(Rodin_Variational_Problem, AssemblePoisson)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  mesh.getConnectivity().compute(1, 2);
  P1 Vh(mesh);
  TrialFunction u(Vh);
  TestFunction v(Vh);

  RealFunction f(1.0);

  Problem problem(u, v);
  problem = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, RealFunction(0.0));
  problem.assemble();

  // Assemble succeeds
  SUCCEED();
}

TEST(Rodin_Variational_Problem, SolvePoisson)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  mesh.getConnectivity().compute(1, 2);
  P1 Vh(mesh);
  TrialFunction u(Vh);
  TestFunction v(Vh);

  RealFunction f(1.0);

  Problem problem(u, v);
  problem = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, RealFunction(0.0));

  Solver::CG(problem).solve();

  // Solution should be non-trivial
  const auto& sol = u.getSolution();
  bool hasNonZero = false;
  for (Index i = 0; i < sol.getSize(); ++i)
  {
    if (std::abs(sol.getData()(i)) > 1e-15)
    {
      hasNonZero = true;
      break;
    }
  }
  EXPECT_TRUE(hasNonZero);
}

TEST(Rodin_Variational_Problem, SolutionSize)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  mesh.getConnectivity().compute(1, 2);
  P1 Vh(mesh);
  TrialFunction u(Vh);
  TestFunction v(Vh);

  RealFunction f(1.0);

  Problem problem(u, v);
  problem = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, RealFunction(0.0));

  Solver::CG(problem).solve();

  const auto& sol = u.getSolution();
  EXPECT_EQ(sol.getSize(), Vh.getSize());
}
