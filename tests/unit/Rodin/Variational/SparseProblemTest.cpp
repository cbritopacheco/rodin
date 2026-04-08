/*
 * Unit tests for SparseProblem.
 *
 * Tests cover:
 *   - Construction with TrialFunction/TestFunction
 *   - Assembly of sparse linear system
 *   - Sparse matrix dimensions match FE space size
 *   - Sparsity (non-zero count much less than dense)
 *   - Solve and verify non-trivial solution
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_SparseProblem, Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 Vh(mesh);

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    SparseProblem problem(u, v);

    SUCCEED();
  }

  TEST(Rodin_Variational_SparseProblem, Assemble)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 Vh(mesh);

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    SparseProblem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            - Integral(RealFunction(1.0), v)
            + DirichletBC(u, RealFunction(0.0));

    problem.assemble();

    const auto& ls = problem.getLinearSystem();
    const auto& A = ls.getOperator();
    EXPECT_EQ(static_cast<size_t>(A.rows()), Vh.getSize());
    EXPECT_EQ(static_cast<size_t>(A.cols()), Vh.getSize());

    const auto& b = ls.getVector();
    EXPECT_EQ(static_cast<size_t>(b.size()), Vh.getSize());
  }

  TEST(Rodin_Variational_SparseProblem, MatrixIsSparse)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P1 Vh(mesh);

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    SparseProblem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            + DirichletBC(u, RealFunction(0.0));

    problem.assemble();

    const auto& A = problem.getLinearSystem().getOperator();
    const size_t n = Vh.getSize();
    const size_t nnz = static_cast<size_t>(A.nonZeros());

    // For a sparse FEM matrix, the number of non-zeros should be
    // much less than n*n (a dense matrix)
    EXPECT_LT(nnz, n * n);
    // But should have at least n non-zeros (diagonal)
    EXPECT_GE(nnz, n);
  }

  TEST(Rodin_Variational_SparseProblem, StiffnessMatrixIsSymmetric)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 Vh(mesh);

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    SparseProblem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            + DirichletBC(u, RealFunction(0.0));

    problem.assemble();

    const auto& A = problem.getLinearSystem().getOperator();
    const auto AT = Eigen::SparseMatrix<double>(A.transpose());

    // Check symmetry via norm of difference
    const double diff = (A - AT).norm();
    EXPECT_NEAR(diff, 0.0, 1e-12);
  }

  TEST(Rodin_Variational_SparseProblem, SolvePoisson)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);
    P1 Vh(mesh);

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    SparseProblem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            - Integral(RealFunction(1.0), v)
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

  TEST(Rodin_Variational_SparseProblem, RHSIsNonZeroWithForcing)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 Vh(mesh);

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    SparseProblem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            - Integral(RealFunction(1.0), v)
            + DirichletBC(u, RealFunction(0.0));

    problem.assemble();

    const auto& b = problem.getLinearSystem().getVector();

    bool hasNonZero = false;
    for (Eigen::Index i = 0; i < b.size(); i++)
    {
      if (std::abs(b(i)) > 1e-14)
      {
        hasNonZero = true;
        break;
      }
    }
    EXPECT_TRUE(hasNonZero);
  }
}
