/*
 * Unit tests for DenseProblem.
 *
 * Tests cover:
 *   - Construction with TrialFunction/TestFunction
 *   - Assembly of dense linear system
 *   - Dense matrix dimensions match FE space size
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Assembly.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_DenseProblem, Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 Vh(mesh);

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    DenseProblem problem(u, v);

    // Should construct without error
    SUCCEED();
  }

  TEST(Rodin_Variational_DenseProblem, Assemble)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 Vh(mesh);

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    DenseProblem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            - Integral(RealFunction(1.0), v)
            + DirichletBC(u, RealFunction(0.0));

    problem.assemble();

    // Dense stiffness matrix should have correct dimensions
    const auto& ls = problem.getLinearSystem();
    const auto& A = ls.getOperator();
    EXPECT_EQ(static_cast<size_t>(A.rows()), Vh.getSize());
    EXPECT_EQ(static_cast<size_t>(A.cols()), Vh.getSize());

    // Dense RHS vector should have correct dimensions
    const auto& b = ls.getVector();
    EXPECT_EQ(static_cast<size_t>(b.size()), Vh.getSize());
  }

  TEST(Rodin_Variational_DenseProblem, StiffnessMatrixIsSymmetric)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 Vh(mesh);

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    DenseProblem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            + DirichletBC(u, RealFunction(0.0));

    problem.assemble();

    const auto& A = problem.getLinearSystem().getOperator();

    // Stiffness matrix for Laplacian should be symmetric
    for (Eigen::Index i = 0; i < A.rows(); i++)
      for (Eigen::Index j = 0; j < A.cols(); j++)
        EXPECT_NEAR(A(i, j), A(j, i), 1e-12);
  }

  TEST(Rodin_Variational_DenseProblem, RHSIsNonZeroWithForcing)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 Vh(mesh);

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    DenseProblem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            - Integral(RealFunction(1.0), v)
            + DirichletBC(u, RealFunction(0.0));

    problem.assemble();

    const auto& b = problem.getLinearSystem().getVector();

    // With f=1, the RHS should have at least some non-zero entries
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
