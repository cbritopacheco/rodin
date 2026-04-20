/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for preassembled BilinearForm/LinearForm block offset assembly
 * in mixed (multi-variable) Problems.
 *
 * The fix under test ensures that when a preassembled BilinearForm or
 * LinearForm is added to a multi-variable Problem, its entries are placed
 * at the correct block offsets in the global matrix and vector.
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Assembly.h>
#include <Rodin/Solver/BiCGSTAB.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  // -----------------------------------------------------------------------
  // Helper: build two mixed systems — one with integrals only, one that
  // replaces one integral term with a preassembled BilinearForm — and
  // verify the resulting global stiffness matrices are identical.
  // -----------------------------------------------------------------------

  /**
   * @brief Mixed P0/P1 problem: compare integral-only vs preassembled BF.
   *
   * System:
   *   Integral(u, v)  – Integral(p, v)  = 0       (test v in P1)
   *   Integral(p, q)                     = (f, q)  (test q in P0)
   *
   * We replace the Integral(p, q) block with a preassembled BilinearForm
   * and verify the global matrix matches.
   */
  TEST(MixedProblem_PreassembledBF, SparseP0P1_MatrixMatchesIntegralOnly)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P0 p0h(mesh);
    P1 p1h(mesh);

    // --- Approach 1: integral-only ---
    TrialFunction u1(p1h);
    TrialFunction p1(p0h);
    TestFunction  v1(p1h);
    TestFunction  q1(p0h);

    Problem integralOnly(u1, v1, p1, q1);
    integralOnly = Integral(u1, v1)
                 - Integral(p1, v1)
                 + Integral(p1, q1)
                 - Integral(RealFunction(1.0), q1);
    integralOnly.assemble();

    // --- Approach 2: replace Integral(p, q) with preassembled BF ---
    TrialFunction u2(p1h);
    TrialFunction p2(p0h);
    TestFunction  v2(p1h);
    TestFunction  q2(p0h);

    // Preassemble the mass matrix on P0
    BilinearForm massP0(p2, q2);
    massP0 = Integral(p2, q2);
    massP0.assemble();

    Problem preassembled(u2, v2, p2, q2);
    preassembled = Integral(u2, v2)
                 - Integral(p2, v2)
                 + massP0
                 - Integral(RealFunction(1.0), q2);
    preassembled.assemble();

    // Compare global matrices
    const auto& A1 = integralOnly.getLinearSystem().getOperator();
    const auto& A2 = preassembled.getLinearSystem().getOperator();

    ASSERT_EQ(A1.rows(), A2.rows());
    ASSERT_EQ(A1.cols(), A2.cols());

    const auto diff = (A1 - A2).norm();
    EXPECT_NEAR(diff, 0.0, 1e-12)
      << "Preassembled BF matrix differs from integral-only matrix by " << diff;
  }

  /**
   * @brief Same test for the RHS vector: compare integral-only vs
   * preassembled LinearForm.
   */
  TEST(MixedProblem_PreassembledLF, SparseP0P1_VectorMatchesIntegralOnly)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P0 p0h(mesh);
    P1 p1h(mesh);

    // --- Approach 1: integral-only ---
    TrialFunction u1(p1h);
    TrialFunction p1(p0h);
    TestFunction  v1(p1h);
    TestFunction  q1(p0h);

    Problem integralOnly(u1, v1, p1, q1);
    integralOnly = Integral(u1, v1)
                 - Integral(p1, v1)
                 + Integral(p1, q1)
                 - Integral(RealFunction(1.0), q1);
    integralOnly.assemble();

    // --- Approach 2: preassemble the load integral on P0 ---
    TrialFunction u2(p1h);
    TrialFunction p2(p0h);
    TestFunction  v2(p1h);
    TestFunction  q2(p0h);

    LinearForm loadP0(q2);
    loadP0 = Integral(RealFunction(1.0), q2);
    loadP0.assemble();

    Problem preassembled(u2, v2, p2, q2);
    preassembled = Integral(u2, v2)
                 - Integral(p2, v2)
                 + Integral(p2, q2)
                 - loadP0;
    preassembled.assemble();

    // Compare global RHS vectors
    const auto& b1 = integralOnly.getLinearSystem().getVector();
    const auto& b2 = preassembled.getLinearSystem().getVector();

    ASSERT_EQ(b1.size(), b2.size());

    const auto diff = (b1 - b2).norm();
    EXPECT_NEAR(diff, 0.0, 1e-12)
      << "Preassembled LF vector differs from integral-only vector by " << diff;
  }

  /**
   * @brief Verify that a solve using preassembled BF in a mixed Problem
   * produces the same numerical solution as the integral-only version.
   */
  TEST(MixedProblem_PreassembledBF, SolveP0P1_SolutionMatches)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 6, 6 });
    mesh.getConnectivity().compute(1, 2);

    P0 p0h(mesh);
    P1 p1h(mesh);

    auto exact = RealFunction(2.0);

    // --- Integral-only ---
    {
      TrialFunction u(p1h);
      TrialFunction p(p0h);
      TestFunction  v(p1h);
      TestFunction  q(p0h);

      Problem mixed(u, v, p, q);
      mixed = Integral(u, v)
            - Integral(p, v)
            + Integral(p, q)
            - Integral(exact, q)
            + DirichletBC(u, exact);

      Solver::BiCGSTAB(mixed).solve();

      GridFunction diff_u(p1h);
      diff_u = Pow(u.getSolution() - exact, 2);
      EXPECT_NEAR(Integral(diff_u).compute(), 0, 1e-12);
    }

    // --- With preassembled P0 mass matrix ---
    {
      TrialFunction u(p1h);
      TrialFunction p(p0h);
      TestFunction  v(p1h);
      TestFunction  q(p0h);

      BilinearForm massP0(p, q);
      massP0 = Integral(p, q);
      massP0.assemble();

      Problem mixed(u, v, p, q);
      mixed = Integral(u, v)
            - Integral(p, v)
            + massP0
            - Integral(exact, q)
            + DirichletBC(u, exact);

      Solver::BiCGSTAB(mixed).solve();

      GridFunction diff_u(p1h);
      diff_u = Pow(u.getSolution() - exact, 2);
      EXPECT_NEAR(Integral(diff_u).compute(), 0, 1e-12);
    }
  }

  /**
   * @brief Verify that a solve using preassembled LF in a mixed Problem
   * produces the same numerical solution as the integral-only version.
   */
  TEST(MixedProblem_PreassembledLF, SolveP0P1_SolutionMatches)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 6, 6 });
    mesh.getConnectivity().compute(1, 2);

    P0 p0h(mesh);
    P1 p1h(mesh);

    auto exact = RealFunction(2.0);

    // --- With preassembled P0 load vector ---
    {
      TrialFunction u(p1h);
      TrialFunction p(p0h);
      TestFunction  v(p1h);
      TestFunction  q(p0h);

      LinearForm loadP0(q);
      loadP0 = Integral(exact, q);
      loadP0.assemble();

      Problem mixed(u, v, p, q);
      mixed = Integral(u, v)
            - Integral(p, v)
            + Integral(p, q)
            - loadP0
            + DirichletBC(u, exact);

      Solver::BiCGSTAB(mixed).solve();

      GridFunction diff_u(p1h);
      diff_u = Pow(u.getSolution() - exact, 2);
      EXPECT_NEAR(Integral(diff_u).compute(), 0, 1e-12);
    }
  }
}
#include <gtest/gtest.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Assembly.h>
#include <iostream>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Debug, PrintMatrixDiff) {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);
    
    P0 p0h(mesh);
    P1 p1h(mesh);
    
    std::cout << "P1 size: " << p1h.getSize() << ", P0 size: " << p0h.getSize() << "\n";

    // Integral-only
    TrialFunction u1(p1h); TrialFunction p1(p0h);
    TestFunction  v1(p1h); TestFunction  q1(p0h);
    Problem integralOnly(u1, v1, p1, q1);
    integralOnly = Integral(u1, v1) - Integral(p1, v1) + Integral(p1, q1) - Integral(RealFunction(1.0), q1);
    integralOnly.assemble();
    const auto& A1 = integralOnly.getLinearSystem().getOperator();
    const auto& b1 = integralOnly.getLinearSystem().getVector();

    // Preassembled BF  
    TrialFunction u2(p1h); TrialFunction p2(p0h);
    TestFunction  v2(p1h); TestFunction  q2(p0h);
    BilinearForm massP0(p2, q2);
    massP0 = Integral(p2, q2);
    massP0.assemble();
    std::cout << "massP0: " << massP0.getOperator().rows() << "x" << massP0.getOperator().cols() 
              << ", nnz=" << massP0.getOperator().nonZeros() << "\n";
    double msum = 0;
    for (int k = 0; k < massP0.getOperator().outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(massP0.getOperator(), k); it; ++it)
            msum += it.value();
    std::cout << "massP0 sum of entries: " << msum << "\n";

    Problem preassembled(u2, v2, p2, q2);
    preassembled = Integral(u2, v2) - Integral(p2, v2) + massP0 - Integral(RealFunction(1.0), q2);
    preassembled.assemble();
    const auto& A2 = preassembled.getLinearSystem().getOperator();
    const auto& b2 = preassembled.getLinearSystem().getVector();

    std::cout << "A1: " << A1.rows() << "x" << A1.cols() << ", nnz=" << A1.nonZeros() << "\n";
    std::cout << "A2: " << A2.rows() << "x" << A2.cols() << ", nnz=" << A2.nonZeros() << "\n";
    std::cout << "b1 norm: " << b1.norm() << "\n";
    std::cout << "b2 norm: " << b2.norm() << "\n";

    auto diffA = (A1 - A2).eval();
    std::cout << "diff norm: " << diffA.norm() << "\n";
    std::cout << "Diff non-zero entries:\n";
    for (int k = 0; k < diffA.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(diffA, k); it; ++it)
            if (std::abs(it.value()) > 1e-14)
                std::cout << "  [" << it.row() << "," << it.col() << "]=" << it.value() << "\n";

    auto diffb = (b1 - b2).eval();
    std::cout << "b diff:\n";
    for (int i = 0; i < diffb.size(); ++i)
        if (std::abs(diffb[i]) > 1e-14)
            std::cout << "  [" << i << "]=" << diffb[i] << "\n";
}
