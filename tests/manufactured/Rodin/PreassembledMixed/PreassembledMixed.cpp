/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Manufactured tests for preassembled BilinearForm / LinearForm in mixed
 * (multi-variable) Problems.
 *
 * These tests verify that substituting an integral term with a preassembled
 * BilinearForm or LinearForm in a mixed Problem produces identical numerical
 * results to the integral-only formulation.
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/BiCGSTAB.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::PreassembledMixed
{
  /**
   * @brief Fixture for preassembled mixed tests, parameterised on polytope.
   */
  template <size_t NX, size_t NY>
  class PreassembledMixed_Test : public ::testing::TestWithParam<Polytope::Type>
  {
  protected:
    void SetUp() override
    {
      const auto geom = GetParam();
      m_mesh = Mesh().UniformGrid(geom, { NX, NY });
      m_mesh.scale(1.0 / (NX - 1));
      m_mesh.getConnectivity().compute(1, 2);
      m_mesh.getConnectivity().compute(2, 1);
      m_mesh.getConnectivity().compute(1, 0);
    }

    const Mesh<Context::Local>& getMesh() const { return m_mesh; }

  private:
    Mesh<Context::Local> m_mesh;
  };

  using PreassembledMixed_10x10 = PreassembledMixed_Test<10, 10>;

  // -----------------------------------------------------------------------
  // Test 1: Constant exact solution with preassembled P0 mass BF.
  //
  // Mixed system:
  //   (u, v) – (p, v) + [M_p0](p, q) – (f, q) + BC(u = exact)
  //
  // where [M_p0] is the P0 mass matrix supplied as a preassembled
  // BilinearForm instead of Integral(p, q).
  //
  // With exact_solution = 2, the exact discrete solution is (u, p) =
  // (exact, exact).
  // -----------------------------------------------------------------------
  TEST_P(PreassembledMixed_10x10,
    P0P1_PreassembledBF_ConstantExactSolution)
  {
    const auto& mesh = getMesh();

    P0 p0h(mesh);
    P1 p1h(mesh);

    auto exact_solution = RealFunction(2.0);

    TrialFunction u(p1h);
    TrialFunction p(p0h);
    TestFunction  v(p1h);
    TestFunction  q(p0h);

    // Preassemble the P0 mass matrix
    BilinearForm massP0(p, q);
    massP0 = Integral(p, q);
    massP0.assemble();

    Problem mixed(u, v, p, q);
    mixed = Integral(u, v)
          - Integral(p, v)
          + massP0
          - Integral(exact_solution, q)
          + DirichletBC(u, exact_solution);

    BiCGSTAB(mixed).solve();

    // Check residual
    auto& A = mixed.getLinearSystem().getOperator();
    auto& b = mixed.getLinearSystem().getVector();
    auto& x = mixed.getLinearSystem().getSolution();

    GridFunction u_exact_coeffs(p1h);
    u_exact_coeffs = exact_solution;
    GridFunction p_exact_coeffs(p0h);
    p_exact_coeffs = exact_solution;

    auto x_exact = x;
    const auto uSize = u_exact_coeffs.getData().size();
    const auto pSize = p_exact_coeffs.getData().size();
    x_exact.head(uSize) = u_exact_coeffs.getData();
    x_exact.tail(pSize) = p_exact_coeffs.getData();

    auto r = A * x - b;
    auto re = A * x_exact - b;

    const Real scale = std::max<Real>(b.norm(), 1);
    EXPECT_NEAR(r.norm() / scale, 0, 1e-10);
    EXPECT_NEAR(re.norm() / scale, 0, 1e-12);

    // Check L2 error of u
    GridFunction diff_u(p1h);
    diff_u = Pow(u.getSolution() - exact_solution, 2);
    EXPECT_NEAR(Integral(diff_u).compute(), 0, 1e-12);

    // Check L2 error of p
    GridFunction diff_p(p0h);
    diff_p = Pow(p.getSolution() - exact_solution, 2);
    EXPECT_NEAR(Integral(diff_p).compute(), 0, 1e-12);
  }

  // -----------------------------------------------------------------------
  // Test 2: Constant exact solution with preassembled P0 load LF.
  //
  // Mixed system:
  //   (u, v) – (p, v) + (p, q) – [L_p0](q) + BC(u = exact)
  //
  // where [L_p0] is the load vector \int f\cdotq supplied as a preassembled
  // LinearForm instead of Integral(f, q).
  // -----------------------------------------------------------------------
  TEST_P(PreassembledMixed_10x10,
    P0P1_PreassembledLF_ConstantExactSolution)
  {
    const auto& mesh = getMesh();

    P0 p0h(mesh);
    P1 p1h(mesh);

    auto exact_solution = RealFunction(2.0);

    TrialFunction u(p1h);
    TrialFunction p(p0h);
    TestFunction  v(p1h);
    TestFunction  q(p0h);

    // Preassemble the P0 load vector
    LinearForm loadP0(q);
    loadP0 = Integral(exact_solution, q);
    loadP0.assemble();

    Problem mixed(u, v, p, q);
    mixed = Integral(u, v)
          - Integral(p, v)
          + Integral(p, q)
          - loadP0
          + DirichletBC(u, exact_solution);

    BiCGSTAB(mixed).solve();

    // Check residual
    auto& A = mixed.getLinearSystem().getOperator();
    auto& b = mixed.getLinearSystem().getVector();
    auto& x = mixed.getLinearSystem().getSolution();

    const Real scale = std::max<Real>(b.norm(), 1);
    auto r = A * x - b;
    EXPECT_NEAR(r.norm() / scale, 0, 1e-10);

    // Check L2 errors
    GridFunction diff_u(p1h);
    diff_u = Pow(u.getSolution() - exact_solution, 2);
    EXPECT_NEAR(Integral(diff_u).compute(), 0, 1e-12);

    GridFunction diff_p(p0h);
    diff_p = Pow(p.getSolution() - exact_solution, 2);
    EXPECT_NEAR(Integral(diff_p).compute(), 0, 1e-12);
  }

  // -----------------------------------------------------------------------
  // Test 3: Both BF and LF preassembled simultaneously.
  // -----------------------------------------------------------------------
  TEST_P(PreassembledMixed_10x10,
    P0P1_PreassembledBFAndLF_ConstantExactSolution)
  {
    const auto& mesh = getMesh();

    P0 p0h(mesh);
    P1 p1h(mesh);

    auto exact_solution = RealFunction(2.0);

    TrialFunction u(p1h);
    TrialFunction p(p0h);
    TestFunction  v(p1h);
    TestFunction  q(p0h);

    // Preassemble both
    BilinearForm massP0(p, q);
    massP0 = Integral(p, q);
    massP0.assemble();

    LinearForm loadP0(q);
    loadP0 = Integral(exact_solution, q);
    loadP0.assemble();

    Problem mixed(u, v, p, q);
    mixed = Integral(u, v)
          - Integral(p, v)
          + massP0
          - loadP0
          + DirichletBC(u, exact_solution);

    BiCGSTAB(mixed).solve();

    // Check L2 errors
    GridFunction diff_u(p1h);
    diff_u = Pow(u.getSolution() - exact_solution, 2);
    EXPECT_NEAR(Integral(diff_u).compute(), 0, 1e-12);

    GridFunction diff_p(p0h);
    diff_p = Pow(p.getSolution() - exact_solution, 2);
    EXPECT_NEAR(Integral(diff_p).compute(), 0, 1e-12);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage2D,
    PreassembledMixed_10x10,
    ::testing::Values(
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral)
  );
}
