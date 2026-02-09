/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <algorithm>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"
#include "Rodin/Solver/BiCGSTAB.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::P0P1Mixed
{
  template <size_t NX, size_t NY, size_t NZ = 1>
  class Manufactured_P0P1_Mixed_Test : public ::testing::TestWithParam<Polytope::Type>
  {
  protected:
    void SetUp() override
    {
      const auto geom = GetParam();
      if (geom == Polytope::Type::Tetrahedron || geom == Polytope::Type::Hexahedron || geom == Polytope::Type::Wedge)
      {
        m_mesh = Mesh().UniformGrid(geom, { NX, NY, NZ });
        m_mesh.scale(1.0 / (NX - 1));
        m_mesh.getConnectivity().compute(2, 3);
        m_mesh.getConnectivity().compute(3, 2);
      }
      else
      {
        m_mesh = Mesh().UniformGrid(geom, { NX, NY });
        m_mesh.scale(1.0 / (NX - 1));
        m_mesh.getConnectivity().compute(1, 2);
        m_mesh.getConnectivity().compute(2, 1);
      }
      m_mesh.getConnectivity().compute(1, 0);
    }

    const Mesh<Context::Local>& getMesh() const { return m_mesh; }

    static auto rhs_polynomial()
    {
      return RealFunction([](const Geometry::Point& p) -> double
      {
        return 2.0 * (p.x() * (1.0 - p.x()) + p.y() * p.y());
      });
    }

    static auto rhs_sine()
    {
      auto pi = Rodin::Math::Constants::pi();
      return RealFunction([pi](const Geometry::Point& p) -> double
      {
        return std::sin(pi * p.x()) * std::sin(pi * p.y()) * std::sin(pi * p.z());
      });
    }

  private:
    Mesh<Context::Local> m_mesh;
  };

  using Manufactured_P0P1_Mixed_Test_10x10 =
    Manufactured_P0P1_Mixed_Test<10, 10>;
  using Manufactured_P0P1_Mixed_Test_6x6x6 =
    Manufactured_P0P1_Mixed_Test<6, 6, 6>;

  TEST_P(Manufactured_P0P1_Mixed_Test_10x10, P0P1_Mixed_ConstantExactSolutionResidual)
  {
    const auto& mesh = getMesh();

    P0 p0h(mesh);
    P1 p1h(mesh);

    auto exact_solution = RealFunction(2.0); // arbitrary non-zero constant exact solution in both spaces
    // Mixed saddle-point system weakly enforces u ≈ p and p ≈ exact_solution; picking
    // exact_solution as the forcing makes (u, p) = (exact_solution, exact_solution)
    // the exact discrete solution when combined with the Dirichlet conditions below.

    TrialFunction u(p1h);
    TrialFunction p(p0h);
    TestFunction  v(p1h);
    TestFunction  q(p0h);

    Problem mixed(u, v, p, q);
    mixed = Integral(u, v)
          - Integral(p, v)
          + Integral(p, q)
          - Integral(exact_solution, q)
          + DirichletBC(u, exact_solution);

    BiCGSTAB(mixed).solve();

    GridFunction u_exact_coeffs(p1h); // grid-function coefficients of the exact solution for u
    u_exact_coeffs = exact_solution;
    GridFunction p_exact_coeffs(p0h); // grid-function coefficients of the exact solution for p
    p_exact_coeffs = exact_solution;

    auto& A = mixed.getLinearSystem().getOperator();
    auto& b = mixed.getLinearSystem().getVector();
    auto& x = mixed.getLinearSystem().getSolution();

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

    GridFunction diff_u(p1h);
    diff_u = Pow(u.getSolution() - exact_solution, 2);
    EXPECT_NEAR(Integral(diff_u).compute(), 0, 1e-12);

    GridFunction diff_p(p0h);
    diff_p = Pow(p.getSolution() - exact_solution, 2);
    EXPECT_NEAR(Integral(diff_p).compute(), 0, 1e-12);
  }

  TEST_P(Manufactured_P0P1_Mixed_Test_10x10, P0P1_MixedProblem_PolynomialRHS)
  {
    const auto& mesh = getMesh();

    P0 p0h(mesh);
    P1 p1h(mesh);

    TrialFunction u0(p1h);
    TrialFunction u(p1h);
    TestFunction  v(p1h);

    TrialFunction p(p0h);
    TestFunction  q(p0h);

    const auto f = rhs_polynomial();

    // First, L2 projection of f onto P0: solve (p, q) = (f, q)
    Problem p_l2(p, q);
    p_l2 = Integral(p, q) - Integral(f, q);
    BiCGSTAB(p_l2).solve();

    // L2 projection of p onto P1: solve (u, v) = (p, v)
    Problem u_l2(u0, v);
    u_l2 = Integral(u0, v) - Integral(p.getSolution(), v);
    BiCGSTAB(u_l2).solve();

    // Reset for mixed solve
    p.getSolution() = 0.0;

    // Mixed system:
    // (u, v) - (p, v) + (p, q) - (f, q) = 0
    Problem mixed(u, v, p, q);
    mixed = Integral(u, v)
          - Integral(p, v)
          + Integral(p, q)
          - Integral(f, q);

    BiCGSTAB(mixed).solve();

    // Check that u matches the L2 projection of f onto P1 (from u_l2)
    GridFunction diff(p1h);
    diff = Pow(u.getSolution() - u0.getSolution(), 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  TEST_P(Manufactured_P0P1_Mixed_Test_6x6x6, P0P1_MixedProblem_SineRHS)
  {
    const auto& mesh = getMesh();

    P0 p0h(mesh);
    P1 p1h(mesh);

    TrialFunction u0(p1h);
    TrialFunction u(p1h);
    TestFunction  v(p1h);

    TrialFunction p(p0h);
    TestFunction  q(p0h);

    const auto f = rhs_sine();

    Problem p_l2(p, q);
    p_l2 = Integral(p, q) - Integral(f, q);
    BiCGSTAB(p_l2).solve();

    Problem u_l2(u0, v);
    u_l2 = Integral(u0, v) - Integral(p.getSolution(), v);
    BiCGSTAB(u_l2).solve();

    p.getSolution() = 0.0;
    u.getSolution() = 0.0;

    Problem mixed(u, v, p, q);
    mixed = Integral(u, v)
          - Integral(p, v)
          + Integral(p, q)
          - Integral(f, q);

    BiCGSTAB(mixed).solve();

    GridFunction diff(p1h);
    diff = Pow(u.getSolution() - u0.getSolution(), 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage2D,
    Manufactured_P0P1_Mixed_Test_10x10,
    ::testing::Values(
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral)
  );

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Manufactured_P0P1_Mixed_Test_6x6x6,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge)
  );
}
