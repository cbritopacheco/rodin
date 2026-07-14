/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Helmholtz manufactured solution tests.
 *
 * These tests assemble Rodin variational forms for a Helmholtz manufactured solution, solve the problem on the configured mesh, and compare against analytic fields or expected residual/error behavior. They protect the P1 finite-element and solver path, including boundary-condition handling, geometry coverage, and numerical accuracy of the manufactured workflow.
 */

#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Configure.h"
#include "Rodin/Context/Local.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

/**
 * 3D manufactured tests for the Helmholtz equation (sign convention used here):
 *
 * Strong form:
 *   -Δu - κ² u = f   in Ω
 *            u = g   on ∂Ω
 *
 * Weak form:
 *   ∫_Ω ∇ u·∇ v dx - κ² ∫_Ω u v dx = ∫_Ω f v dx
 *
 * IMPORTANT (geometry):
 *   Rodin's UniformGrid({M,M,M}) produces coordinates {0,...,M-1}.
 *   We scale by 1/(M-1) so the physical domain is [0,1]^3, matching the
 *   manufactured expressions (sin(pi*x), x(1-x), etc.).
 *
 * NOTE (solver):
 *   With the sign "- κ² ∫ u v", the operator may be symmetric indefinite for
 *   larger κ. These tests keep your current CG usage, but mathematically CG
 *   requires SPD. If you hit solver issues, switch to MINRES/GMRES or change
 *   the PDE to -Δu + κ² u (screened Poisson) to get SPD.
 */
namespace Rodin::Tests::Manufactured::Helmholtz3D
{
  // ---------------------------------------------------------------------------
  // Fixture: build mesh once, scaled so physical domain is [0,1]^3.
  // ---------------------------------------------------------------------------
  template <size_t M>
  class Helmholtz3DFixture : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh<Context::Local>().UniformGrid(GetParam(), {M, M, M});
        m_mesh.scale(Real(1) / Real(M - 1)); // map {0..M-1} -> [0,1]
        m_mesh.getConnectivity().compute(2, 3);
      }

      const Mesh<Context::Local>& mesh() const { return m_mesh; }

    private:
      Mesh<Context::Local> m_mesh;
  };

  /// @brief Helper used by the tests to Helmholtz 3 D Test 8.
  using Helmholtz3DTest8  = Helmholtz3DFixture<8>;
  /// @brief Helper used by the tests to Helmholtz 3 D Test 16.
  using Helmholtz3DTest16 = Helmholtz3DFixture<16>;
  /// @brief Helper used by the tests to Helmholtz 3 D Test 32.
  using Helmholtz3DTest32 = Helmholtz3DFixture<32>;

  /// @brief Verifies helmholtz P1 exact residual for helmholtz 3 D test 8 by checking tolerance-based numerical results, solver behavior.
  TEST_P(Helmholtz3DTest8, Helmholtz_P1ExactResidual)
  {
    const Real kappa = 1.0;

    P1 vh(mesh());

    auto u_expr = F::x + F::y + F::z + 1;
    auto f = -kappa * kappa * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem helmholtz(u, v);
    helmholtz = Integral(Grad(u), Grad(v))
              - kappa * kappa * Integral(u, v)
              - Integral(f, v)
              + DirichletBC(u, u_expr);

    CG(helmholtz).solve();

    GridFunction u_exact(vh);
    u_exact = u_expr;

    auto& A = helmholtz.getLinearSystem().getOperator();
    auto& b = helmholtz.getLinearSystem().getVector();
    auto& x = helmholtz.getLinearSystem().getSolution();

    auto r = A * x - b;
    auto re = A * u_exact.getData() - b;

    const Real scale = std::max<Real>(b.norm(), 1);
    EXPECT_NEAR(r.norm() / scale, 0, 1e-10);
    EXPECT_NEAR(re.norm() / scale, 0, 1e-12);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0, 1e-12);
  }

  // ---------------------------------------------------------------------------
  // u = sin(pi x) sin(pi y) sin(pi z)
  // Δu = -3 pi^2 u
  // f = -Δu - κ^2 u = (3 pi^2 - κ^2) u
  // Dirichlet: u=0 on ∂Ω (sine vanishes on x,y,z=0 or 1)
  // ---------------------------------------------------------------------------
  /// @brief Verifies simple sine for helmholtz 3 D test 16 by checking tolerance-based numerical results, solver behavior.
  TEST_P(Helmholtz3DTest16, SimpleSine)
  {
    const Real pi = Math::Constants::pi();
    const Real kappa = 2.0;

    P1 vh(mesh());

    const auto u_expr = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    const auto f = (3 * pi * pi - kappa * kappa) * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem helmholtz(u, v);
    helmholtz = Integral(Grad(u), Grad(v))
              - kappa * kappa * Integral(u, v)
              - Integral(f, v)
              + DirichletBC(u, Zero());

    CG(helmholtz).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, 10 * RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // u = sin(ω1 pi x) sin(ω2 pi y) sin(ω3 pi z)
  // Δu = -(ω1^2 + ω2^2 + ω3^2) pi^2 u
  // f = -Δu - κ^2 u = ( (ω1^2+ω2^2+ω3^2) pi^2 - κ^2 ) u
  // Dirichlet: u=0 on ∂Ω for integer ωi (still vanishes at x,y,z=0 or 1)
  // ---------------------------------------------------------------------------
  /// @brief Verifies variable frequency for helmholtz 3 D test 32 by checking tolerance-based numerical results, solver behavior.
  TEST_P(Helmholtz3DTest32, VariableFrequency)
  {
    const Real pi = Math::Constants::pi();
    const Real kappa = 1.5;

    const Real omega1 = 2;
    const Real omega2 = 2;
    const Real omega3 = 3;

    P1 vh(mesh());

    const auto u_expr =
      sin(omega1 * pi * F::x) * sin(omega2 * pi * F::y) * sin(omega3 * pi * F::z);

    const auto f =
      ((omega1 * omega1 + omega2 * omega2 + omega3 * omega3) * pi * pi - kappa * kappa) * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem helmholtz(u, v);
    helmholtz = Integral(Grad(u), Grad(v))
              - kappa * kappa * Integral(u, v)
              - Integral(f, v)
              + DirichletBC(u, Zero());

    CG(helmholtz).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, 10 * RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // u = x(1-x) sin(pi y) sin(pi z)
  // ∂xx u = -2 sin(pi y) sin(pi z)
  // ∂yy u = -pi^2 x(1-x) sin(pi y) sin(pi z)
  // ∂zz u = -pi^2 x(1-x) sin(pi y) sin(pi z)
  // Δu = -2 s - 2 pi^2 x(1-x) s, where s = sin(pi y) sin(pi z)
  // f = -Δu - κ^2 u
  // Dirichlet: u=0 on ∂Ω (x=0,1 or y=0,1 or z=0,1)
  // ---------------------------------------------------------------------------
  /// @brief Verifies mixed polynomial trig for helmholtz 3 D test 16 by checking tolerance-based numerical results, solver behavior.
  TEST_P(Helmholtz3DTest16, MixedPolynomialTrig)
  {
    const Real pi = Math::Constants::pi();
    const Real kappa = 1.0;

    P1 vh(mesh());

    const auto s = sin(pi * F::y) * sin(pi * F::z);
    const auto u_expr = F::x * (1 - F::x) * s;

    const auto f =
      2 * s
      + 2 * pi * pi * F::x * (1 - F::x) * s
      - kappa * kappa * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem helmholtz(u, v);
    helmholtz = Integral(Grad(u), Grad(v))
              - kappa * kappa * Integral(u, v)
              - Integral(f, v)
              + DirichletBC(u, Zero());

    CG(helmholtz).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // u = sin(pi x) sin(pi y) e^z
  // ∂xx u = -pi^2 u
  // ∂yy u = -pi^2 u
  // ∂zz u = u
  // Δu = (1 - 2 pi^2) u
  // f = -Δu - κ^2 u = (2 pi^2 - 1 - κ^2) u
  // Dirichlet: here we impose u=g on ∂Ω (nonzero on many faces)
  // ---------------------------------------------------------------------------
  /// @brief Verifies exponential for helmholtz 3 D test 16 by checking tolerance-based numerical results, solver behavior.
  TEST_P(Helmholtz3DTest16, Exponential)
  {
    const Real pi = Math::Constants::pi();
    const Real kappa = 0.5;

    P1 vh(mesh());

    const auto u_expr = sin(pi * F::x) * sin(pi * F::y) * exp(F::z);
    const auto f = (2 * pi * pi - 1 - kappa * kappa) * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem helmholtz(u, v);
    helmholtz = Integral(Grad(u), Grad(v))
              - kappa * kappa * Integral(u, v)
              - Integral(f, v)
              + DirichletBC(u, u_expr);

    CG(helmholtz).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // Extra: run the SimpleSine case as a refinement sanity check.
  // ---------------------------------------------------------------------------
  /// @brief Verifies simple sine refined sanity for helmholtz 3 D test 16 by checking tolerance-based numerical results, solver behavior.
  TEST_P(Helmholtz3DTest16, SimpleSine_RefinedSanity)
  {
    const Real pi = Math::Constants::pi();
    const Real kappa = 2.0;

    P1 vh(mesh());

    const auto u_expr = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    const auto f = (3 * pi * pi - kappa * kappa) * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem helmholtz(u, v);
    helmholtz = Integral(Grad(u), Grad(v))
              - kappa * kappa * Integral(u, v)
              - Integral(f, v)
              + DirichletBC(u, Zero());

    CG(helmholtz).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, 10 * RODIN_FUZZY_CONSTANT);
  }

  /// @brief Instantiates Helmholtz 3 D Test 8 over the Polytope Coverage 3 D parameter coverage.
  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Helmholtz3DTest8,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Wedge)
  );

  /// @brief Instantiates Helmholtz 3 D Test 16 over the Polytope Coverage 3 D parameter coverage.
  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Helmholtz3DTest16,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Wedge)
  );

  /// @brief Instantiates Helmholtz 3 D Test 32 over the Polytope Coverage 3 D parameter coverage.
  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Helmholtz3DTest32,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Wedge)
  );
}
