/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
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
 * 3D manufactured tests for the reaction–diffusion equation:
 *
 * Strong form:
 *   -Δu + α u = f   in Ω
 *             u = g   on ∂Ω
 *
 * Weak form:
 *   ∫_Ω ∇ u·∇ v dx + α ∫_Ω u v dx = ∫_Ω f v dx
 *
 * Geometry note:
 *   UniformGrid({M,M,M}) creates coordinates {0,…,M−1}.
 *   We scale by 1/(M−1) so the physical domain is [0,1]^3,
 *   matching the manufactured expressions.
 */

namespace Rodin::Tests::Manufactured::ReactionDiffusion3D
{
  // ---------------------------------------------------------------------------
  // Fixture: shared mesh with correct physical scaling
  // ---------------------------------------------------------------------------
  template <size_t M>
  class ReactionDiffusion3DFixture : public ::testing::TestWithParam<Polytope::Type>
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

  using ReactionDiffusion3DTest8  = ReactionDiffusion3DFixture<8>;
  using ReactionDiffusion3DTest16 = ReactionDiffusion3DFixture<16>;

  TEST_P(ReactionDiffusion3DTest8, ReactionDiffusion_P1ExactResidual)
  {
    const Real alpha = 1.0;

    P1 vh(mesh());

    auto u_expr = F::x + F::y + F::z + 1;
    auto f = alpha * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem rd(u, v);
    rd = Integral(Grad(u), Grad(v))
       + alpha * Integral(u, v)
       - Integral(f, v)
       + DirichletBC(u, u_expr);

    CG(rd).solve();

    GridFunction u_exact(vh);
    u_exact = u_expr;

    auto& A = rd.getLinearSystem().getOperator();
    auto& b = rd.getLinearSystem().getVector();
    auto& x = rd.getLinearSystem().getSolution();

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
  // f = -Δu + α u = (3 pi^2 + α) u
  // ---------------------------------------------------------------------------
  TEST_P(ReactionDiffusion3DTest16, SimpleSine)
  {
    const Real pi = Math::Constants::pi();
    const Real alpha = 1.0;

    P1 vh(mesh());

    const auto u_expr =
      sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    const auto f = (3 * pi * pi + alpha) * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem rd(u, v);
    rd = Integral(Grad(u), Grad(v))
       + alpha * Integral(u, v)
       - Integral(f, v)
       + DirichletBC(u, Zero());

    CG(rd).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0.0, 2 * RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // Polynomial:
  // u = x(1-x)y(1-y)z(1-z)
  // Δu = 2y(1-y)z(1-z) + 2x(1-x)z(1-z) + 2x(1-x)y(1-y)
  // f = -Δu + α u
  // ---------------------------------------------------------------------------
  TEST_P(ReactionDiffusion3DTest16, Polynomial)
  {
    const Real alpha = 2.0;

    P1 vh(mesh());

    const auto u_expr =
      F::x * (1 - F::x) *
      F::y * (1 - F::y) *
      F::z * (1 - F::z);

    const auto lap_u =
      -2 * F::y * (1 - F::y) * F::z * (1 - F::z)
      -2 * F::x * (1 - F::x) * F::z * (1 - F::z)
      -2 * F::x * (1 - F::x) * F::y * (1 - F::y);

    const auto f = -lap_u + alpha * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem rd(u, v);
    rd = Integral(Grad(u), Grad(v))
       + alpha * Integral(u, v)
       - Integral(f, v)
       + DirichletBC(u, Zero());

    CG(rd).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // Mixed polynomial–trigonometric
  // u = x(1-x) sin(pi y) sin(pi z)
  // ---------------------------------------------------------------------------
  TEST_P(ReactionDiffusion3DTest16, MixedPolynomialTrig)
  {
    const Real pi = Math::Constants::pi();
    const Real alpha = 0.5;

    P1 vh(mesh());

    const auto s = sin(pi * F::y) * sin(pi * F::z);
    const auto u_expr = F::x * (1 - F::x) * s;

    const auto lap_u =
      -2 * s
      -2 * pi * pi * F::x * (1 - F::x) * s;

    const auto f = -lap_u + alpha * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem rd(u, v);
    rd = Integral(Grad(u), Grad(v))
       + alpha * Integral(u, v)
       - Integral(f, v)
       + DirichletBC(u, Zero());

    CG(rd).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // Exponential with non-homogeneous Dirichlet
  // u = cos(pi x) cos(pi y) e^z
  // Δu = (1 - 2 pi^2) u
  // ---------------------------------------------------------------------------
  TEST_P(ReactionDiffusion3DTest16, Exponential)
  {
    const Real pi = Math::Constants::pi();
    const Real alpha = 1.5;

    P1 vh(mesh());

    const auto u_expr =
      cos(pi * F::x) * cos(pi * F::y) * exp(F::z);

    const auto lap_u = (1 - 2 * pi * pi) * u_expr;
    const auto f = -lap_u + alpha * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem rd(u, v);
    rd = Integral(Grad(u), Grad(v))
       + alpha * Integral(u, v)
       - Integral(f, v)
       + DirichletBC(u, u_expr);

    CG(rd).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0.0, 2 * RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // Extra: refinement sanity check
  // ---------------------------------------------------------------------------
  TEST_P(ReactionDiffusion3DTest16, SimpleSine_RefinedSanity)
  {
    const Real pi = Math::Constants::pi();
    const Real alpha = 1.0;

    P1 vh(mesh());

    const auto u_expr =
      sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    const auto f = (3 * pi * pi + alpha) * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem rd(u, v);
    rd = Integral(Grad(u), Grad(v))
       + alpha * Integral(u, v)
       - Integral(f, v)
       + DirichletBC(u, Zero());

    CG(rd).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0.0, 2 * RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    ReactionDiffusion3DTest8,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Wedge)
  );

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    ReactionDiffusion3DTest16,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Wedge)
  );
}
