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
 *   ∫Ω ∇u·∇v dx + α ∫Ω u v dx = ∫Ω f v dx
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
  template <Polytope::Type G, size_t M>
  class ReactionDiffusion3DFixture : public ::testing::Test
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh<Context::Local>().UniformGrid(G, {M, M, M});
        m_mesh.scale(Real(1) / Real(M - 1)); // map {0..M-1} -> [0,1]
        m_mesh.getConnectivity().compute(2, 3);
      }

      const Mesh<Context::Local>& mesh() const { return m_mesh; }

    private:
      Mesh<Context::Local> m_mesh;
  };

  using Tetra8  = ReactionDiffusion3DFixture<Polytope::Type::Tetrahedron, 8>;
  using Hex8    = ReactionDiffusion3DFixture<Polytope::Type::Hexahedron,   8>;
  using Tetra16 = ReactionDiffusion3DFixture<Polytope::Type::Tetrahedron, 16>;
  using Hex16   = ReactionDiffusion3DFixture<Polytope::Type::Hexahedron,  16>;
  using Tetra32 = ReactionDiffusion3DFixture<Polytope::Type::Tetrahedron, 32>;

  // ---------------------------------------------------------------------------
  // u = sin(pi x) sin(pi y) sin(pi z)
  // Δu = -3 pi^2 u
  // f = -Δu + α u = (3 pi^2 + α) u
  // ---------------------------------------------------------------------------
  TEST_F(Tetra16, SimpleSine_Tetrahedron)
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
    EXPECT_NEAR(Integral(diff).compute(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST_F(Hex16, SimpleSine_Hexahedron)
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
    EXPECT_NEAR(0.7 * Integral(diff).compute(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // Polynomial:
  // u = x(1-x)y(1-y)z(1-z)
  // Δu = 2y(1-y)z(1-z) + 2x(1-x)z(1-z) + 2x(1-x)y(1-y)
  // f = -Δu + α u
  // ---------------------------------------------------------------------------
  TEST_F(Tetra16, Polynomial_Tetrahedron)
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
  TEST_F(Hex16, MixedPolynomialTrig_Hexahedron)
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
  TEST_F(Tetra16, Exponential_Tetrahedron)
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
    EXPECT_NEAR(0.7 * Integral(diff).compute(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // Extra: Hex16 refinement sanity check
  // ---------------------------------------------------------------------------
  TEST_F(Hex16, SimpleSine_Hexahedron_16)
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
    EXPECT_NEAR(0.7 * Integral(diff).compute(), 0.0, RODIN_FUZZY_CONSTANT);
  }
}
