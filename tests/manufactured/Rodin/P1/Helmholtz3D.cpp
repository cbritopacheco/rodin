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
 * 3D manufactured tests for the Helmholtz equation (sign convention used here):
 *
 * Strong form:
 *   -Δu - κ² u = f   in Ω
 *            u = g   on ∂Ω
 *
 * Weak form:
 *   ∫Ω ∇u·∇v dx - κ² ∫Ω u v dx = ∫Ω f v dx
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
  template <Polytope::Type G, size_t M>
  class Helmholtz3DFixture : public ::testing::Test
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

  using Tetra8  = Helmholtz3DFixture<Polytope::Type::Tetrahedron, 8>;
  using Hex8    = Helmholtz3DFixture<Polytope::Type::Hexahedron,   8>;
  using Tetra16 = Helmholtz3DFixture<Polytope::Type::Tetrahedron, 16>;
  using Hex16   = Helmholtz3DFixture<Polytope::Type::Hexahedron,  16>;
  using Tetra32 = Helmholtz3DFixture<Polytope::Type::Tetrahedron, 32>;

  // ---------------------------------------------------------------------------
  // u = sin(pi x) sin(pi y) sin(pi z)
  // Δu = -3 pi^2 u
  // f = -Δu - κ^2 u = (3 pi^2 - κ^2) u
  // Dirichlet: u=0 on ∂Ω (sine vanishes on x,y,z=0 or 1)
  // ---------------------------------------------------------------------------
  TEST_F(Tetra16, SimpleSine_Tetrahedron)
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
    EXPECT_NEAR(0.1 * error, 0, RODIN_FUZZY_CONSTANT);
  }

  TEST_F(Hex16, SimpleSine_Hexahedron)
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
    EXPECT_NEAR(0.7 * error, 0, RODIN_FUZZY_CONSTANT);
  }

  // ---------------------------------------------------------------------------
  // u = sin(ω1 pi x) sin(ω2 pi y) sin(ω3 pi z)
  // Δu = -(ω1^2 + ω2^2 + ω3^2) pi^2 u
  // f = -Δu - κ^2 u = ( (ω1^2+ω2^2+ω3^2) pi^2 - κ^2 ) u
  // Dirichlet: u=0 on ∂Ω for integer ωi (still vanishes at x,y,z=0 or 1)
  // ---------------------------------------------------------------------------
  TEST_F(Tetra32, VariableFrequency_Tetrahedron)
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
    EXPECT_NEAR(0.1 * error, 0, RODIN_FUZZY_CONSTANT);
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
  TEST_F(Tetra16, MixedPolynomialTrig_Tetrahedron)
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
  TEST_F(Hex16, Exponential_Hexahedron)
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
  // Extra: run the SimpleSine case on Hex16 as requested.
  // ---------------------------------------------------------------------------
  TEST_F(Hex16, SimpleSine_Hexahedron_16)
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
    EXPECT_NEAR(0.7 * error, 0, RODIN_FUZZY_CONSTANT);
  }
}
