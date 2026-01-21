/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

/**
 * @brief 3D Manufactured solutions for the Helmholtz equation using P1 spaces.
 *
 * The Helmholtz equation in 3D:
 * @f[
 *   -\Delta u - \kappa^2 u = f \quad \text{in } \Omega \subset \mathbb{R}^3
 * @f]
 * with Dirichlet boundary conditions
 * @f[
 *   u = g \quad \text{on } \partial\Omega.
 * @f]
 *
 * Weak formulation: Find @f$ u \in H^1_0(\Omega) @f$ such that
 * @f[
 *   \int_\Omega \nabla u \cdot \nabla v \,dx - \kappa^2 \int_\Omega u v \,dx = \int_\Omega f v \,dx,
 * @f]
 * for all @f$ v \in H^1_0(\Omega) @f$.
 */
namespace Rodin::Tests::Manufactured::Helmholtz3D
{
  /**
   * @brief Test using SimpleSine manufactured solution
   * u(x,y,z) = sin(πx) sin(πy) sin(πz)
   */
  TEST(Helmholtz3D, SimpleSine_Tetrahedron)
  {
    auto pi = Math::Constants::pi();
    const Real kappa = 2.0; // wavenumber

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto f = (3 * pi * pi - kappa * kappa) * u_expr;

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

  /**
   * @brief Test using SimpleSine manufactured solution on Hexahedron mesh
   */
  TEST(Helmholtz3D, SimpleSine_Hexahedron)
  {
    auto pi = Math::Constants::pi();
    const Real kappa = 2.0;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Hexahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto f = (3 * pi * pi - kappa * kappa) * u_expr;

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

  /**
   * @brief Test with variable frequency in each direction
   * u(x,y,z) = sin(ω₁πx) sin(ω₂πy) sin(ω₃πz)
   */
  TEST(Helmholtz3D, VariableFrequency_Tetrahedron)
  {
    auto pi = Math::Constants::pi();
    const Real kappa = 1.5;

    const Real omega1 = 2;
    const Real omega2 = 2;
    const Real omega3 = 3;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {16, 16, 16});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = sin(omega1 * pi * F::x) * sin(omega2 * pi * F::y) * sin(omega3 * pi * F::z);
    auto f = ((omega1 * omega1 + omega2 * omega2 + omega3 * omega3) * pi * pi - kappa * kappa) * u_expr;

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

  /**
   * @brief Test with polynomial-trigonometric mixed solution
   * u(x,y,z) = x(1-x) sin(πy) sin(πz)
   */
  TEST(Helmholtz3D, MixedPolynomialTrig_Tetrahedron)
  {
    auto pi = Math::Constants::pi();
    const Real kappa = 1.0;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = F::x * (1 - F::x) * sin(pi * F::y) * sin(pi * F::z);
    // f = -Δu - κ²u
    // Δu = -2 sin(πy) sin(πz) - x(1-x)(π² sin(πy) sin(πz) + π² sin(πy) sin(πz))
    //    = -2 sin(πy) sin(πz) - 2π² x(1-x) sin(πy) sin(πz)
    auto f = 2 * sin(pi * F::y) * sin(pi * F::z)
           + 2 * pi * pi * F::x * (1 - F::x) * sin(pi * F::y) * sin(pi * F::z)
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

  /**
   * @brief Test with exponential component
   * u(x,y,z) = sin(πx) sin(πy) e^z
   */
  TEST(Helmholtz3D, Exponential_Hexahedron)
  {
    auto pi = Math::Constants::pi();
    const Real kappa = 0.5;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Hexahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = sin(pi * F::x) * sin(pi * F::y) * exp(F::z);
    // Δu = -2π² sin(πx) sin(πy) e^z + sin(πx) sin(πy) e^z
    //    = (1 - 2π²) sin(πx) sin(πy) e^z
    auto f = (2 * pi * pi - 1 - kappa * kappa) * u_expr;

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
}
