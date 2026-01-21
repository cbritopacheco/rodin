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
 * @brief 3D Manufactured solutions for reaction-diffusion equation using P1 spaces.
 *
 * The reaction-diffusion equation in 3D:
 * @f[
 *   -\Delta u + \alpha u = f \quad \text{in } \Omega \subset \mathbb{R}^3
 * @f]
 * with Dirichlet boundary conditions
 * @f[
 *   u = g \quad \text{on } \partial\Omega.
 * @f]
 *
 * Weak formulation: Find @f$ u \in H^1_0(\Omega) @f$ such that
 * @f[
 *   \int_\Omega \nabla u \cdot \nabla v \,dx + \alpha \int_\Omega u v \,dx = \int_\Omega f v \,dx,
 * @f]
 * for all @f$ v \in H^1_0(\Omega) @f$.
 */
namespace Rodin::Tests::Manufactured::ReactionDiffusion3D
{
  /**
   * @brief Simple sine manufactured solution
   * u(x,y,z) = sin(πx) sin(πy) sin(πz)
   */
  TEST(ReactionDiffusion3D, SimpleSine_Tetrahedron)
  {
    auto pi = Math::Constants::pi();
    const Real alpha = 1.0; // reaction coefficient

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    // f = -Δu + αu = 3π²u + αu
    auto f = (3 * pi * pi + alpha) * u_expr;

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
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Simple sine on Hexahedron mesh
   */
  TEST(ReactionDiffusion3D, SimpleSine_Hexahedron)
  {
    auto pi = Math::Constants::pi();
    const Real alpha = 1.0;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Hexahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto f = (3 * pi * pi + alpha) * u_expr;

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
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Polynomial manufactured solution
   * u(x,y,z) = x(1-x) y(1-y) z(1-z)
   */
  TEST(ReactionDiffusion3D, Polynomial_Tetrahedron)
  {
    const Real alpha = 2.0;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = F::x * (1 - F::x) * F::y * (1 - F::y) * F::z * (1 - F::z);
    // Δu = 2y(1-y)z(1-z) + 2x(1-x)z(1-z) + 2x(1-x)y(1-y)
    auto laplacian_u = 2 * F::y * (1 - F::y) * F::z * (1 - F::z)
                     + 2 * F::x * (1 - F::x) * F::z * (1 - F::z)
                     + 2 * F::x * (1 - F::x) * F::y * (1 - F::y);
    auto f = -laplacian_u + alpha * u_expr;

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
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief High reaction coefficient test
   * Tests stability with large reaction term
   */
  TEST(ReactionDiffusion3D, HighReaction_Tetrahedron)
  {
    auto pi = Math::Constants::pi();
    const Real alpha = 10.0; // high reaction coefficient

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {16, 16, 16});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto f = (3 * pi * pi + alpha) * u_expr;

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
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Mixed polynomial-trigonometric solution
   * u(x,y,z) = x(1-x) sin(πy) sin(πz)
   */
  TEST(ReactionDiffusion3D, MixedSolution_Hexahedron)
  {
    auto pi = Math::Constants::pi();
    const Real alpha = 0.5;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Hexahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = F::x * (1 - F::x) * sin(pi * F::y) * sin(pi * F::z);
    // Δu = -2 sin(πy) sin(πz) - 2π² x(1-x) sin(πy) sin(πz)
    auto laplacian_u = -2 * sin(pi * F::y) * sin(pi * F::z)
                     - 2 * pi * pi * F::x * (1 - F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto f = -laplacian_u + alpha * u_expr;

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
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Exponential solution with nonhomogeneous boundary
   * u(x,y,z) = cos(πx) cos(πy) e^z
   */
  TEST(ReactionDiffusion3D, Exponential_Tetrahedron)
  {
    auto pi = Math::Constants::pi();
    const Real alpha = 1.5;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    P1 vh(mesh);

    auto u_expr = cos(pi * F::x) * cos(pi * F::y) * exp(F::z);
    // Δu = -2π² cos(πx) cos(πy) e^z + cos(πx) cos(πy) e^z
    //    = (1 - 2π²) cos(πx) cos(πy) e^z
    auto laplacian_u = (1 - 2 * pi * pi) * u_expr;
    auto f = -laplacian_u + alpha * u_expr;

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
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }
}
