/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief DG Poisson example using the symmetric interior-penalty (Nitsche) method.
 *
 * Solves the Poisson problem with Dirichlet boundary conditions using a
 * Nitsche (SIPG) DG formulation on a P1 space:
 *
 * @f[
 *   -\Delta u = f \quad \text{in } \Omega = [0,1]^2, \qquad
 *   u = g \quad \text{on } \partial\Omega.
 * @f]
 *
 * The Nitsche weak formulation reads: find @f$ u \in H^1(\Omega) @f$ such that
 * @f[
 *   \underbrace{\int_\Omega \nabla u \cdot \nabla v \, dx}_{\text{stiffness}}
 *   - \underbrace{\int_{\partial\Omega} (\nabla u \cdot \mathbf{n})\, v \, ds}_{\text{consistency}}
 *   - \underbrace{\int_{\partial\Omega} u\, (\nabla v \cdot \mathbf{n}) \, ds}_{\text{symmetry}}
 *   + \underbrace{\frac{\sigma}{h} \int_{\partial\Omega} u \, v \, ds}_{\text{penalty}}
 *   = \int_\Omega f \, v \, dx
 *   - \int_{\partial\Omega} g\, (\nabla v \cdot \mathbf{n}) \, ds
 *   + \frac{\sigma}{h} \int_{\partial\Omega} g \, v \, ds,
 * @f]
 * for all @f$ v \in H^1(\Omega) @f$.
 *
 * The penalty parameter @f$ \sigma / h @f$ must be chosen large enough for
 * coercivity.  A safe choice is @f$ \sigma \geq C_0 \, k^2 @f$ (here @f$ k=1
 * @f$ for P1) times the inverse mesh size.
 */
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  // ---- Mesh -----------------------------------------------------------------
  static constexpr size_t N = 32;
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, {N, N});
  mesh.scale(1.0 / (N - 1));
  mesh.getConnectivity().compute(1, 2);

  // ---- Manufactured solution: u = sin(π x) sin(π y), g = 0 ---------------
  const Real pi = Math::Constants::pi();
  auto f = 2.0 * pi * pi * sin(pi * F::x) * sin(pi * F::y);
  auto g = Zero();

  // ---- FE space and functions -----------------------------------------------
  P1 vh(mesh);
  TrialFunction u(vh);
  TestFunction  v(vh);

  BoundaryNormal n(mesh);

  // Penalty parameter: σ/h with σ = 10 (coercivity constant) and h = 1/(N-1)
  const Real h     = 1.0 / static_cast<Real>(N - 1);
  const Real sigma = 10.0 / h;

  // ---- Nitsche / SIPG bilinear form -----------------------------------------
  //
  //  a(u,v) = ∫_Ω ∇u·∇v dx
  //           − ∫_∂Ω (∇u·n) v ds       [consistency]
  //           − ∫_∂Ω u (∇v·n) ds       [symmetry]
  //           + σ/h ∫_∂Ω u v ds         [penalty]
  //
  //  l(v)   = ∫_Ω f v dx
  //           − ∫_∂Ω g (∇v·n) ds       [RHS consistency with BC]
  //           + σ/h ∫_∂Ω g v ds         [RHS penalty with BC]
  //
  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - BoundaryIntegral(Dot(Grad(u), n), v)
          - BoundaryIntegral(u, Dot(Grad(v), n))
          + BoundaryIntegral(sigma * u, v)
          - Integral(f, v)
          - BoundaryIntegral(g, Dot(Grad(v), n))
          + BoundaryIntegral(sigma * g, v);

  // ---- Solve ----------------------------------------------------------------
  Solver::SparseLU(poisson).solve();

  // ---- Save -----------------------------------------------------------------
  u.getSolution().save("DGPoisson.gf", IO::FileFormat::MFEM);
  mesh.save("DGPoisson.mesh", IO::FileFormat::MFEM);

  return 0;
}
