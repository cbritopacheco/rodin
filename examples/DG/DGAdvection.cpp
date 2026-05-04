/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief DG advection example using an upwind DG formulation.
 *
 * Solves the steady linear advection problem
 * @f[
 *   \boldsymbol{\beta} \cdot \nabla u = f \quad \text{in } \Omega = [0,1]^2,
 *   \qquad u = g \quad \text{on } \partial\Omega_{\mathrm{in}},
 * @f]
 * where @f$ \partial\Omega_{\mathrm{in}} = \{ x \in \partial\Omega :
 * \boldsymbol{\beta} \cdot \mathbf{n} < 0 \} @f$ is the inflow boundary.
 *
 * The DG weak formulation is obtained by integration by parts over the whole
 * domain and upwinding at the boundary:
 * @f[
 *   - \int_\Omega u\,(\boldsymbol{\beta} \cdot \nabla v)\,dx
 *   + \int_{\partial\Omega} (\boldsymbol{\beta} \cdot \mathbf{n})^+ \, u \, v \, ds
 *   = \int_\Omega f \, v \, dx
 *   - \int_{\partial\Omega} (\boldsymbol{\beta} \cdot \mathbf{n})^- \, g \, v \, ds,
 * @f]
 * for all @f$ v @f$, where @f$ (\cdot)^+ = \max(\cdot,\,0) @f$ (outflow) and
 * @f$ (\cdot)^- = \min(\cdot,\,0) @f$ (inflow).  The upwind splitting
 * naturally handles inflow/outflow without explicit attribute restrictions.
 *
 * Manufactured solution: @f$ u = x(1-x)\,y(1-y) @f$ with
 * @f$ \boldsymbol{\beta} = (1,0)^T @f$.
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

  // ---- Advection velocity: β = (1, 0) --------------------------------------
  Math::Vector<Real> betaVec(2);
  betaVec << 1.0, 0.0;
  VectorFunction beta(betaVec);

  // ---- Manufactured solution: u = x(1-x) y(1-y) ----------------------------
  //   β·∇u = ∂u/∂x = (1 − 2x) y(1−y)
  auto f = (1.0 - 2.0 * F::x) * F::y * (1.0 - F::y);

  // Inflow BC: g = u_exact|_{∂Ω_in}.  For β=(1,0) the inflow is x=0,
  // where u_exact = 0.  The Max/Min splitting below enforces g = 0 on the
  // inflow naturally.
  auto g = Zero();

  // ---- FE space and functions -----------------------------------------------
  P1 vh(mesh);
  TrialFunction u(vh);
  TestFunction  v(vh);

  BoundaryNormal n(mesh);

  // ---- Upwind DG bilinear form (IBP) ----------------------------------------
  //
  //  a(u,v) = −∫_Ω u (β·∇v) dx
  //           + ∫_∂Ω max(β·n, 0) u v ds    [outflow upwind]
  //
  //  l(v)   = ∫_Ω f v dx
  //           − ∫_∂Ω min(β·n, 0) g v ds    [inflow BC]
  //
  // On the outflow boundary β·n > 0:  max(β·n,0) contributes to LHS (u unknown)
  // On the inflow boundary  β·n < 0:  min(β·n,0) contributes to RHS (g known)
  // On glancing boundaries  β·n = 0:  both terms vanish (no contribution)
  //
  Problem advection(u, v);
  advection = -Integral(u, Dot(beta, Grad(v)))
            + BoundaryIntegral(Max(Dot(beta, n), Zero()) * u, v)
            - Integral(f, v)
            - BoundaryIntegral(Min(Dot(beta, n), Zero()) * g, v);

  // ---- Solve ----------------------------------------------------------------
  Solver::SparseLU(advection).solve();

  // ---- Save -----------------------------------------------------------------
  u.getSolution().save("DGAdvection.gf", IO::FileFormat::MFEM);
  mesh.save("DGAdvection.mesh", IO::FileFormat::MFEM);

  return 0;
}
