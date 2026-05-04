/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief DG advection example using an upwind P0 DG scheme.
 *
 * Solves the steady linear advection problem
 * @f[
 *   \boldsymbol{\beta} \cdot \nabla u = f \quad \text{in } \Omega = [0,1]^2,
 *   \qquad u = g \quad \text{on } \partial\Omega_{\mathrm{in}},
 * @f]
 * where @f$ \partial\Omega_{\mathrm{in}} = \{ x \in \partial\Omega :
 * \boldsymbol{\beta} \cdot \mathbf{n} < 0 \} @f$, using a piecewise-constant
 * (P0) finite element space.
 *
 * Since @f$ \nabla v|_K = 0 @f$ for P0, the DG weak form is assembled
 * element-by-element via interface integrals (no volume gradient term):
 * @f[
 *   \sum_{E \in \Gamma_{\mathrm{int}}} \int_E
 *     \Bigl[ (\boldsymbol{\beta}\cdot\mathbf{n})\,\{\!\{u\}\!\}
 *           + \tfrac{1}{2}|\boldsymbol{\beta}\cdot\mathbf{n}|\,[\![u]\!]
 *     \Bigr] [\![v]\!]\,ds
 *   + \int_{\partial\Omega} (\boldsymbol{\beta}\cdot\mathbf{n})^+\, u\, v\, ds
 *   = \int_\Omega f\,v\,dx
 *   - \int_{\partial\Omega} (\boldsymbol{\beta}\cdot\mathbf{n})^-\, g\, v\, ds,
 * @f]
 * where @f$ a^+ = \max(a,0) @f$, @f$ a^- = \min(a,0) @f$, @f$
 * \{\!\{u\}\!\} @f$ and @f$ [\![u]\!] @f$ are the average and jump across
 * the face, and @f$ \mathbf{n} @f$ is the face normal oriented outward from
 * the first adjacent cell.
 *
 * The face-normal orientation is fixed by assigning attribute 1 to all cells
 * and calling @c FaceNormal.traceOf(1), which consistently picks the outward
 * direction from the first (attribute-1) adjacent cell.
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

  // Assign attribute 1 to every cell so that FaceNormal.traceOf(1) gives a
  // consistent orientation on all interior faces.
  for (auto it = mesh.getCell(); it; ++it)
    mesh.setAttribute({ mesh.getDimension(), it->getIndex() }, 1);

  // ---- Advection velocity: β = (1, 0) --------------------------------------
  Math::Vector<Real> betaVec(2);
  betaVec << 1.0, 0.0;
  VectorFunction beta(betaVec);

  // ---- Manufactured solution: u = x(1-x) y(1-y) ----------------------------
  //   β·∇u = (1 − 2x) y(1−y)
  auto f = (1.0 - 2.0 * F::x) * F::y * (1.0 - F::y);
  auto g = Zero();

  // ---- FE space (P0 = discontinuous piecewise constant) --------------------
  P0 vh(mesh);
  TrialFunction u(vh);
  TestFunction  v(vh);

  // Face normals:
  //   fn — interior faces, oriented outward from cells with attribute 1
  //   bn — boundary faces, always oriented outward from the domain
  FaceNormal     fn(mesh);
  BoundaryNormal bn(mesh);

  auto fn1   = fn.traceOf(1);
  auto betaN = Dot(beta, fn1);
  auto betaBN = Dot(beta, bn);

  // ---- Upwind DG bilinear form (interface-integral formulation) -------------
  //
  //  a(u,v) = ∫_Γ_int  (β·n) {{u}} [[v]] ds       [advective average flux]
  //         + ∫_Γ_int  |β·n|/2 [[u]] [[v]] ds      [upwind dissipation]
  //         + ∫_∂Ω  (β·n)⁺ u v ds                 [outflow boundary]
  //
  //  l(v)   = ∫_Ω  f v dx
  //         − ∫_∂Ω (β·n)⁻ g v ds                  [inflow BC]
  //
  Problem advection(u, v);
  advection = InterfaceIntegral(betaN * Average(u), Jump(v))
            + InterfaceIntegral(Abs(betaN) / 2.0 * Jump(u), Jump(v))
            + BoundaryIntegral(Max(betaBN, Zero()) * u, v)
            - Integral(f, v)
            - BoundaryIntegral(Min(betaBN, Zero()) * g, v);

  // ---- Solve ----------------------------------------------------------------
  Solver::SparseLU(advection).solve();

  // ---- Save -----------------------------------------------------------------
  u.getSolution().save("DGAdvection.gf", IO::FileFormat::MFEM);
  mesh.save("DGAdvection.mesh", IO::FileFormat::MFEM);

  return 0;
}

