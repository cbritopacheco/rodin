/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief DG Poisson example using a P0 interior-penalty scheme.
 *
 * Solves the Poisson problem
 * @f[
 *   -\Delta u = f \quad \text{in } \Omega = [0,1]^2, \qquad
 *   u = g \quad \text{on } \partial\Omega,
 * @f]
 * using a discontinuous Galerkin (DG) penalty method on a piecewise-constant
 * (P0) finite element space.
 *
 * Since @f$ \nabla u|_K = 0 @f$ for P0, the DG bilinear form reduces to
 * a pure interior-penalty (graph-Laplacian / TPFA-type) scheme:
 * @f[
 *   a_h(u,v) =
 *   \sigma \int_{\Gamma_{\mathrm{int}}} [\![u]\!][\![v]\!]\,ds
 *   + \sigma \int_{\partial\Omega} u\,v\,ds
 *   = \int_\Omega f\,v\,dx
 *   + \sigma \int_{\partial\Omega} g\,v\,ds,
 * @f]
 * where @f$ \sigma = 1/h @f$ is the penalty parameter and
 * @f$ h @f$ is the mesh size.  This recovers the two-point flux
 * approximation (TPFA) finite-volume scheme on triangular meshes.
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

  // ---- Manufactured solution: u = sin(π x) sin(π y), f = 2π² u, g = 0 -----
  const Real pi = Math::Constants::pi();
  auto f = 2.0 * pi * pi * sin(pi * F::x) * sin(pi * F::y);
  auto g = Zero();

  // ---- FE space (P0 = discontinuous piecewise constant) --------------------
  P0 vh(mesh);
  TrialFunction u(vh);
  TestFunction  v(vh);

  // Penalty parameter σ = 1/h (TPFA-like scaling)
  const Real h     = 1.0 / static_cast<Real>(N - 1);
  const Real sigma = 1.0 / h;

  // ---- DG interior-penalty bilinear form ------------------------------------
  //
  //  a(u,v) = σ ∫_Γ_int [[u]][[v]] ds           [interior coupling]
  //           + σ ∫_∂Ω  u v ds                   [boundary coupling]
  //
  //  l(v)   = ∫_Ω f v dx                         [source]
  //           + σ ∫_∂Ω g v ds                    [Dirichlet RHS]
  //
  Problem poisson(u, v);
  poisson = InterfaceIntegral(sigma * Jump(u), Jump(v))
          + BoundaryIntegral(sigma * u, v)
          - Integral(f, v)
          - BoundaryIntegral(sigma * g, v);

  // ---- Solve ----------------------------------------------------------------
  Solver::SparseLU(poisson).solve();

  // ---- Save -----------------------------------------------------------------
  u.getSolution().save("DGPoisson.gf", IO::FileFormat::MFEM);
  mesh.save("DGPoisson.mesh", IO::FileFormat::MFEM);

  return 0;
}

