/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file TemperatureOptimization.cpp
 * @brief Density-based topology optimization with Hilbert gradient regularization.
 *
 * @details
 * This example implements a density-based topology optimization problem for a
 * scalar Poisson equation using the Rodin finite element framework.
 *
 * The optimization problem consists in minimizing a compliance-like functional
 * penalized by material usage:
 *
 * @f[
 *   J(\gamma, u) =
 *   \frac{1}{|\Omega|} \int_\Omega u \, dx
 *   + \ell \int_\Omega \gamma \, dx,
 * @f]
 *
 * subject to the state equation:
 *
 * @f[
 *   -\nabla \cdot \big(a(\gamma)\nabla u\big) = 1 \quad \text{in } \Omega,
 * @f]
 *
 * with Dirichlet boundary conditions on @f$\Gamma_D@f$. The material
 * interpolation is defined by:
 *
 * @f[
 *   a(\gamma) = g_{\min} + (g_{\max} - g_{\min}) \gamma^3.
 * @f]
 *
 * The design variable @f$\gamma \in [0,1]@f$ represents a material density.
 *
 * ---
 *
 * ## Adjoint formulation
 *
 * The adjoint variable @f$p@f$ solves:
 *
 * @f[
 *   -\nabla \cdot \big(a(\gamma)\nabla p\big) = -\frac{1}{|\Omega|}.
 * @f]
 *
 * The gradient of the objective with respect to @f$\gamma@f$ is:
 *
 * @f[
 *   \frac{\partial J}{\partial \gamma}
 *   =
 *   \ell
 *   + 3(g_{\max}-g_{\min}) \gamma^2 \nabla u \cdot \nabla p.
 * @f]
 *
 * ---
 *
 * ## Hilbert gradient regularization
 *
 * To obtain a smooth descent direction, the gradient is filtered through a
 * Hilbert space operator:
 *
 * @f[
 *   \alpha (\nabla g, \nabla w) + (g, w)
 *   =
 *   \left(
 *     \ell
 *     + 3(g_{\max}-g_{\min}) \gamma^2 \nabla u \cdot \nabla p,
 *     w
 *   \right).
 * @f]
 *
 * The resulting function @f$g@f$ is a regularized gradient.
 *
 * ---
 *
 * ## Optimization algorithm
 *
 * The design is updated via a projected gradient descent:
 *
 * @f[
 *   \gamma^{k+1}
 *   =
 *   \Pi_{[0,1]}\left(
 *     \gamma^k - \mu g
 *   \right),
 * @f]
 *
 * where @f$\Pi_{[0,1]}@f$ denotes pointwise projection onto the admissible
 * interval.
 *
 * ---
 *
 * ## Numerical setup
 *
 * - Domain: unit square @f$[0,1]^2@f$
 * - Mesh: structured triangular grid
 * - FE spaces: @f$P_1@f$ Lagrange for state, adjoint, and design
 * - Solver: sparse direct solver (LU)
 *
 * ---
 *
 * ## Output
 *
 * Results are written in XDMF format for visualization in ParaView:
 *
 * - `u`     : state solution
 * - `p`     : adjoint solution
 * - `g`     : regularized gradient
 * - `gamma` : density field
 * - `step`  : descent step
 *
 * ---
 *
 * @note
 * This formulation uses a penalized objective instead of a strict volume
 * constraint. Replacing the term @f$\ell \int \gamma@f$ by a constraint
 * @f$\int \gamma \leq V_0@f$ would require a different optimization strategy
 * (e.g., augmented Lagrangian or projection methods).
 */

#include "Rodin/Context/Local.h"
#include <Rodin/IO.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// Boundary attributes
static constexpr int Gamma0 = 1;
static constexpr int GammaD = 2;

// Optimization parameters
static constexpr Real ell = 4;
static constexpr Real mu = 0.01;
static constexpr Real gmin = 0.0001;
static constexpr Real gmax = 1;
static constexpr Real alpha = 0.05;
static constexpr size_t maxIterations = 5000;
static constexpr Real radius = 0.1;

int main(int, char**)
{
  constexpr size_t n = 32;
  constexpr Geometry::Polytope::Type g = Geometry::Polytope::Type::Triangle;

  Alert::Info() << "Generating uniform grid..." << Alert::Raise;
  Mesh mesh = Mesh<Context::Local>::UniformGrid(g, { n, n });

  mesh.getConnectivity().compute(1, 2);
  mesh.scale(1.0 / (n - 1));

  // Boundary labeling
  for (auto it = mesh.getBoundary(); it; ++it)
  {
    bool onGammaD = false;
    for (auto vit = it->getVertex(); vit; ++vit)
    {
      const auto& vertex = *vit;
      const Real x = vertex.x();
      const Real y = vertex.y();

      if (std::abs(x) < 1e-12 && y > 0.5 - radius && y < 0.5 + radius)
      {
        onGammaD = true;
        break;
      }
    }
    mesh.setAttribute(it.key(), onGammaD ? GammaD : Gamma0);
  }

  IO::XDMF xdmf("TemperatureMinimization");
  xdmf.grid().setMesh(mesh);

  P1 vh(mesh);
  P1 ph(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  TrialFunction p(vh);
  TestFunction  q(vh);

  TrialFunction gfun(vh);
  TestFunction  w(vh);

  GridFunction gamma(ph);
  gamma = 0.9;

  GridFunction step(ph);
  step = 0.0;

  xdmf.add("u", u.getSolution());
  xdmf.add("p", p.getSolution());
  xdmf.add("g", gfun.getSolution());
  xdmf.add("gamma", gamma);
  xdmf.add("step", step);

  const Real vol = mesh.getMeasure(mesh.getDimension());

  for (size_t i = 0; i < maxIterations; i++)
  {
    Alert::Info() << "Iteration: " << i << Alert::Raise;

    const RealFunction f(1.0);
    const auto a = gmin + (gmax - gmin) * Pow(gamma, 3);

    Problem poisson(u, v);
    poisson = Integral(a * Grad(u), Grad(v))
            - Integral(f * v)
            + DirichletBC(u, RealFunction(0.0)).on(GammaD);
    Solver::SparseLU(poisson).solve();

    const Real Ju = Integral(u.getSolution()).compute();
    const Real Jg = Integral(gamma).compute();
    const Real J = Ju / vol + ell * Jg;

    Alert::Info() << "Objective: " << J << Alert::Raise;

    Problem adjoint(p, q);
    adjoint = Integral(a * Grad(p), Grad(q))
            + Integral(RealFunction(1.0 / vol), q)
            + DirichletBC(p, RealFunction(0.0)).on(GammaD);
    Solver::SparseLU(adjoint).solve();

    Problem hilbert(gfun, w);
    hilbert = Integral(alpha * Grad(gfun), Grad(w))
            + Integral(gfun, w)
            - Integral(
                ell
                + 3 * (gmax - gmin) * Pow(gamma, 2)
                  * Dot(Grad(u.getSolution()), Grad(p.getSolution())),
                w)
            + DirichletBC(gfun, RealFunction(0.0)).on(GammaD);
    Solver::SparseLU(hilbert).solve();

    step = mu * gfun.getSolution();

    gamma -= step;
    gamma = Min(1.0, Max(0.0, gamma));

    xdmf.write().flush();
  }

  xdmf.close();
  return 0;
}
