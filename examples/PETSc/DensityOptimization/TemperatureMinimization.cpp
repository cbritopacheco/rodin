/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file PETSc/TemperatureOptimization.cpp
 * @brief Distributed density-based topology optimization using PETSc and MPI.
 *
 * @details
 * This example implements a distributed topology optimization algorithm for a
 * Poisson problem using the Rodin framework with PETSc and MPI backends.
 *
 * The optimization minimizes a compliance-like objective penalized by material
 * usage:
 *
 * @f[
 *   J(\gamma, u)
 *   =
 *   \frac{1}{|\Omega|}
 *   \int_\Omega u \, dx
 *   + \ell \int_\Omega \gamma \, dx,
 * @f]
 *
 * subject to the state equation:
 *
 * @f[
 *   -\nabla \cdot \big(a(\gamma)\nabla u\big) = 1 \quad \text{in } \Omega,
 * @f]
 *
 * with Dirichlet boundary conditions on @f$\Gamma_D@f$.
 *
 * The material interpolation is given by:
 *
 * @f[
 *   a(\gamma) = g_{\min} + (g_{\max} - g_{\min}) \gamma^3.
 * @f]
 *
 * The design variable @f$\gamma \in [0,1]@f$ represents a density field.
 *
 * ---
 *
 * ## Adjoint formulation
 *
 * The adjoint variable @f$p@f$ satisfies:
 *
 * @f[
 *   -\nabla \cdot \big(a(\gamma)\nabla p\big)
 *   =
 *   -\frac{1}{|\Omega|}.
 * @f]
 *
 * The gradient of the objective with respect to @f$\gamma@f$ is:
 *
 * @f[
 *   \frac{\partial J}{\partial \gamma}
 *   =
 *   \ell
 *   + 3(g_{\max} - g_{\min}) \gamma^2 \nabla u \cdot \nabla p.
 * @f]
 *
 * ---
 *
 * ## Hilbert gradient regularization
 *
 * The raw gradient is regularized through a Helmholtz-type operator:
 *
 * @f[
 *   \alpha (\nabla g, \nabla w) + (g, w)
 *   =
 *   \left(
 *     \ell
 *     + 3(g_{\max} - g_{\min}) \gamma^2 \nabla u \cdot \nabla p,
 *     w
 *   \right).
 * @f]
 *
 * This produces a smooth descent direction @f$g@f$.
 *
 * ---
 *
 * ## Optimization algorithm
 *
 * The design is updated iteratively via a projected gradient method:
 *
 * @f[
 *   \gamma^{k+1}
 *   =
 *   \Pi_{[0,1]}\big(\gamma^k \pm \mu g\big).
 * @f]
 *
 * ---
 *
 * ## Parallel implementation
 *
 * - Mesh is distributed across MPI ranks.
 * - PETSc vectors (`Vec`) store degrees of freedom.
 * - Assembly is performed locally with global synchronization.
 * - Linear systems are solved using PETSc KSP solvers.
 *
 * ---
 *
 * ## Numerical setup
 *
 * - Domain: unit square @f$[0,1]^2@f$
 * - Mesh: distributed structured triangular grid
 * - FE spaces: @f$P_1@f$ Lagrange elements
 * - Solver: PETSc Krylov methods (KSP)
 *
 * ---
 *
 * ## Output
 *
 * Results are written in parallel XDMF format:
 *
 * - `u`     : state solution
 * - `p`     : adjoint solution
 * - `g`     : regularized gradient
 * - `gamma` : density field
 * - `step`  : update step
 *
 * ---
 *
 * @note
 * This formulation uses a penalization approach instead of a strict volume
 * constraint. Introducing a constraint @f$\int_\Omega \gamma \le V_0@f$
 * would require a different optimization strategy.
 */

#include "Rodin/Alert/Success.h"
#include <Rodin/IO.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <petscmat.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// Boundary attributes
static constexpr int Gamma0 = 1, GammaD = 2;

// Optimization parameters
static constexpr Real ell = 4;
static constexpr Real mu = 0.01;
static constexpr Real gmin = 0.0001;
static constexpr Real gmax = 1;
static constexpr Real alpha = 0.05;
static constexpr size_t maxIterations = 1e4;
static constexpr Real radius = 0.1;

int main(int argc, char** argv)
{
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
  assert(ierr == PETSC_SUCCESS);

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);

  constexpr size_t n = 512;
  constexpr Geometry::Polytope::Type g = Geometry::Polytope::Type::Quadrilateral;

  if (world.rank() == 0)
    Alert::Info() << "Generating uniform grid..." << Alert::Raise;

  Context::MPI mpi(env, world);
  auto mesh = Mesh<Context::MPI>::UniformGrid(mpi, g, { n, n });

  mesh.getConnectivity().compute(1, 2);
  mesh.reconcile(1);

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

  const size_t nv = mesh.getVertexCount();
  const size_t nc = mesh.getCellCount();

  if (world.rank() == 0)
  {
    Alert::Success() << "Mesh generated." << Alert::Raise;
    Alert::Info() << "Vertices: " << nv << Alert::Raise;
    Alert::Info() << "Cells: " << nc << Alert::Raise;
  }

  IO::XDMF xdmf(world, "TemperatureMinimization");
  xdmf.grid().setMesh(mesh);

  P1 vh(mesh);

  {
    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);
    xdmf.add("u", u.getSolution());

    PETSc::Variational::TrialFunction p(vh);
    PETSc::Variational::TestFunction  q(vh);
    xdmf.add("p", p.getSolution());

    PETSc::Variational::TrialFunction gfun(vh);
    PETSc::Variational::TestFunction  w(vh);
    xdmf.add("g", gfun.getSolution());

    PETSc::Variational::GridFunction gamma(vh);
    gamma = 0.9;
    xdmf.add("gamma", gamma);

    PETSc::Variational::GridFunction step(vh);
    xdmf.add("step", step);

    const Real vol = mesh.getMeasure(mesh.getDimension());

    for (size_t i = 0; i < maxIterations; i++)
    {
      if (world.rank() == 0)
        Alert::Info() << "Iteration: " << i << Alert::Raise;

      const RealFunction f(1.0);
      const auto a = gmin + (gmax - gmin) * Pow(gamma, 3);

      Problem poisson(u, v);
      poisson = Integral(a * Grad(u), Grad(v))
              - Integral(f * v)
              + DirichletBC(u, RealFunction(0.0)).on(GammaD);

      if (world.rank() == 0)
        Alert::Info() << "Assembling state equation..." << Alert::Raise;

      poisson.assemble();

      if (world.rank() == 0)
        Alert::Info() << "State equation assembled." << Alert::Raise;
      Solver::KSP(poisson).solve();

      if (world.rank() == 0)
        Alert::Info() << "State equation solved." << Alert::Raise;

      // const Real Ju = Integral(u.getSolution()).compute();
      // const Real Jg = Integral(gamma).compute();
      // const Real J = Ju / vol + ell * Jg;

      // if (world.rank() == 0)
      //   Alert::Info() << "Objective: " << J << Alert::Raise;

      Problem adjoint(p, q);
      adjoint = Integral(a * Grad(p), Grad(q))
              + Integral(RealFunction(1.0 / vol), q)
              + DirichletBC(p, RealFunction(0.0)).on(GammaD);
      Solver::KSP(adjoint).solve();

      Problem hilbert(gfun, w);
      hilbert = Integral(alpha * Grad(gfun), Grad(w))
              + Integral(gfun, w)
              - Integral(
                  ell
                  + 3 * (gmax - gmin) * Pow(gamma, 2)
                    * Dot(Grad(u.getSolution()), Grad(p.getSolution())),
                  w)
              + DirichletBC(gfun, RealFunction(0.0)).on(GammaD);
      Solver::KSP(hilbert).solve();

      step = mu * gfun.getSolution();

      gamma -= step;

      gamma = Min(1.0, Max(0.0, gamma));

      xdmf.write().flush();
    }
  }

  Alert::Success() << "Optimization completed." << Alert::Raise;

  xdmf.close();
  PetscFinalize();

  return 0;
}
