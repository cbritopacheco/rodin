/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file PETSc/TemperatureMinimization3D.cpp
 * @brief Distributed 3D density-based topology optimization using PETSc and MPI.
 *
 * @details
 * This example implements a distributed topology optimization algorithm for a
 * Poisson problem in a cube using the Rodin framework with PETSc and MPI
 * backends.
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
 * - Domain: unit cube @f$[0,1]^3@f$
 * - Mesh: distributed structured tetrahedral grid
 * - Dirichlet boundary: centered patch on the bottom face @f$z = 0@f$
 * - FE spaces: @f$P_1@f$ Lagrange elements
 * - Solver: PETSc Krylov methods (KSP)
 *
 * Command-line options:
 *
 * @code
 * --n=<vertices-per-axis>
 * --iters=<optimization-iterations>
 * --ell=<volume-penalty>
 * --alpha=<Hilbert-regularization>
 * --mu=<gradient-step-size>
 * --radius=<bottom-patch-radius>
 * --output-every=<iteration-stride>
 * @endcode
 *
 * ---
 *
 * ## Output
 *
 * Results are written under @c RODIN_EXAMPLE_OUTPUT_DIR in parallel XDMF
 * format:
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

#include <filesystem>
#include <string>
#include <vector>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

#ifndef RODIN_EXAMPLE_OUTPUT_DIR
#define RODIN_EXAMPLE_OUTPUT_DIR "tmp"
#endif

// Boundary attributes
static constexpr int Gamma0 = 1, GammaD = 2;

namespace
{
  size_t parseSizeTOption(int argc, char** argv, const std::string& name, size_t fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; i++)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) == 0)
        return static_cast<size_t>(std::stoul(arg.substr(prefix.size())));
    }
    return fallback;
  }

  Real parseRealOption(int argc, char** argv, const std::string& name, Real fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; i++)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) == 0)
        return static_cast<Real>(std::stod(arg.substr(prefix.size())));
    }
    return fallback;
  }

  bool isRodinOption(const std::string& arg)
  {
    return arg.rfind("--n=", 0) == 0 || arg.rfind("--iters=", 0) == 0 ||
      arg.rfind("--ell=", 0) == 0 || arg.rfind("--alpha=", 0) == 0 ||
      arg.rfind("--mu=", 0) == 0 || arg.rfind("--radius=", 0) == 0 ||
      arg.rfind("--output-every=", 0) == 0;
  }
}

// Optimization parameter defaults
static constexpr Real defaultEll = 4;
static constexpr Real defaultMu = 0.01;
static constexpr Real gmin = 0.0001;
static constexpr Real gmax = 1;
static constexpr Real defaultAlpha = 0.05;
static constexpr size_t defaultMaxIterations = 1e4;
static constexpr size_t defaultN = 64;
static constexpr Real defaultRadius = 0.1;
static constexpr Real eps = 1e-12;

int main(int argc, char** argv)
{
  const size_t requestedN = parseSizeTOption(argc, argv, "n", defaultN);
  const size_t n = requestedN < 2 ? 2 : requestedN;
  const size_t maxIterations =
    parseSizeTOption(argc, argv, "iters", defaultMaxIterations);
  const size_t requestedOutputEvery = parseSizeTOption(argc, argv, "output-every", 1);
  const size_t outputEvery = requestedOutputEvery == 0 ? 1 : requestedOutputEvery;
  const Real ell = parseRealOption(argc, argv, "ell", defaultEll);
  const Real alpha = parseRealOption(argc, argv, "alpha", defaultAlpha);
  const Real mu = parseRealOption(argc, argv, "mu", defaultMu);
  const Real radius = parseRealOption(argc, argv, "radius", defaultRadius);

  std::vector<char*> petscArgv;
  petscArgv.reserve(static_cast<size_t>(argc));
  petscArgv.push_back(argv[0]);
  for (int i = 1; i < argc; i++)
  {
    if (!isRodinOption(argv[i]))
      petscArgv.push_back(argv[i]);
  }
  int petscArgc = static_cast<int>(petscArgv.size());
  char** petscArgvData = petscArgv.data();

  PetscErrorCode ierr;
  ierr = PetscInitialize(&petscArgc, &petscArgvData, PETSC_NULLPTR, PETSC_NULLPTR);
  assert(ierr == PETSC_SUCCESS);

  boost::mpi::environment env(petscArgc, petscArgvData);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);

  constexpr Geometry::Polytope::Type g = Geometry::Polytope::Type::Tetrahedron;

  if (world.rank() == 0)
  {
    Alert::Info() << "Generating uniform grid..." << Alert::Raise;
    Alert::Info() << "n=" << n << ", iters=" << maxIterations << ", ell=" << ell
                  << ", alpha=" << alpha << ", mu=" << mu << ", radius=" << radius
                  << ", output-every=" << outputEvery << Alert::Raise;
  }

  Context::MPI mpi(env, world);
  auto mesh = Mesh<Context::MPI>::UniformGrid(mpi, g, {n, n, n});

  mesh.getConnectivity().compute(2, 3);
  mesh.reconcile(2);

  mesh.scale(Real(1) / static_cast<Real>(n - 1));

  // Boundary labeling
  for (auto it = mesh.getBoundary(); it; ++it)
  {
    bool onBottom = true;
    Real xMean = 0;
    Real yMean = 0;
    size_t vertexCount = 0;
    for (auto vit = it->getVertex(); vit; ++vit)
    {
      const auto& vertex = *vit;
      if (std::abs(vertex.z()) >= eps)
        onBottom = false;
      xMean += vertex.x();
      yMean += vertex.y();
      vertexCount++;
    }
    xMean /= vertexCount;
    yMean /= vertexCount;

    const Real dx = xMean - 0.5;
    const Real dy = yMean - 0.5;
    const bool onGammaD = onBottom && dx * dx + dy * dy < radius * radius;
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

  const std::filesystem::path outputDirectory =
    std::filesystem::path(RODIN_EXAMPLE_OUTPUT_DIR) / "TemperatureMinimization3D";
  if (world.rank() == 0)
    std::filesystem::create_directories(outputDirectory);
  world.barrier();

  IO::XDMF xdmf(world, (outputDirectory / "TemperatureMinimization3D").string());
  xdmf.grid().setMesh(mesh);

  P1 vh(mesh);

  {
    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction v(vh);
    xdmf.add("u", u.getSolution());

    PETSc::Variational::TrialFunction p(vh);
    PETSc::Variational::TestFunction q(vh);
    xdmf.add("p", p.getSolution());

    PETSc::Variational::TrialFunction gfun(vh);
    PETSc::Variational::TestFunction w(vh);
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
      poisson = Integral(a * Grad(u), Grad(v)) - Integral(f * v) +
        DirichletBC(u, RealFunction(0.0)).on(GammaD);

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
      adjoint = Integral(a * Grad(p), Grad(q)) + Integral(RealFunction(1.0 / vol), q) +
        DirichletBC(p, RealFunction(0.0)).on(GammaD);
      Solver::KSP(adjoint).solve();

      Problem hilbert(gfun, w);
      hilbert = Integral(alpha * Grad(gfun), Grad(w)) + Integral(gfun, w) -
        Integral(ell +
            3 * (gmax - gmin) * Pow(gamma, 2) *
              Dot(Grad(u.getSolution()), Grad(p.getSolution())),
          w) +
        DirichletBC(gfun, RealFunction(0.0)).on(GammaD);
      Solver::KSP(hilbert).solve();

      step = mu * gfun.getSolution();

      gamma -= step;

      gamma = Min(1.0, Max(0.0, gamma));

      if (i % outputEvery == 0)
        xdmf.write().flush();
    }
  }

  Alert::Success() << "Optimization completed." << Alert::Raise;

  xdmf.close();
  PetscFinalize();

  return 0;
}
