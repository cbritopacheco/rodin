/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file TemperatureMinimization.cpp
 * @brief State, adjoint and projected regularized density updates for heat conduction.
 *
 * The conductivity is @f$a(\gamma)=10^{-4}+(1-10^{-4})\gamma^3@f$.
 * The objective combines mean temperature with the fixed material penalty
 * @f$4\int_\Omega\gamma\,dx@f$. State, adjoint and descent fields are exported
 * for the same density before each projected update.
 *
 * @see examples-density_optimization-simp_temperature_minimization
 */

//! [density-source]
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
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

int main(int argc, char** argv)
{
  size_t iterations = maxIterations;
  if (argc > 2)
  {
    std::cerr << "Usage: " << argv[0] << " [positive-iteration-count]" << std::endl;
    return EXIT_FAILURE;
  }
  if (argc == 2)
  {
    const auto end = argv[1] + std::strlen(argv[1]);
    const auto result = std::from_chars(argv[1], end, iterations);
    if (result.ec != std::errc{} || result.ptr != end || iterations == 0)
      return EXIT_FAILURE;
  }
  //! [density-mesh]
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

  //! [density-mesh]

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

  std::ofstream history("TemperatureMinimization.csv");
  history.precision(17);
  history << "iteration,mean_temperature,material,objective\n";
  for (size_t i = 0; i < iterations; i++)
  {
    Alert::Info() << "Iteration: " << i << Alert::Raise;

    const RealFunction f(1.0);
    const auto a = gmin + (gmax - gmin) * Pow(gamma, 3);

    //! [density-state]
    Problem poisson(u, v);
    poisson = Integral(a * Grad(u), Grad(v))
            - Integral(f * v)
            + DirichletBC(u, RealFunction(0.0)).on(GammaD);
    Solver::SparseLU stateSolver(poisson);
    stateSolver.solve();
    if (!stateSolver.success())
      return EXIT_FAILURE;

    //! [density-state]

    const Real Ju = Integral(u.getSolution()).compute();
    const Real Jg = Integral(gamma).compute();
    const Real J = Ju / vol + ell * Jg;

    if (!std::isfinite(J))
      return EXIT_FAILURE;
    history << i << ',' << Ju / vol << ',' << Jg << ',' << J << '\n';
    Alert::Info() << "Objective: " << J << Alert::Raise;

    //! [density-adjoint]
    Problem adjoint(p, q);
    adjoint = Integral(a * Grad(p), Grad(q))
            + Integral(RealFunction(1.0 / vol), q)
            + DirichletBC(p, RealFunction(0.0)).on(GammaD);
    Solver::SparseLU adjointSolver(adjoint);
    adjointSolver.solve();
    if (!adjointSolver.success())
      return EXIT_FAILURE;

    //! [density-adjoint]

    //! [density-gradient]
    Problem hilbert(gfun, w);
    hilbert = Integral(alpha * Grad(gfun), Grad(w))
            + Integral(gfun, w)
            - Integral(
                ell
                + 3 * (gmax - gmin) * Pow(gamma, 2)
                  * Dot(Grad(u.getSolution()), Grad(p.getSolution())),
                w)
            + DirichletBC(gfun, RealFunction(0.0)).on(GammaD);
    Solver::SparseLU gradientSolver(hilbert);
    gradientSolver.solve();
    if (!gradientSolver.success())
      return EXIT_FAILURE;

    //! [density-gradient]

    //! [density-update]
    step = mu * gfun.getSolution();

    // Snapshot the state, adjoint and descent computed for this same density.
    xdmf.write(Real(i)).flush();
    gamma -= step;
    gamma = Min(1.0, Max(0.0, gamma));
    //! [density-update]
  }

  xdmf.close();
  history.flush();
  return history ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [density-source]
