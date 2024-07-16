/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// Define boundary attributes
static constexpr int Gamma0 = 1, GammaD = 2;

// Optimization parameters
static constexpr double ell = 4;
static constexpr double mu = 0.2;
static constexpr double gmin = 0.0001;
static constexpr double gmax = 1;
static constexpr double alpha = 0.05;
static constexpr size_t maxIterations = 5000;

int main(int, char**)
{
  const char* meshFile = "../resources/examples/DensityOptimization/TemperatureMinimization.mfem.mesh";

  // Load mesh
  External::MMG::Mesh Omega;
  Omega.load(meshFile);
  External::MMG::Optimizer().setHMax(0.005).setHMin(0.001).optimize(Omega);

  Omega.save("outDensity/gamma.mesh");

  // Build finite element space
  P1 vh(Omega);
  P1 ph(Omega);

  GridFunction gamma(ph);
  gamma = 0.9;

  double vol = Omega.getVolume();

  // Optimization loop
  for (size_t i = 0; i < maxIterations; i++)
  {
    Alert::Info() << "Iteration: " << i << Alert::Raise;

    // Poisson problem
    RealFunction f(1.0);

    TrialFunction u(vh);
    TestFunction  v(vh);
    Problem poisson(u, v);
    poisson = Integral((gmin + (gmax - gmin) * Pow(gamma, 3)) * Grad(u), Grad(v))
            - Integral(f * v)
            + DirichletBC(u, RealFunction(0.0)).on(GammaD);
    Solver::SparseLU(poisson).solve();

    // Adjoint problem
    TrialFunction p(vh);
    TestFunction  q(vh);
    Problem adjoint(p, q);
    adjoint = Integral((gmin + (gmax - gmin) * Pow(gamma, 3)) * Grad(p), Grad(q))
            + Integral(RealFunction(1.0 / vol), q)
            + DirichletBC(p, RealFunction(0.0)).on(GammaD);
    Solver::SparseLU(adjoint).solve();

    // Hilbert extension-regularization
    TrialFunction g(vh);
    TestFunction  w(vh);
    Problem hilbert(g, w);
    hilbert = Integral(alpha * Grad(g), Grad(w))
            + Integral(g, w)
            - Integral(
                ell + 3 * (gmax - gmin) * Pow(gamma, 2) *
                Dot(Grad(u.getSolution()), Grad(p.getSolution())), w)
            + DirichletBC(g, RealFunction(0.0)).on(GammaD);
    Solver::SparseLU(hilbert).solve();

    GridFunction step(ph);
    step = mu * g.getSolution();

    gamma -= step;
    gamma = Min(1.0, Max(0.0, gamma));

    Omega.save("gamma.mesh");
    gamma.save("gamma.gf");
    gamma.save("outDensity/gamma." + std::to_string(i) + ".gf");
  }

  Omega.save("gamma.mesh");
  gamma.save("gamma.gf");

  return 0;
}
