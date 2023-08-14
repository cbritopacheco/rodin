/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// Define boundary attributes
static constexpr int Gamma0 = 1, GammaD = 2;

// Optimization parameters
static constexpr double ell = 1;
static constexpr double mu = 0.2;
static constexpr double hmin = 0.01;
static constexpr double hmax = 10;
static constexpr double alpha = 0.05;
static constexpr size_t maxIterations = 1000;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/temperature-minimization-example.mesh";

  // Load mesh
  Mesh Omega;
  Omega.load(meshFile);

  Omega.save("outDensity/gamma.mesh");

  // Build finite element space
  H1 Vh(Omega);
  L2 Ph(Omega);

  GridFunction gamma(Ph);
  gamma = 0.9;

  double vol = Omega.getVolume();

  Solver::CG solver;
  solver.printIterations(false);

  // Optimization loop
  for (size_t i = 0; i < maxIterations; i++)
  {
   Alert::Info() << "Iteration: " << i << Alert::Raise;

   // Poisson problem
   ScalarFunction f(1.0);

   TrialFunction u(Vh);
   TestFunction  v(Vh);
   Problem poisson(u, v);
   poisson = Integral((hmin + (hmax - hmin) * Pow(gamma, 3)) * Grad(u), Grad(v))
        - Integral(f * v)
        + DirichletBC(u, ScalarFunction(0.0)).on(GammaD);
   poisson.solve(solver);

   // Adjoint problem
   TrialFunction p(Vh);
   TestFunction  q(Vh);
   Problem adjoint(p, q);
   adjoint = Integral((hmin + (hmax - hmin) * Pow(gamma, 3)) * Grad(p), Grad(q))
        + Integral(ScalarFunction(1.0 / vol), q)
        + DirichletBC(p, ScalarFunction(0.0)).on(GammaD);
   adjoint.solve(solver);

   // Hilbert extension-regularization
   TrialFunction g(Vh);
   TestFunction  w(Vh);
   Problem hilbert(g, w);
   hilbert = Integral(alpha * Grad(g), Grad(w))
        + Integral(g, w)
        - Integral(
           ell + 3 * (hmax - hmin) * Pow(gamma, 2) *
            Dot(Grad(u.getSolution()), Grad(p.getSolution())), w)
        + DirichletBC(g, ScalarFunction(0.0)).on(GammaD);
   hilbert.solve(solver);

   GridFunction step(Ph);
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
