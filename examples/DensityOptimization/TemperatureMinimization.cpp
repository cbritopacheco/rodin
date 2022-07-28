#include "Rodin/Mesh.h"
#include "Rodin/Solver.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Variational;

// Define boundary attributes
static constexpr int Gamma0 = 1, GammaD = 2;

// Optimization parameters
static constexpr double ell = 1;
static constexpr double mu = 0.1;
static constexpr double min = 0.0001;
static constexpr double max = 1;
static constexpr double alpha = 0.1;
static constexpr size_t maxIterations = 25;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/temperature-minimization-example.mesh";

  // Load mesh
  Mesh Omega;
  Omega.load(meshFile);
  // Omega.refine();

  Omega.save("Omega.mesh");

  // Build finite element space
  FiniteElementSpace<H1> Vh(Omega);

  GridFunction gamma(Vh);
  gamma = 0.5;

  auto solver = Solver::UMFPack();

  // Optimization loop
  for (size_t i = 0; i < maxIterations; i++)
  {
    Alert::Info() << "Iteration: " << i << Alert::Raise;

    // Poisson problem
    ScalarFunction f(1.0);

    TrialFunction u(Vh);
    TestFunction  v(Vh);
    Problem poisson(u, v);
    poisson = Integral((min + (max - min) * Pow(gamma, 3)) * Grad(u), Grad(v))
            - Integral(f * v)
            + DirichletBC(u, ScalarFunction(0.0)).on(GammaD);
    solver.solve(poisson);
    u.getGridFunction().save("u.gf");

    // Adjoint problem
    TrialFunction p(Vh);
    TestFunction  q(Vh);
    Problem adjoint(p, q);
    adjoint = Integral((min + (max - min) * Pow(gamma, 3)) * Grad(p), Grad(q))
            + Integral(ScalarFunction(1.0 / Omega.getVolume()), q)
            + DirichletBC(p, ScalarFunction(0.0)).on(GammaD);
    solver.solve(adjoint);
    p.getGridFunction().save("p.gf");

    // Hilbert extension-regularization
    TrialFunction g(Vh);
    TestFunction  w(Vh);
    Problem hilbert(g, w);
    hilbert = Integral(alpha * Grad(g), Grad(w))
            + Integral(g, w)
            - Integral(
                ell + 3 * (max - min) * Pow(gamma, 2) *
                  Dot(Grad(u.getGridFunction()), Grad(p.getGridFunction())), w)
            + DirichletBC(g, ScalarFunction(0.0)).on(GammaD);
    solver.solve(hilbert);
    p.getGridFunction().save("g.gf");

    GridFunction step(Vh);
    step = mu * g.getGridFunction();
    gamma.getHandle() -= step.getHandle();

    gamma.save("gamma.gf");
  }

  return 0;
}
