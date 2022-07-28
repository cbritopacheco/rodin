/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Mesh.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Variational;

// Define boundary attributes
static constexpr int Gamma0 = 1; // Traction free boundary
static constexpr int GammaD = 2; // Homogenous Dirichlet boundary
static constexpr int GammaN = 3; // Inhomogenous Neumann boundary

// Lamé coefficients
static constexpr double mu     = 0.3846;
static constexpr double lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 40;
static constexpr double eps   = 1e-6;
static constexpr double hmax  = 0.1;
static constexpr double ell   = 5.0;
static constexpr double alpha = 4 * hmax * hmax;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/simple-cantilever2d-example.mesh";

  // Load mesh
  Mesh Omega;
  Omega.load(meshFile);
  Omega.save("Omega0.mesh");
  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  // Compliance
  auto compliance = [&](GridFunction<H1>& w)
  {
    auto& Vh = w.getFiniteElementSpace();
    TrialFunction u(Vh);
    TestFunction  v(Vh);
    BilinearForm  bf(u, v);
    bf = Integral(lambda * Div(u), Div(v))
       + Integral(
           mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()));
    return bf(w, w);
  };

  // Conjugate Gradient solver
  auto cg = Solver::CG().setMaxIterations(500)
                        .setRelativeTolerance(1e-6);

  // Optimization loop
  std::vector<double> obj;
  for (size_t i = 0; i < maxIt; i++)
  {
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    // Finite element spaces
    int d = 2;
    FiniteElementSpace<H1> Vh(Omega, d);

    // Pull-down force
    auto f = VectorFunction{0, -1};

    // Elasticity equation
    TrialFunction u(Vh);
    TestFunction  v(Vh);
    Problem elasticity(u, v);
    elasticity = Integral(lambda * Div(u), Div(v))
               + Integral(
                   mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()))
               - BoundaryIntegral(f, v).over(GammaN)
               + DirichletBC(u, VectorFunction{0, 0}).on(GammaD);
    cg.solve(elasticity);

    // Hilbert extension-regularization procedure
    TrialFunction g(Vh);
    TestFunction  w(Vh);

    auto e = 0.5 * (Jacobian(u.getGridFunction()) + Jacobian(u.getGridFunction()).T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    auto n = Normal(d);

    Problem hilbert(g, w);
    hilbert = Integral(alpha * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - BoundaryIntegral(Dot(Ae, e) - ell, Dot(w, n)).over(Gamma0)
            + DirichletBC(g, VectorFunction{0, 0}).on({GammaD, GammaN});
    cg.solve(hilbert);

    // Update objective
    obj.push_back(compliance(u.getGridFunction()) + ell * Omega.getVolume());

    Alert::Info() << "    | Objective: " << obj[i] << Alert::Raise;

    // Test for convergence
    if (i > 0 && abs(obj[i] - obj[i - 1]) < eps)
    {
      Alert::Info() << "Convergence!" << Alert::Raise;
      break;
    }

    // Make the displacement
    double dt = Omega.getMaximumDisplacement(g.getGridFunction());
    if (dt < 1e-4)
    {
      // If the maximum displacement is too small the mesh will degenerate
      break;
    }
    else
    {
      g.getGridFunction() *= hmax * dt;
      Omega.displace(g.getGridFunction());

      // Refine the mesh using MMG
      MMG::MeshOptimizer().setHMax(hmax).optimize(Omega);
    }

    // Save mesh
    Omega.save("Omega.mesh");
  }

  Alert::Info() << "Saved final mesh to Omega.mesh" << Alert::Raise;

  return 0;
}
