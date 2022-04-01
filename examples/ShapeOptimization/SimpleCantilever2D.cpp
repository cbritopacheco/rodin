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

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/simple-cantilever2d-example.mesh";

  // Define boundary attributes
  int Gamma0 = 1, GammaD = 2, GammaN = 3;

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);
  Omega.save("Omega0.mesh");
  std::cout << "Saved initial mesh to Omega0.mesh" << std::endl;

  // Lamé coefficients
  auto mu     = ScalarCoefficient(0.3846),
       lambda = ScalarCoefficient(0.5769);

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

  // Optimization parameters
  size_t maxIt = 40;
  double eps = 1e-6;
  double hmax = 0.1;
  auto ell = ScalarCoefficient(5);
  auto alpha = ScalarCoefficient(4 * hmax * hmax);

  std::vector<double> obj;

  // Optimization loop
  for (size_t i = 0; i < maxIt; i++)
  {
    // Finite element spaces
    int d = 2;
    H1 Vh(Omega, d);

    // Pull-down force
    auto f = VectorCoefficient{0, -1};

    // Elasticity equation
    TrialFunction u(Vh);
    TestFunction  v(Vh);
    Problem elasticity(u, v);
    elasticity = Integral(lambda * Div(u), Div(v))
               + Integral(
                   mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()))
               - BoundaryIntegral(f, v).over(GammaN)
               + DirichletBC(u, VectorCoefficient{0, 0}).on(GammaD);
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
            + DirichletBC(g, VectorCoefficient{0, 0}).on({GammaD, GammaN});
    cg.solve(hilbert);

    // Update objective
    obj.push_back(compliance(u.getGridFunction()) + ell.getValue() * Omega.getVolume());

    std::cout << "[" << i << "] Objective: " << obj[i] << std::endl;

    // Test for convergence
    if (i > 0 && abs(obj[i] - obj[i - 1]) < eps)
      break;

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
      auto mmgMesh = Cast(Omega).to<MMG::Mesh2D>();
      MMG::MeshOptimizer2D().setHMax(hmax).optimize(mmgMesh);
      Omega = Cast(mmgMesh).to<Mesh>();
    }

    // Save mesh
    Omega.save("Omega.mesh");
  }

  std::cout << "Saved final mesh to Omega.mesh" << std::endl;

  std::ofstream plt("obj.txt");
  for (size_t i = 0; i < obj.size(); i++)
    plt << obj[i] << "\n";

  std::cout << "Saved objective history to obj.txt" << std::endl;

  return 0;
}
