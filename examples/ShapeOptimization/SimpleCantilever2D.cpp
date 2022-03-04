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
  const char* meshFile = "../resources/mfem/meshes/simple-cantilever-example.mesh";

  // Define boundary attributes
  int Gamma0 = 1, GammaD = 2, GammaN = 3;

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);
  Omega.save("Omega0.mesh");
  std::cout << "Saved initial mesh to Omega0.mesh" << std::endl;

  // Lamé coefficients
  auto mu     = ScalarCoefficient(0.3846),
       lambda = ScalarCoefficient(0.5769);

  // Conjugate Gradient solver
  auto cg = Solver::CG().setMaxIterations(500)
                        .setRelativeTolerance(1e-6);

  // Optimization parameters
  size_t maxIt = 30;
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

    // Compliance
    auto compliance = [&](GridFunction<H1>& v)
    {
      BilinearForm bf(Vh);
      bf = ElasticityIntegrator(lambda, mu);
      return bf(v, v);
    };

    // Elasticity equation
    TrialFunction u(Vh);
    TestFunction  v(Vh);
    Problem elasticity(u, v);
    elasticity = ElasticityIntegrator(lambda, mu)
               - BoundaryIntegral(VectorCoefficient{0, -1} * v).over(GammaN)
               + DirichletBC(u, VectorCoefficient{0, 0}).on(GammaD);
    cg.solve(elasticity);

    u.getGridFunction().save("u.gf");
    Omega.save("Omega.mesh");

    // Hilbert extension-regularization procedure
    TrialFunction g(Vh);
    TestFunction  w(Vh);
    auto e = Mult(0.5, Jacobian(u.getGridFunction()) + Jacobian(u.getGridFunction()).T());
    auto Ae = Mult(Mult(2.0, mu), e) + Mult(lambda, Mult(Trace(e), IdentityMatrix(d)));

    Problem hilbert(g, w);
    hilbert = VectorDiffusionIntegrator(alpha)
            + VectorMassIntegrator()
            - VectorBoundaryFluxLFIntegrator(Dot(Ae, e) - ell).over(Gamma0)
            + DirichletBC(g, VectorCoefficient{0, 0}).on(GammaD)
            + DirichletBC(g, VectorCoefficient{0, 0}).on(GammaN);
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
