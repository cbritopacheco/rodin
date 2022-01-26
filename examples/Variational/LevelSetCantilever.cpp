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
  const char* meshFile = "Omega.mesh";

  // Define interior and exterior for level set discretization
  int Interior = 1, Exterior = 2;

  // Define boundary attributes
  int Gamma0 = 1, GammaD = 2, GammaN = 3, Gamma = 4;

  // Load mesh
  Mesh D = Mesh::load(meshFile);
  D.save("Omega0.mesh");
  std::cout << "Saved initial mesh to Omega0.mesh" << std::endl;

  // Lamé coefficients
  auto mu     = ScalarCoefficient(0.3846),
       lambda = ScalarCoefficient(0.5769);

  // Conjugate Gradient solver
  auto cg = Solver::CG().setMaxIterations(500)
                        .setRelativeTolerance(1e-12);

  // Optimization parameters
  size_t maxIt = 30;
  double eps = 1e-6;
  double oldObj, newObj;
  auto ell = ScalarCoefficient(5);
  auto alpha = ScalarCoefficient(0.1);

  // Optimization loop
  for (size_t i = 0; i < maxIt; i++)
  {
    // Finite element space
    int d = 2;
    H1 Vh(D, d);

    // Compliance
    auto compliance = [&](GridFunction<H1>& v)
    {
      BilinearForm bf(Vh);
      bf = ElasticityIntegrator(lambda, mu).over(Interior);
      return bf(v, v);
    };

    // Elasticity equation
    GridFunction u(Vh);
    Problem elasticity(u);
    elasticity = ElasticityIntegrator(lambda, mu).over(Interior)
               + DirichletBC(GammaD, VectorCoefficient{0, 0})
               + NeumannBC(GammaN, VectorCoefficient{0, -1});
    cg.solve(elasticity);

    // Hilbert extension-regularization procedure
    GridFunction theta(Vh);
    auto e = ScalarCoefficient(0.5) * (Jacobian(u) + Jacobian(u).T());
    auto Ae = ScalarCoefficient(2.0) * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    Problem hilbert(theta);
    hilbert = VectorDiffusionIntegrator(alpha)
            + VectorMassIntegrator()
            - VectorBoundaryFluxLFIntegrator(Dot(Ae, e) - ell).over(Gamma)
            + DirichletBC(GammaN, VectorCoefficient{0, 0});
    cg.solve(hilbert);

    // Update objective
    oldObj = newObj;
    newObj = compliance(u) + ell.getValue() * D.getVolume(Interior);

    std::cout << "[" << i << "] Objective: " << newObj << std::endl;

    // Test for convergence
    if (i > 0 && abs(oldObj - newObj) < eps)
      break;

    // Convert data types to mmg types
    auto mmgMesh = Cast(D).to<MMG::Mesh2D>();
    auto mmgDisp = Cast(theta).to<MMG::IncompleteVectorSolution2D>().setMesh(mmgMesh);

    // Generate signed distance function
    auto ls =
      MMG::Distancer2D().setInteriorDomains({ Interior }).distance(mmgMesh);

    // Advect the level set function
    MMG::Advect2D(ls, mmgDisp).step(0.01);

    // Recover the implicit domain
    auto [mmgImplicit, _] =
      MMG::ImplicitDomainMesher2D().split(Interior, {Interior, Exterior})
                                   .setBoundaryReference(Gamma)
                                   .discretize(ls);

    // Convert back to Rodin data type
    D = Cast(mmgImplicit).to<Rodin::Mesh>();

    // Save current mesh
    D.save("Omega.mesh");
  }

  std::cout << "Saved final mesh to Omega.mesh" << std::endl;

  return 0;
}

