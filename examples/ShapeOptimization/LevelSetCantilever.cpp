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
  // const char* meshFile = "../resources/mfem/meshes/levelset-cantilever-example.mesh";
  const char* meshFile = "Omega.mesh";

  // Define interior and exterior for level set discretization
  int Interior = 1, Exterior = 2;

  // Define boundary attributes
  int Gamma0 = 1, GammaD = 2, GammaN = 3, Gamma = 4;

  // Lamé coefficients
  auto mu     = ScalarCoefficient(0.3846),
       lambda = ScalarCoefficient(0.5769);

  // Compliance
  auto compliance = [&](GridFunction<H1>& v)
  {
    BilinearForm bf(v.getFiniteElementSpace());
    bf = ElasticityIntegrator(lambda, mu).over(Interior);
    // Set to zero where the function might not defined
    v[Rodin::isNaN(v) || Rodin::isInf(v)] = 0.0;
    return bf(v, v);
  };

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);
  Omega.save("Omega0.mesh");
  std::cout << "Saved initial mesh to Omega0.mesh" << std::endl;

  // UMFPack
  auto solver = Solver::UMFPack();

  // Optimization parameters
  size_t maxIt = 200;
  double eps = 1e-6;
  double hmax = 0.05;
  auto ell = ScalarCoefficient(0.5);
  auto alpha = ScalarCoefficient(hmax * hmax);

  std::vector<double> obj;

  // Optimization loop
  for (size_t i = 0; i < maxIt; i++)
  {
    // Finite element space
    int d = 2;
    H1 Vh(Omega, d);

    // Elasticity equation
    GridFunction u(Vh);
    Problem elasticity(u);
    elasticity = ElasticityIntegrator(lambda, mu).over(Interior)
               + DirichletBC(GammaD, VectorCoefficient{0, 0})
               + NeumannBC(GammaN, VectorCoefficient{0, -1});
    solver.solve(elasticity);

    // Hilbert extension-regularization procedure
    GridFunction theta(Vh);
    auto e = ScalarCoefficient(0.5) * (Jacobian(u) + Jacobian(u).T());
    auto Ae = ScalarCoefficient(2.0) * mu * e + lambda * Trace(e) * IdentityMatrix(d);

    H1 Ph(Omega);
    GridFunction w(Ph);
    w = (Dot(Ae, e) - ell).restrictedTo(Interior);

    Problem hilbert(theta);
    hilbert = VectorDiffusionIntegrator(alpha)
            + VectorMassIntegrator()
            - DivergenceIntegrator<Domain, Linear>(Dot(Ae, e) - ell).over(Interior)
            - VectorDomainLFIntegrator(Gradient(w)).over(Interior)
            + DirichletBC(GammaN, VectorCoefficient{0, 0});
    solver.solve(hilbert);

    Omega.save("Omegai.mesh");
    theta.save("theta.gf");

    if (i == 0)
    {
      obj.push_back(
          compliance(u) + ell.getValue() * Omega.getVolume(Interior));
    }

    // Linear search
    int maxSearchIt = 5;
    int searchIt = 0;

    double min = theta.min(),
           max = theta.max();
    double linf = abs(min) > max ? abs(min) : max;
    double dt = hmax / linf;
    while (searchIt++ < maxSearchIt)
    {
      double newObj = obj.back(), oldObj;
      Mesh OmegaCandidate(Omega);

      // Convert data types to mmg types
      auto mmgMesh = Cast(OmegaCandidate).to<MMG::Mesh2D>();
      auto mmgVel = Cast(theta).to<MMG::IncompleteVectorSolution2D>().setMesh(mmgMesh);

      // Generate signed distance function
      auto mmgLs =
        MMG::Distancer2D().setInteriorDomains({ Interior }).distance(mmgMesh);

      // Advect the level set function
      MMG::Advect2D advect(mmgLs, mmgVel);
      advect.enableExtrapolation(false).step(dt);

      // Recover the implicit domain
      auto [mmgImplicit, _] =
        MMG::ImplicitDomainMesher2D().split(Interior, {Interior, Exterior})
                                     .split(Exterior, {Interior, Exterior})
                                     .setRMC(1e-3)
                                     .setHMax(hmax)
                                     .setBoundaryReference(Gamma)
                                     .discretize(mmgLs);

      // Convert back to Rodin data type
      OmegaCandidate = Cast(mmgImplicit).to<Rodin::Mesh>();

      // Compute new objective
      H1 Ch(OmegaCandidate, d);
      GridFunction uc(Ch);
      Problem elasticityc(uc);
      elasticityc = ElasticityIntegrator(lambda, mu).over(Interior)
                  + DirichletBC(GammaD, VectorCoefficient{0, 0})
                  + NeumannBC(GammaN, VectorCoefficient{0, -1});
      solver.solve(elasticityc);

      oldObj = newObj;
      newObj = compliance(uc) + ell.getValue() * OmegaCandidate.getVolume(Interior);

      if (oldObj < newObj)
      {
        // Reject iteration
        dt /= 2.0;
      } else
      {
        // Accept iteration
        Omega = std::move(OmegaCandidate);

        // Update objective
        obj.push_back(newObj);
        std::cout << "[" << i << "] Objective: " << obj.back() << std::endl;
        std::cout << "dt = " << dt << std::endl;
        break;
      }
    }

    Omega.save("Omega.mesh");

    // Test for convergence
    if ((i > 0 && abs(obj[i] - obj[i - 1]) < eps))
    {
      std::cout << "Convergence!" << std::endl;
      break;
    }

    std::ofstream plt("obj.txt", std::ios::trunc);
    for (size_t i = 0; i < obj.size(); i++)
      plt << i << "," << obj[i] << "\n";
  }

  std::cout << "Saved final mesh to Omega.mesh" << std::endl;

  return 0;
}

