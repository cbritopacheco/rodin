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
  const char* meshFile = "../resources/mfem/levelset-arch2d-example.mesh";

  // Define interior and exterior for level set discretization
  int Interior = 1, Exterior = 2;

  // Define boundary attributes
  int Gamma0 = 1,  // Traction free boundary
      GammaD = 2,  // Homogenous Dirichlet
      GammaN = 3,  // Inhomogenous Neumann
      Gamma  = 4;  // Shape boundary

  // Lamé coefficients
  auto mu     = ScalarCoefficient(0.3846),
       lambda = ScalarCoefficient(0.5769);

  // Compliance
  auto compliance = [&](GridFunction<H1>& w)
  {
    auto& Vh = w.getFiniteElementSpace();
    TrialFunction u(Vh);
    TestFunction  v(Vh);
    BilinearForm bf(v.getFiniteElementSpace());
    bf = Integral(lambda * Div(u), Div(v)).over(Interior)
       + Integral(
           mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T())).over(Interior);
    return bf(w, w);
  };

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);
  Omega.save("Omega0.mesh");
  std::cout << "Saved initial mesh to Omega0.mesh" << std::endl;

  // UMFPack
  auto solver = Solver::UMFPack();

  // Optimization parameters
  size_t maxIt = 250;
  double eps = 1e-6;
  double hmax = 0.05;
  auto ell = ScalarCoefficient(2.0);
  auto alpha = ScalarCoefficient(4 * hmax * hmax);

  std::vector<double> obj;

  // Optimization loop
  for (size_t i = 0; i < maxIt; i++)
  {
    // Vector field finite element space over the whole domain
    int d = 2;
    H1 Vh(Omega, d);

    // Trim the exterior part of the mesh to solve the elasticity system
    SubMesh trimmed = Omega.trim(Exterior, Gamma);

    // Build a finite element space over the trimmed mesh
    H1 VhInt(trimmed, d);

    // Elasticity equation
    auto f = VectorCoefficient{0, -1};
    TrialFunction uInt(VhInt);
    TestFunction  vInt(VhInt);
    Problem elasticity(uInt, vInt);
    elasticity = Integral(lambda * Div(uInt), Div(vInt))
               + Integral(
                   mu * (Jacobian(uInt) + Jacobian(uInt).T()), 0.5 * (Jacobian(vInt) + Jacobian(vInt).T()))
               - BoundaryIntegral(f, vInt).over(GammaN)
               + DirichletBC(uInt, VectorCoefficient{0, 0}).on(GammaD);
    solver.solve(elasticity);

    // Transfer solution back to original domain
    GridFunction u(Vh);
    uInt.getGridFunction().transfer(u);

    // Hilbert extension-regularization procedure
    auto e = 0.5 * (Jacobian(u).traceOf(Interior) + Jacobian(u).traceOf(Interior).T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    auto n = Normal(d);

    TrialFunction g(Vh);
    TestFunction  v(Vh);
    Problem hilbert(g, v);
    hilbert = Integral(alpha * Jacobian(g), Jacobian(v))
            + Integral(g, v)
            + BoundaryIntegral(Dot(Ae, e) - ell, Dot(v, n)).over(Gamma)
            + DirichletBC(g, VectorCoefficient{0, 0}).on(GammaN);
    solver.solve(hilbert);

    // Update objective
    obj.push_back(
        compliance(u) + ell.getValue() * Omega.getVolume(Interior));
    std::cout << "[" << i << "] Objective: " << obj.back() << std::endl;

    // Convert data types to mmg types
    auto mmgMesh = Cast(Omega).to<MMG::Mesh2D>();
    auto mmgVel = Cast(g.getGridFunction()).to<MMG::IncompleteVectorSolution2D>().setMesh(mmgMesh);

    // Generate signed distance function
    auto mmgLs = MMG::Distancer2D().setInteriorDomain(Interior).distance(mmgMesh);

    // Advect the level set function
    double gInf = std::max(g.getGridFunction().max(), -g.getGridFunction().min());
    double dt = 2 * hmax / gInf;
    MMG::Advect2D(mmgLs, mmgVel).step(dt);

    // Recover the implicit domain
    auto [mmgImplicit, _] =
      MMG::ImplicitDomainMesher2D().split(Interior, {Interior, Exterior})
                                   .split(Exterior, {Interior, Exterior})
                                   .setRMC(1e-3)
                                   .setHMax(hmax)
                                   .setBoundaryReference(Gamma)
                                   .discretize(mmgLs);

    // Convert back to Rodin data type
    Omega = Cast(mmgImplicit).to<Rodin::Mesh>();

    // Save mesh
    Omega.save("Omega.mesh");

    // Test for convergence
    if (obj.size() >= 2 && abs(obj[i] - obj[i - 1]) < eps)
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

