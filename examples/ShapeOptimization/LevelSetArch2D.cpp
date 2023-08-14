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
using namespace Rodin::External;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// Define interior and exterior for level set discretization
static constexpr int Interior = 1, Exterior = 2;

// Define boundary attributes
static constexpr int Gamma0 = 1;  // Traction free boundary
static constexpr int GammaD = 2;  // Homogenous Dirichlet
static constexpr int GammaN = 3;  // Inhomogenous Neumann
static constexpr int Gamma  = 4;  // Shape boundary

// Lamé coefficients
static constexpr double mu    = 0.3846;
static constexpr double lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 250;
static constexpr double eps  = 1e-6;
static constexpr double hmax  = 0.05;
static constexpr double ell  = 0.2;
static constexpr double alpha = 4 * hmax * hmax;

// Compliance
double compliance(GridFunction<H1<Context::Serial>>& w);

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/levelset-arch2d-example.mesh";

  // Load mesh
  MMG::Mesh Omega;
  Omega.load(meshFile);

  MMG::MeshOptimizer().setHMax(hmax / 2.0).optimize(Omega);

  Omega.save("Omega0.mesh");

  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  Solver::UMFPack solver;

  // Optimization loop
  std::vector<double> obj;
  for (size_t i = 0; i < maxIt; i++)
  {
   Alert::Info() << "----- Iteration: " << i << Alert::Raise;

   // Vector field finite element space over the whole domain
   int d = 2;
   H1 Vh(Omega, d);

   // Trim the exterior part of the mesh to solve the elasticity system
   SubMesh trimmed = Omega.trim(Exterior);

   // Build a finite element space over the trimmed mesh
   H1 VhInt(trimmed, d);

   // Elasticity equation
   auto f = VectorFunction{0, -1};
   TrialFunction uInt(VhInt);
   TestFunction  vInt(VhInt);
   Problem elasticity(uInt, vInt);
   elasticity = Integral(lambda * Div(uInt), Div(vInt))
          + Integral(
             mu * (Jacobian(uInt) + Jacobian(uInt).T()), 0.5 * (Jacobian(vInt) + Jacobian(vInt).T()))
          - BoundaryIntegral(f, vInt).over(GammaN)
          + DirichletBC(uInt, VectorFunction{0, 0}).on(GammaD);
   elasticity.solve(solver);

   auto jac = Jacobian(uInt.getSolution());
   jac.traceOf(Interior);

   // Hilbert extension-regularization procedure
   auto e = 0.5 * (jac + jac.T());
   auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
   auto n = Normal(d);
   n.traceOf(Interior);

   TrialFunction g(Vh);
   TestFunction  v(Vh);
   Problem hilbert(g, v);
   hilbert = Integral(alpha * Jacobian(g), Jacobian(v))
        + Integral(g, v)
        - FaceIntegral(Dot(Ae, e) - ell, Dot(n, v)).over(Gamma)
        + DirichletBC(g, VectorFunction{0, 0}).on(GammaN);
   hilbert.solve(solver);
   g.getSolution().save("g.gf");
   Omega.save("g.mesh");

   // Update objective
   obj.push_back(
      compliance(uInt.getSolution()) + ell * Omega.getVolume(Interior));
   Alert::Info() << "   | Objective: " << obj[i] << Alert::Raise;

   // Generate signed distance function
   H1 Dh(Omega);
   auto dist = MMG::Distancer(Dh).setInteriorDomain(Interior)
                       .distance(Omega);

   // Advect the level set function
   Alert::Info() << "   | Advecting the distance function." << Alert::Raise;
   GridFunction gNorm(Dh);
   gNorm = ScalarFunction(
      [&](const Point& v) -> double
      {
       Math::Vector val = g.getSolution()(v);
       return val.norm();
      });
   double gInf = gNorm.max();
   double dt = 4 * hmax / gInf;
   MMG::Advect(dist, g.getSolution()).step(dt);

   // Recover the implicit domain
   Omega = MMG::ImplicitDomainMesher().split(Interior, {Interior, Exterior})
                          .split(Exterior, {Interior, Exterior})
                          .setRMC(1e-3)
                          .setAngleDetection(false)
                          .setBoundaryReference(Gamma)
                          .setBaseReferences(GammaD)
                          .discretize(dist);
   MMG::MeshOptimizer().setHMax(hmax).optimize(Omega);

   // Save mesh
   Omega.save("Omega.mesh");

   // Test for convergence
   if (obj.size() >= 2 && abs(obj[i] - obj[i - 1]) < eps)
   {
    Alert::Info() << "Convergence!" << Alert::Raise;
    break;
   }
  }

  Alert::Info() << "Saved final mesh to Omega.mesh" << Alert::Raise;

  return 0;
}

double compliance(GridFunction<H1<Context::Serial>>& w)
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

