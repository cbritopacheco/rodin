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
using namespace Rodin::External;
using namespace Rodin::Variational;


// Define interior and exterior for level set discretization
static constexpr int Interior = 1, Exterior = 2;

// Define boundary attributes
static constexpr int Gamma0 = 1, GammaD = 2, GammaN = 3, Gamma = 4;

// Lamé coefficients
static constexpr double mu = 0.3846;
static constexpr double lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 250;
static constexpr double eps = 1e-6;
static constexpr double hmax = 0.05;
static constexpr double ell = 0.4;
static constexpr double alpha = 4 * hmax * hmax;

// Compliance
double compliance(GridFunction<H1<Context::Serial>>& w);

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/levelset-cantilever2d-example.mesh";

  // Load mesh
  Mesh Omega;
  Omega.load(meshFile);

  Omega.save("Omega0.mesh");
  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  // UMFPack
  auto solver = Solver::UMFPack();

  // Optimization loop
  std::vector<double> obj;
  std::ofstream fObj("obj.txt");
  for (size_t i = 0; i < maxIt; i++)
  {
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    Alert::Info() << "    | Trimming mesh." << Alert::Raise;
    SubMesh trimmed = Omega.trim(Exterior);

    Alert::Info() << "    | Building finite element spaces." << Alert::Raise;
    int d = 2;
    H1 Vh(Omega, d);
    H1 VhInt(trimmed, d);

    Alert::Info() << "    | Solving state equation." << Alert::Raise;
    auto f = VectorFunction{0, -1};
    TrialFunction uInt(VhInt);
    TestFunction  vInt(VhInt);

    // Elasticity equation
    Problem elasticity(uInt, vInt);
    elasticity = Integral(lambda * Div(uInt), Div(vInt))
               + Integral(
                   mu * (Jacobian(uInt) + Jacobian(uInt).T()), 0.5 * (Jacobian(vInt) + Jacobian(vInt).T()))
               - BoundaryIntegral(f, vInt).over(GammaN)
               + DirichletBC(uInt, VectorFunction{0, 0}).on(GammaD);
    solver.solve(elasticity);

    // Transfer solution back to original domain
    Alert::Info() << "    | Computing shape gradient." << Alert::Raise;
    auto e = 0.5 * (Jacobian(uInt.getGridFunction()) + Jacobian(uInt.getGridFunction()).T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    auto n = Normal(d);

    GridFunction normal(Vh);
    normal.projectOnBoundary(n, Gamma);
    Omega.save("normal.mesh");
    normal.save("normal.gf");

    // Hilbert extension-regularization procedure
    TrialFunction g(Vh);
    TestFunction  v(Vh);
    Problem hilbert(g, v);
    hilbert = Integral(alpha * Jacobian(g), Jacobian(v))
            + Integral(g, v)
            - BoundaryIntegral(Dot(Ae, e) - ell, Dot(n, v)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0}).on(GammaN);
    solver.solve(hilbert);

    trimmed.save("out/trimmed." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);

    // Update objective
    double objective = compliance(uInt.getGridFunction()) + ell * Omega.getVolume(Interior);
    obj.push_back(objective);
    fObj << objective << "\n";
    fObj.flush();
    Alert::Info() << "    | Objective: " << obj.back() << Alert::Raise;

    Alert::Info() << "    | Distancing domain." << Alert::Raise;
    H1 Dh(Omega);
    auto dist = MMG::Distancer(Dh).setInteriorDomain(Interior)
                                  .distance(Omega);
    Omega.save("dist.mesh");
    g.getGridFunction().save("g.gf");
    dist.save("dist.gf");

    // Advect the level set function
    Alert::Info() << "    | Advecting the distance function." << Alert::Raise;
    GridFunction gNorm(Dh);
    gNorm = ScalarFunction(
        [&](const Vertex& v) -> double
        {
          mfem::Vector val = g.getGridFunction()(v);
          return val.Norml2();
        });
    double gInf = gNorm.max();
    double dt = 4 * hmax / gInf;
    MMG::Advect(dist, g.getGridFunction()).step(dt);

    // Recover the implicit domain
    Alert::Info() << "    | Meshing the domain." << Alert::Raise;
    Omega = MMG::ImplicitDomainMesher().split(Interior, {Interior, Exterior})
                                       .split(Exterior, {Interior, Exterior})
                                       .setRMC(1e-3)
                                       .setBoundaryReference(Gamma)
                                       .discretize(dist);
    MMG::MeshOptimizer().setHMax(hmax).optimize(Omega);

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

