/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/LinearElasticity.h>>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::External;
using namespace Rodin::Variational;

using FES = VectorP1<Context::Serial>;

// Define interior and exterior for level set discretization
static constexpr Attribute Interior = 1, Exterior = 2;

// Define boundary attributes
static constexpr Attribute Gamma0 = 1, GammaD = 2, GammaN = 3, Gamma = 4;

// Lamé coefficients
static constexpr double mu = 0.3846;
static constexpr double lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 300;
static constexpr double eps = 1e-6;
static constexpr double hmax = 0.05;
static constexpr double ell = 0.4;
static constexpr double alpha = 4 * hmax * hmax;

// Compliance
inline Scalar compliance(const GridFunction<FES>& w)
{
  auto& vh = w.getFiniteElementSpace();
  TrialFunction u(vh);
  TestFunction  v(vh);
  BilinearForm  bf(u, v);
  bf = Integral(lambda * Div(u), Div(v))
     + Integral(mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()));
  return bf(w, w);
};

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/levelset-cantilever2d-example.mesh";

  // Load mesh
  MMG::Mesh Omega;
  Omega.load(meshFile);

  MMG::Optimizer().setHMax(hmax / 2.0).optimize(Omega);

  Omega.save("Omega0.mesh");
  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  // Solver
  Solver::SparseLU solver;

  // Optimization loop
  std::vector<double> obj;
  std::ofstream fObj("obj.txt");
  for (size_t i = 0; i < maxIt; i++)
  {
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    Alert::Info() << "   | Trimming mesh." << Alert::Raise;
    SubMesh trimmed = Omega.trim(Exterior);
    trimmed.save("trimmed.mesh");
    std::exit(1);

    Alert::Info() << "   | Building finite element spaces." << Alert::Raise;
    int d = 2;
    P1 vh(Omega, d);
    P1 shInt(trimmed);
    P1 vhInt(trimmed, d);
    Omega.save("miaow.mesh");

    Alert::Info() << "   | Solving state equation." << Alert::Raise;
    auto f = VectorFunction{0, -1};
    TrialFunction uInt(vhInt);
    TestFunction  vInt(vhInt);

    // Elasticity equation
    Problem elasticity(uInt, vInt);
    elasticity = LinearElasticityIntegral(uInt, vInt)(lambda, mu)
               - BoundaryIntegral(f, vInt).over(GammaN)
               + DirichletBC(uInt, VectorFunction{0, 0}).on(GammaD);
    elasticity.solve(solver);

    Alert::Info() << "   | Computing shape gradient." << Alert::Raise;

    auto jac = Jacobian(uInt.getSolution());
    jac.traceOf(Interior);

    auto e = 0.5 * (jac + jac.T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    auto n = FaceNormal(Omega);
    n.traceOf(Interior);

    // Hilbert extension-regularization procedure
    TrialFunction g(vh);
    TestFunction  v(vh);
    Problem hilbert(g, v);
    hilbert = Integral(alpha * Jacobian(g), Jacobian(v))
            + Integral(g, v)
            - FaceIntegral(Dot(Ae, e) - ell, Dot(n, v)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0}).on(GammaN);
    hilbert.solve(solver);

    // Update objective
    double objective = compliance(uInt.getSolution()) + ell * Omega.getVolume(Interior);
    obj.push_back(objective);
    fObj << objective << "\n";
    fObj.flush();
    Alert::Info() << "   | Objective: " << obj.back() << Alert::Raise;
    Alert::Info() << "   | Distancing domain." << Alert::Raise;
    P1 Dh(Omega);
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
    Alert::Info() << "   | Meshing the domain." << Alert::Raise;

    Omega = MMG::ImplicitDomainMesher().split(Interior, {Interior, Exterior})
                                       .split(Exterior, {Interior, Exterior})
                                       .setRMC(1e-3)
                                       .setAngleDetection(false)
                                       .setBoundaryReference(Gamma)
                                       .setBaseReferences(GammaD)
                                       .discretize(dist);

    MMG::Optimizer().setHMax(hmax).optimize(Omega);

    Omega.save("Omega.mesh");
  }

  Alert::Info() << "Saved final mesh to Omega.mesh" << Alert::Raise;

  return 0;
}

