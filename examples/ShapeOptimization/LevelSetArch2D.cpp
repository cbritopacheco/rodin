/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Geometry/Region.h"
#include "Rodin/Models/Advection/Lagrangian.h"
#include "Rodin/Models/Distance/Eikonal.h"
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/LinearElasticity.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::External;
using namespace Rodin::Variational;

// Define interior and exterior for level set discretization
static constexpr Attribute interior = 1, exterior = 2;

// Define boundary attributes
static constexpr Attribute Gamma0 = 1, GammaD = 2, GammaN = 3, Gamma = 4;

// Lamé coefficients
static constexpr double mu = 0.3846;
static constexpr double lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 250;
static constexpr double eps  = 1e-6;
static constexpr double hmax  = 0.05;
static constexpr double hmin = 0.1 * hmax;
static constexpr double hausd = 0.5 * hmin;
static constexpr double ell  = 0.2;
const constexpr Real dt = 0.5 * (hmax - hmin);
static constexpr double alpha = 0.1;

using FES = VectorP1<Mesh<Context::Local>>;

template <class Data>
Real compliance(const GridFunction<FES, Data>& w)
{
  auto& vh = w.getFiniteElementSpace();
  TrialFunction u(vh);
  TestFunction  v(vh);
  BilinearForm  bf(u, v);
  bf = LinearElasticityIntegral(u, v)(lambda, mu);
  bf.assemble();
  return bf(w, w);
}

int main(int, char**)
{
  const char* meshFile = "../resources/examples/ShapeOptimization/LevelSetArch2D.mfem.mesh";

  // Load mesh
  MMG::Mesh th;
  th.load(meshFile);
  MMG::Optimizer().setHMax(hmax).setHMin(hmin).optimize(th);

  th.save("Omega0.mesh", IO::FileFormat::MEDIT);

  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  // Optimization loop
  std::vector<double> obj;
  std::ofstream fObj("obj.txt");
  for (size_t i = 0; i < maxIt; i++)
  {
    th.getConnectivity().compute(1, 2);
    th.getConnectivity().compute(0, 0);
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    SubMesh trimmed = th.trim(exterior);
    trimmed.save("Omega.mesh");

    const size_t d = th.getSpaceDimension();
    P1 sh(th);
    P1 vh(th, d);
    P1 vhInt(trimmed, d);

    auto f = VectorFunction{0, -1};
    TrialFunction u(vhInt);
    TestFunction  v(vhInt);
    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - BoundaryIntegral(f, v).over(GammaN)
               + DirichletBC(u, VectorFunction{0, 0}).on(GammaD);
    Solver::CG(elasticity).solve();

    auto jac = Jacobian(u.getSolution());
    jac.traceOf(interior);
    auto e = 0.5 * (jac + jac.T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    auto n = FaceNormal(th);
    n.traceOf(interior);

    TrialFunction g(vh);
    TestFunction w(vh);
    Problem hilbert(g, w);
    hilbert = Integral(alpha * alpha * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - FaceIntegral(Dot(Ae, e) - ell, Dot(n, w)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaN);
    Solver::CG(hilbert).solve();

    auto& dJ = g.getSolution();
    const double objective = compliance(u.getSolution()) + ell * th.getArea(interior);
    obj.push_back(objective);
    fObj << objective << "\n";
    fObj.flush();
    Alert::Info() << "   | Objective: " << objective << Alert::Raise;

    GridFunction dist(sh);
    Models::Distance::Eikonal(dist).setInterior(interior)
                                   .setInterface(Gamma)
                                   .solve()
                                   .sign();

    GridFunction norm(sh);
    norm = Frobenius(dJ);
    dJ /= norm.max();

    TrialFunction advect(sh);
    TestFunction test(sh);
    Models::Advection::Lagrangian(advect, test, dist, dJ).step(dt);

    th = MMG::ImplicitDomainMesher().split(interior, {interior, exterior})
                                    .split(exterior, {interior, exterior})
                                    .setRMC(1e-6)
                                    .setHMax(hmax)
                                    .setHMin(hmin)
                                    .setHausdorff(hausd)
                                    .setAngleDetection(false)
                                    .setBoundaryReference(Gamma)
                                    .setBaseReferences(GammaD)
                                    .discretize(advect.getSolution());

    MMG::Optimizer().setHMax(hmax)
                    .setHMin(hmin)
                    .setHausdorff(hausd)
                    .setAngleDetection(false)
                    .optimize(th);

    th.save("out/Omega." + std::to_string(i) + ".mesh");

    if (obj.size() >= 2 && std::abs(obj[i] - obj[i - 1]) < eps)
    {
      Alert::Info() << "Convergence!" << Alert::Raise;
      break;
    }
  }

  Alert::Info() << "Saved final meshes to out/Omega.*.mesh" << Alert::Raise;

  return 0;
}
