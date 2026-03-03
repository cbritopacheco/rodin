/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;


// Define interior and exterior for level set discretization
static constexpr Attribute Interior = 3, Exterior = 2;

// Define boundary attributes
static constexpr Attribute GammaD = 4, GammaN = 7, Gamma = 10, Gamma0 = 3;

// Lamé coefficients
static constexpr Real mu = 0.3846;
static constexpr Real lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 300;
static constexpr Real eps = 1e-12;
static constexpr Real hgrad = 1.6;
static constexpr Real ell = 0.1;
static Real elementStep = 0.5;
static Real hmax = 0.5;
static Real hmin = 0.1 * hmax;
static Real hausd = 0.2 * hmin;
static size_t hmaxIt = maxIt / 2;
const Real k = 0.5;
const Real dt = k * (hmax - hmin) / 2;
static Real alpha = dt;

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
};

int main(int, char**)
{
  const char* meshFile = "../resources/examples/ShapeOptimization/LevelSetCantilever3D.medit.mesh";

  // Load mesh
  MMG::Mesh th;
  th.load(meshFile, IO::FileFormat::MEDIT);

  P1 sh(th);
  GridFunction dist(sh);

  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;


  th.save("Omega0.mesh", IO::FileFormat::MEDIT);

  // Optimization loop
  std::vector<double> obj;
  std::ofstream fObj("obj.txt");
  for (size_t i = 0; i < maxIt; i++)
  {
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    Alert::Info() << "   | Optimizing the domain..." << Alert::Raise;
    MMG::Optimizer().setHMax(hmax)
                    .setHMin(hmin)
                    // .setGradation(1.1)
                    .setHausdorff(hausd)
                    .setAngleDetection(false)
                    .optimize(th);

    th.save("Optimized.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "   | Computing connectivity." << Alert::Raise;
    auto& conn = th.getConnectivity();

    conn.discover(3, 2);
    conn.discover(3, 1);

    conn.restrict(1, 0);
    conn.restrict(2, 0);

    conn.restrict(2, 3);

    conn.discover(0, 0);

    Alert::Info() << "   | Trimming mesh." << Alert::Raise;
    SubMesh trimmed = th.trim(Exterior);
    trimmed.save("Trimmed.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "   | Building finite element spaces." << Alert::Raise;
    const size_t d = th.getSpaceDimension();
    P1 sh(th);
    P1 vh(th, d);

    Alert::Info() << "   | Distancing domain." << Alert::Raise;

    GridFunction dist(sh);
    Distance::Eikonal(dist).setInterior(Interior)
                           .setInterface(Gamma)
                           .solve()
                           .sign();

    dist.getFiniteElementSpace().getMesh().save("Distance.mesh");
    dist.save("Distance.gf");

    P1 shInt(trimmed);
    P1 vhInt(trimmed, d);

    Alert::Info() << "   | Solving state equation." << Alert::Raise;
    auto f = VectorFunction{0, 0, -1};
    TrialFunction u(vhInt);
    TestFunction  v(vhInt);

    // Elasticity equation
    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - BoundaryIntegral(f, v).over(GammaN)
               + DirichletBC(u, VectorFunction{0, 0, 0}).on(GammaD);
    auto cg = Solver::CG(elasticity);
    cg.solve();

    u.getSolution().save("State.gf");
    u.getSolution().getFiniteElementSpace().getMesh().save("State.mesh");

    Alert::Info() << "   | Computing shape gradient." << Alert::Raise;
    auto jac = Jacobian(u.getSolution());
    jac.traceOf(Interior);
    auto e = 0.5 * (jac + jac.T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    auto n = FaceNormal(th);
    n.traceOf(Interior);

    // Hilbert extension-regularization procedure
    TrialFunction g(vh);
    TestFunction  w(vh);
    Problem hilbert(g, w);
    hilbert = Integral(alpha * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - FaceIntegral(Dot(Ae, e) - ell, Dot(n, w)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaN);
    Solver::CG(hilbert).solve();

    auto& dJ = g.getSolution();

    vh.getMesh().save("dJ.mesh");
    dJ.save("dJ.gf");

    // Update objective
    double objective = compliance(u.getSolution()) + ell * th.getVolume(Interior);
    obj.push_back(objective);
    fObj << objective << "\n";
    fObj.flush();
    Alert::Info() << "   | Objective: " << obj.back() << Alert::Raise;

    // Advect the level set function
    Alert::Info() << "   | Advecting the distance function." << Alert::Raise;
    GridFunction norm(sh);
    norm = Frobenius(dJ);
    dJ /= norm.max();

    TrialFunction advect(sh);
    TestFunction test(sh);

    Advection::Lagrangian(advect, test, dist, dJ).step(dt);

    th.save("Advect.mesh");
    advect.getSolution().save("Advect.gf");

    // Recover the implicit domain
    Alert::Info() << "   | Meshing the domain." << Alert::Raise;
    th = MMG::LevelSetDiscretizer().setHMax(hmax)
                                    .setHMin(hmin)
                                    .setHausdorff(hausd)
                                    .setAngleDetection(false)
                                    .setRMC(1e-5)
                                    .setBaseReferences(GammaD)
                                    .setBoundaryReference(Gamma)
                                    .discretize(advect.getSolution());

    th.save("Omega.mesh", IO::FileFormat::MEDIT);
  }

  Alert::Info() << "Saved final mesh to Omega.mesh" << Alert::Raise;

  return 0;
}

