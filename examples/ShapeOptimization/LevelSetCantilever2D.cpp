/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry/Region.h>
#include <Rodin/IO/ForwardDecls.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Variational/ForwardDecls.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/LinearElasticity.h>

#include <Rodin/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using FES = VectorP1<Mesh<Context::Local>>;

// Define interior and exterior for level set discretization
static constexpr Attribute interior = 1, exterior = 2;

// Define boundary attributes
static constexpr Attribute Gamma0 = 1, GammaD = 2, GammaN = 3, Gamma = 4;

// Lamé coefficients
static constexpr Real mu = 0.3846;
static constexpr Real lambda = 0.5769;

// Optimization parameters
static size_t maxIt = 300;
static Real hmax0 = 0.05;
static Real hmax = 0.05;
static Real hmin = 0.1 * hmax;
static Real hausd = 0.5 * hmin;
static Real ell = 0.4;
const Real dt = 0.5 * (hmax - hmin);
static Real alpha = 0.1;

// Compliance
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
  const char* meshFile = "../resources/examples/ShapeOptimization/LevelSetCantilever2D.mfem.mesh";

  // Load mesh
  MMG::Mesh th;
  th.load(meshFile, IO::FileFormat::MFEM);

  MMG::Optimizer().setHMax(hmax).setHMin(hmin).optimize(th);

  th.save("Omega0.mesh", IO::FileFormat::MEDIT);
  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  // XDMF output
  IO::XDMF xdmf("out/ShapeOptimization");

  auto domain = xdmf.grid("domain");
  domain.setMesh(th, IO::XDMF::MeshPolicy::Transient);

  auto state  = xdmf.grid("state");

  // Optimization loop
  std::vector<double> obj;
  std::ofstream fObj("obj.txt");
  size_t i = 0;
  while (i < maxIt)
  {
    try
    {
      MMG::Optimizer().setHMax(hmax)
                      .setHMin(hmin)
                      .setHausdorff(hausd)
                      .setAngleDetection(false)
                      .optimize(th);

      hmax = hmax0;
      hmin = 0.1 * hmax;
    }
    catch (const Alert::Exception& e)
    {
      hmax /= 2;
      hmin = 0.1 * hmax;
      Alert::Warning() << "Meshing failed at iteration " << i
        << ". Reducing hmax to " << hmax << " and retrying." << Alert::Raise;
      continue;
    }

    th.getConnectivity().compute(1, 2);
    th.getConnectivity().compute(0, 0);

    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    Alert::Info() << "   | Trimming mesh." << Alert::Raise;
    SubMesh trimmed = th.trim(exterior);
    trimmed.save("Omega.mesh", IO::FileFormat::MFEM);

    Alert::Info() << "   | Building finite element spaces." << Alert::Raise;
    const size_t d = th.getSpaceDimension();
    P1 sh(th);
    P1 vh(th, d);

    P1 shInt(trimmed);
    P1 vhInt(trimmed, d);

    Alert::Info() << "   | Solving state equation." << Alert::Raise;
    auto f = VectorFunction{0, -1};
    TrialFunction u(vhInt);
    TestFunction  v(vhInt);

    // Elasticity equation
    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - BoundaryIntegral(f, v).over(GammaN)
               + DirichletBC(u, VectorFunction{0, 0}).on(GammaD);
    Solver::CG(elasticity).solve();

    u.getSolution().save("State.gf", IO::FileFormat::MFEM);
    trimmed.save("State.mesh", IO::FileFormat::MFEM);

    Alert::Info() << "   | Computing shape gradient." << Alert::Raise;
    auto jac = Jacobian(u.getSolution());
    jac.traceOf(interior);
    auto e = 0.5 * (jac + jac.T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    auto n = FaceNormal(th);
    n.traceOf(interior);

    // Hilbert extension-regularization procedure
    TrialFunction g(vh);
    TestFunction w(vh);
    Problem hilbert(g, w);
    hilbert = Integral(alpha * alpha * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - FaceIntegral(Dot(Ae, e) - ell, Dot(n, w)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaN);
    Solver::CG(hilbert).solve();

    auto& dJ = g.getSolution();
    dJ.save("dJ.gf", IO::FileFormat::MFEM);
    vh.getMesh().save("dJ.mesh", IO::FileFormat::MFEM);

    // Update objective
    double objective = compliance(u.getSolution()) + ell * th.getArea(interior);
    obj.push_back(objective);
    fObj << objective << "\n";
    fObj.flush();
    Alert::Info() << "   | Objective: " << obj.back() << Alert::Raise;
    Alert::Info() << "   | Distancing domain." << Alert::Raise;

    GridFunction dist(sh);
    Distance::Eikonal(dist).setInterior(interior)
                           .setInterface(Gamma)
                           .solve()
                           .sign();

    th.save("Distance.mesh", IO::FileFormat::MFEM);
    dist.save("Distance.gf", IO::FileFormat::MFEM);

    // Advect the level set function
    Alert::Info() << "   | Advecting the distance function." << Alert::Raise;
    GridFunction norm(sh);
    norm = Frobenius(dJ);
    dJ /= norm.max();

    TrialFunction advect(sh);
    TestFunction test(sh);

    Advection::Lagrangian(advect, test, dist, dJ).step(dt);

    th.save("Advect.mesh", IO::FileFormat::MFEM);
    advect.getSolution().save("Advect.gf", IO::FileFormat::MFEM);

    // ------------------------------------------------------------------------
    // XDMF snapshot for the current iteration
    // ------------------------------------------------------------------------

    // Full evolving mesh 'th'
    domain.clear();
    domain.add("dJ", dJ, IO::XDMF::Center::Node);
    domain.add("dist", dist, IO::XDMF::Center::Node);
    domain.add("advect", advect.getSolution(), IO::XDMF::Center::Node);

    // Trimmed mesh and state field
    state.clear();
    state.setMesh(trimmed, IO::XDMF::MeshPolicy::Transient);
    state.add("u", u.getSolution(), IO::XDMF::Center::Node);

    xdmf.write(i).flush();

    // Recover the implicit domain
    Alert::Info() << "   | Meshing the domain." << Alert::Raise;

    try
    {
      th = MMG::LevelSetDiscretizer().split(interior, {interior, exterior})
                                     .split(exterior, {interior, exterior})
                                     .setRMC(1e-6)
                                     .setHMax(hmax)
                                     .setHMin(hmin)
                                     .setHausdorff(hausd)
                                     .setAngleDetection(false)
                                     .setBoundaryReference(Gamma)
                                     .setBaseReferences(GammaD)
                                     .discretize(advect.getSolution());

      hmax = hmax0;
      hmin = 0.1 * hmax;
    }
    catch (const Alert::Exception& e)
    {
      hmax /= 2;
      hmin = 0.1 * hmax;
      Alert::Warning() << "Meshing failed at iteration " << i
        << ". Reducing hmax to " << hmax << " and retrying." << Alert::Raise;
      continue;
    }

    i++;
    th.save("out/Omega." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
  }

  xdmf.close();

  Alert::Success() << "Saved final mesh to Omega.mesh" << Alert::Raise;

  return 0;
}
