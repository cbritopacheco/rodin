/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/LinearElasticity.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::External;
using namespace Rodin::Variational;

using FES = VectorP1<Context::Sequential>;

// Define interior and exterior for level set discretization
static constexpr Attribute Interior = 1, Exterior = 2;

// Define boundary attributes
static constexpr Attribute Gamma0 = 1, GammaD = 2, GammaN = 3, Gamma = 4;

// Lamé coefficients
static constexpr double mu = 0.3846;
static constexpr double lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 300;
static constexpr double hmax = 0.05;
static constexpr double hmin = 0.1 * hmax;
static constexpr double hausd = 0.5 * hmin;
static constexpr double ell = 0.4;
const constexpr Scalar dt = 4 * (hmax - hmin);
static constexpr double alpha = dt;

// Compliance
inline Scalar compliance(const GridFunction<FES>& w)
{
  auto& vh = w.getFiniteElementSpace();
  TrialFunction u(vh);
  TestFunction  v(vh);
  BilinearForm  bf(u, v);
  bf = LinearElasticityIntegral(u, v)(lambda, mu);
  return bf(w, w);
};

int main(int, char**)
{
  const char* meshFile = "../resources/examples/ShapeOptimization/LevelSetCantilever2D.mfem.mesh";

  // Load mesh
  MMG::Mesh th;
  th.load(meshFile);

  MMG::Optimizer().setHMax(hmax).setHMin(hmin).optimize(th);

  th.save("Omega0.mesh", IO::FileFormat::MEDIT);
  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  // Solver
  Solver::CG solver;

  // Optimization loop
 std::vector<double> obj;
  std::ofstream fObj("obj.txt");
  for (size_t i = 0; i < maxIt; i++)
  {
    th.getConnectivity().compute(1, 2);

    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    Alert::Info() << "   | Trimming mesh." << Alert::Raise;
    SubMesh trimmed = th.trim(Exterior);
    trimmed.save("Omega.mesh");

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
    elasticity.solve(solver);

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
    hilbert = Integral(alpha * alpha * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - FaceIntegral(Dot(Ae, e) - ell, Dot(n, w)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaN);
    hilbert.solve(solver);
    auto& dJ = g.getSolution();
    dJ.save("dJ.gf");
    vh.getMesh().save("dJ.mesh");

    // Update objective
    double objective = compliance(u.getSolution()) + ell * th.getVolume(Interior);
    obj.push_back(objective);
    fObj << objective << "\n";
    fObj.flush();
    Alert::Info() << "   | Objective: " << obj.back() << Alert::Raise;
    Alert::Info() << "   | Distancing domain." << Alert::Raise;

    P1 dh(th);
    auto dist = MMG::Distancer(dh).setInteriorDomain(Interior)
                                  .distance(th);

    // Advect the level set function
    Alert::Info() << "   | Advecting the distance function." << Alert::Raise;
    GridFunction norm(sh);
    norm = Frobenius(dJ);
    dJ /= norm.max();

    MMG::Advect(dist, dJ).step(dt);

    // Recover the implicit domain
    Alert::Info() << "   | Meshing the domain." << Alert::Raise;

    th = MMG::ImplicitDomainMesher().split(Interior, {Interior, Exterior})
                                    .split(Exterior, {Interior, Exterior})
                                    .setRMC(1e-6)
                                    .setHMax(hmax)
                                    .setHMin(hmin)
                                    .setHausdorff(hausd)
                                    .setAngleDetection(false)
                                    .setBoundaryReference(Gamma)
                                    .setBaseReferences(GammaD)
                                    .discretize(dist);

    MMG::Optimizer().setHMax(hmax)
                    .setHMin(hmin)
                    .setHausdorff(hausd)
                    .setAngleDetection(false)
                    .optimize(th);
  }

  Alert::Info() << "Saved final mesh to Omega.mesh" << Alert::Raise;

  return 0;
}

