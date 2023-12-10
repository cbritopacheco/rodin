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
static constexpr int Interior = 3, Exterior = 2;

// Define boundary attributes
static constexpr int GammaD = 999, GammaN = 2, Gamma = 10, Gamma0 = 666;

// Lamé coefficients
static constexpr double mu = 0.3846;
static constexpr double lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 300;
static constexpr double eps = 1e-12;
static constexpr double hgrad = 1.6;
static constexpr double ell = 0.1;
static double elementStep = 0.5;
static double hmax = 0.2;
static double hmin = 0.1 * hmax;
static double hausd = 0.5 * hmin;
static size_t hmaxIt = maxIt / 2;
static double alpha = 4 * hmax;

using FES = VectorP1<Context::Sequential>;

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
  const char* meshFile = "step.0.mesh";

  // Load mesh
  MMG::Mesh th;
  th.load(meshFile, IO::FileFormat::MEDIT);

  MMG::Optimizer().setHMax(hmax)
                  .setHMin(hmin)
                  .setHausdorff(hausd)
                  .setAngleDetection(false)
                  .optimize(th);

  th.save("Omega0.mesh");
  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  // Solver
  // Eigen::ConjugateGradient<Math::SparseMatrix> ecg;
  // Solver::EigenSolver<
  //   Eigen::ConjugateGradient<Math::SparseMatrix>, Math::SparseMatrix, Math::Vector> solver(ecg);
  Solver::CG solver;
  solver.setMaxIterations(200);

  // Optimization loop
  std::vector<double> obj;
  std::ofstream fObj("obj.txt");
  for (size_t i = 0; i < maxIt; i++)
  {
    th.getConnectivity().compute(th.getDimension() - 1, th.getDimension());
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    Alert::Info() << "   | Trimming mesh." << Alert::Raise;
    SubMesh trimmed = th.trim(Exterior);
    trimmed.save("trimmed.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "   | Building finite element spaces." << Alert::Raise;
    const size_t d = th.getSpaceDimension();
    P1 sh(th);
    P1 vh(th, d);

    Alert::Info() << "   | Distancing domain." << Alert::Raise;
    auto dist = MMG::Distancer(sh).setInteriorDomain(Interior)
                                  .distance(th);

    P1 shInt(trimmed);
    P1 vhInt(trimmed, d);

    Alert::Info() << "   | Solving state equation." << Alert::Raise;
    auto f = VectorFunction{0, 0, -1};
    TrialFunction u(vhInt);
    TestFunction  v(vhInt);

    // Elasticity equation
    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - FaceIntegral(f, v).over(GammaN)
               + DirichletBC(u, VectorFunction{0, 0, 0}).on(GammaD);
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
    hilbert = Integral(alpha * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - FaceIntegral(Dot(Ae, e) - ell, Dot(n, w)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaN);
    hilbert.solve(solver);

    auto& dJ = g.getSolution();

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
    th.save("dJ.mesh");
    dJ.save("dJ.gf");
    const Scalar k = hmin;
    const Scalar dt = k;
    MMG::Advect(dist, dJ).step(dt);

    // Recover the implicit domain
    Alert::Info() << "   | Meshing the domain." << Alert::Raise;

    th = MMG::ImplicitDomainMesher().split(Interior, {Interior, Exterior})
                                    .split(Exterior, {Interior, Exterior})
                                    .setRMC(1e-5)
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

    th.save("Omega.mesh", IO::FileFormat::MEDIT);
  }

  Alert::Info() << "Saved final mesh to Omega.mesh" << Alert::Raise;

  return 0;
}

