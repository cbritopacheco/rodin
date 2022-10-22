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


// Define interior and exterior for level set discretization
static constexpr int Interior = 3, Exterior = 2;

// Define boundary attributes
static constexpr int GammaD = 1, GammaN = 2, Gamma = 10, Gamma0 = 3;

// Lamé coefficients
static constexpr double mu = 0.3846;
static constexpr double lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 300;
static constexpr double eps = 1e-12;
static constexpr double hgrad = 1.6;
static constexpr double ell = 0.2;
static double elementStep = 0.5;
static double hausd = 0.01;
static double hmax = 0.1;
static double hmin = 0.1 * hmax;
static size_t hmaxIt = maxIt / 2;

// Compliance
double compliance(GridFunction<H1<Context::Serial>>& w);

int main(int, char**)
{
  const char* meshFile = "step.0.o.mesh";
  // const char* meshFile = "Omega.mesh";

  // Load mesh
  Mesh Omega;
  Omega.load(meshFile, IO::FileFormat::MEDIT);

  Omega.save("Omega0.mesh");
  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  // CG
  Solver::CG solver;
  solver.setMaxIterations(1000).setRelativeTolerance(1e-6).setAbsoluteTolerance(0.0);

  // Optimization loop
  std::vector<double> obj;
  std::ofstream fObj("obj.txt", std::ios_base::app);
  size_t i = 0;
  while (i < maxIt)
  {
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    hmax = 0.1 - (0.1 - 0.04) * static_cast<double>(std::min(i, hmaxIt)) / hmaxIt;
    hmin = std::max(0.01 * hmax, hmin);
    hausd = std::max(0.0001, hausd);

    Alert::Info() << "    | Parameters:\n"
                  << "          | ell: " << ell << '\n'
                  << "          | hmax: " << hmax << '\n'
                  << "          | hmin: " << hmin << '\n'
                  << "          | hausd: " << hausd << '\n'
                  << "          | hgrad: " << hgrad << '\n'
                  << "          | elementStep: " << elementStep
                  << Alert::Raise;

    try
    {
      Alert::Info() << "    | Optimizing the domain." << Alert::Raise;
      // MMG::MeshOptimizer().setHMax(hmax)
      //                     .setHMin(hmin)
      //                     .setGradation(hgrad)
      //                     .setHausdorff(hausd)
      //                     .optimize(Omega);
      Omega.save("optimized.mesh", IO::FileFormat::MEDIT);
    }
    catch (Rodin::Alert::Exception& err)
    {
      Alert::Warning() << " | Failed, skipping iteration." << Alert::Raise;
      // hmin /= 1.6;
      // hausd /= 1.2;
      // elementStep /= 1.6;
      // continue;
    }

    Alert::Info() << "    | Trimming mesh." << Alert::Raise;
    SubMesh trimmed = Omega.keep(Interior);

    Alert::Info() << "    | Building finite element spaces." << Alert::Raise;
    int d = Omega.getSpaceDimension();
    H1 Vh(Omega, d);
    H1 VhInt(trimmed, d);

    Alert::Info() << "    | Solving state equation." << Alert::Raise;
    // auto f = VectorFunction{0, 0, -1};

    auto f = 10 * VectorFunction(
        0,
        [](const Vertex& v) -> double
        {
           return v.z() - 0.5;
        },
        [](const Vertex& v) -> double
        {
          return -v.y() + 0.5;
        }
        );

    GridFunction pull(VhInt);
    pull.projectOnBoundary(f, GammaN);
    trimmed.save("pull.mesh");
    pull.save("pull.gf");

    TrialFunction uInt(VhInt);
    TestFunction  vInt(VhInt);

    // Elasticity equation
    Problem elasticity(uInt, vInt);
    elasticity = Integral(lambda * Div(uInt), Div(vInt))
               + Integral(
                   mu * (Jacobian(uInt) + Jacobian(uInt).T()), 0.5 * (Jacobian(vInt) + Jacobian(vInt).T()))
               - BoundaryIntegral(f, vInt).over(GammaN)
               + DirichletBC(uInt, VectorFunction{0, 0, 0}).on(GammaD);
    solver.solve(elasticity);

    // Transfer solution back to original domain
    GridFunction u(Vh);
    uInt.getGridFunction().transfer(u);

    uInt.getGridFunction().save("u.gf");
    trimmed.save("u.mesh");

    Alert::Info() << "    | Computing shape gradient." << Alert::Raise;
    auto e = 0.5 * (Jacobian(u).traceOf(Interior) + Jacobian(u).traceOf(Interior).T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    auto n = -Normal(d); // Bug: For some reason orientation does not match up
                         // with that of the 2D case.

    // Hilbert extension-regularization procedure
    double alpha = elementStep * hmax;
    TrialFunction g(Vh);
    TestFunction  v(Vh);
    Problem hilbert(g, v);
    hilbert = Integral(alpha * Jacobian(g), Jacobian(v))
            + Integral(g, v)
            - BoundaryIntegral(Dot(Ae, e) - ell, Dot(n, v)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaN)
            ;
    solver.solve(hilbert);

    g.getGridFunction().save("g.gf");
    Omega.save("g.mesh");

    // Update objective
    double objective = compliance(u) + ell * Omega.getVolume(Interior);
    obj.push_back(objective);
    Alert::Info() << "    | Objective: " << obj.back() << Alert::Raise;

    Alert::Info() << "    | Distancing domain." << Alert::Raise;
    H1 Dh(Omega);
    auto dist = MMG::Distancer(Dh).setInteriorDomain(Interior)
                                  .distance(Omega);

    // Advect the level set function
    GridFunction gNorm(Dh);
    gNorm = ScalarFunction(
        [&](const Vertex& v) -> double
        {
          mfem::Vector val = g.getGridFunction()(v);
          return val.Norml2();
        });
    double linfty = gNorm.max();
    double dt = elementStep * hmin / gNorm.max();
    Alert::Info() << "    | Advecting the distance function. || g || = " << linfty
      << Alert::Raise;

    try
    {
      MMG::Advect(dist, g.getGridFunction()).step(dt);
    }
    catch (Rodin::Alert::Exception& e)
    {
      Alert::Warning() << " | Failed, skipping iteration." << Alert::Raise;
      continue;
    }

    // Recover the implicit domain
    try
    {
      Alert::Info() << "    | Meshing the domain." << Alert::Raise;
      Omega = MMG::ImplicitDomainMesher().setRMC(1e-5)
                                         .setHMax(hmax)
                                         .setHMin(hmin)
                                         .setHausdorff(hausd)
                                         .setGradation(hgrad)
                                         .setBaseReferences(GammaD)
                                         .setBoundaryReference(Gamma)
                                         .discretize(dist)
                                          // .split(Interior, {Interior, Exterior})
                                         // .split(Exterior, {Interior, Exterior})
                                         ;
      elementStep = 1.1 * elementStep > 2.0 ? 2.0 : elementStep * 1.1;
      hausd = 1.1 * hausd > 0.01 ? 0.01 : 1.1 * hausd;
      hmin = 1.1 * hmin > 0.2 * hmax ? 0.2 * hmax : 1.1 * hmin;
    }
    catch (Rodin::Alert::Exception& err)
    {
      Alert::Warning() << " | Failed, skipping iteration." << Alert::Raise;
      hmin /= 1.6;
      hausd /= 1.4;
      elementStep /= 1.6;
      continue;
    }

    Alert::Info() << "    | Saving results." << Alert::Raise;
    Omega.save("Omega.mesh", IO::FileFormat::MEDIT);
    Omega.save("out/Omega." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
    auto shape = Omega.keep(Interior);
    shape.save("out/trimmed." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
    shape.save("trimmed.mesh");

    // Test for convergence
    // if (obj.size() >= 2 && abs(obj[i] - obj[i - 1]) < eps)
    // {
    //   Alert::Info() << "Convergence!" << Alert::Raise;
    //   break;
    // }

    fObj << objective << "\n";
    fObj.flush();
    i++;
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
  bf = Integral(lambda * Div(u), Div(v)).over(Interior)
     + Integral(
         mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T())).over(Interior);
  return bf(w, w);
};


