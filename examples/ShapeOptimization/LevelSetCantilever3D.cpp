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
static double hmax = 0.2;
static double hmin = 0.1 * hmax;
static size_t hmaxIt = maxIt / 2;

// Compliance
double compliance(GridFunction<H1<Context::Serial>>& w);

int main(int, char**)
{
  const char* meshFile = "step.0.mesh";
  // const char* meshFile = "Omega.mesh";

  // Load mesh
  MMG::Mesh Omega;
  Omega.load(meshFile, IO::FileFormat::MEDIT);
  Omega.save("miaow.mesh", IO::FileFormat::MEDIT);

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

   hmax = 0.2 - (0.2 - 0.05) * static_cast<double>(std::min(i, hmaxIt)) / hmaxIt;
   hmin = std::max(0.1 * hmax, hmin);
   hausd = std::max(0.0001, hausd);

   Alert::Info() << "   | Parameters:\n"
            << "       | ell: " << ell << '\n'
            << "       | hmax: " << hmax << '\n'
            << "       | hmin: " << hmin << '\n'
            << "       | hausd: " << hausd << '\n'
            << "       | hgrad: " << hgrad << '\n'
            << "       | elementStep: " << elementStep
            << Alert::Raise;

   try
   {
    if (i == 0)
    {
      Alert::Info() << "   | Optimizing the domain." << Alert::Raise;
      MMG::MeshOptimizer().setHMax(hmax)
                   .setHMin(hmin)
                   .setGradation(hgrad)
                   .setHausdorff(hausd)
                   .setAngleDetection(i == 0)
                   .optimize(Omega);
      Omega.save("optimized.mesh", IO::FileFormat::MEDIT);
    }
   }
   catch (Rodin::Alert::Exception& err)
   {
    Alert::Warning() << " | Failed, skipping iteration." << Alert::Raise;
   }

   Alert::Info() << "   | Trimming mesh." << Alert::Raise;
   SubMesh trimmed = Omega.keep(Interior);

   Alert::Info() << "   | Building finite element spaces." << Alert::Raise;
   int d = Omega.getSpaceDimension();
   H1 Vh(Omega, d);
   H1 VhInt(trimmed, d);

   Alert::Info() << "   | Solving state equation." << Alert::Raise;
   auto f = VectorFunction{0, 0, -1};

   TrialFunction uInt(VhInt);
   TestFunction  vInt(VhInt);

   // Elasticity equation
   Problem elasticity(uInt, vInt);
   elasticity = Integral(lambda * Div(uInt), Div(vInt))
          + Integral(
             mu * (Jacobian(uInt) + Jacobian(uInt).T()), 0.5 * (Jacobian(vInt) + Jacobian(vInt).T()))
          - BoundaryIntegral(f, vInt).over(GammaN)
          + DirichletBC(uInt, VectorFunction{0, 0, 0}).on(GammaD);
   elasticity.solve(solver);

   Alert::Info() << "   | Computing shape gradient." << Alert::Raise;
   auto jac = Jacobian(uInt.getSolution());
   jac.traceOf(Interior);

   auto e = 0.5 * (jac + jac.T());
   auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
   auto n = Normal(d);
   n.traceOf(Interior);

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
   hilbert.solve(solver);

   g.getSolution().save("g.gf");
   Omega.save("g.mesh");

   // Update objective
   double objective = compliance(uInt.getSolution()) + ell * Omega.getVolume(Interior);
   obj.push_back(objective);
   Alert::Info() << "   | Objective: " << obj.back() << Alert::Raise;

   Alert::Info() << "   | Distancing domain." << Alert::Raise;
   H1 Dh(Omega);
   auto dist = MMG::Distancer(Dh).setInteriorDomain(Interior)
                       .distance(Omega);

   // Compute norm of gradient
   GridFunction gNorm(Dh);
   gNorm = ScalarFunction(
      [&](const Point& v) -> double
      {
        Math::Vector val = g.getSolution()(v);
        return val.norm();
      });
   double linfty = gNorm.max();
   double dt = elementStep * hmax / gNorm.max();

   try
   {
    Alert::Info()
      << "   | Advecting the distance function. || g || = " << linfty
      << Alert::Raise;
    MMG::Advect(dist, g.getSolution()).step(dt);
   }
   catch (Rodin::Alert::Exception& e)
   {
    Alert::Warning() << " | Failed, skipping iteration." << Alert::Raise;
    continue;
   }

   // Recover the implicit domain
   try
   {
    Alert::Info() << "   | Meshing the domain." << Alert::Raise;
    Omega = MMG::ImplicitDomainMesher().setRMC(1e-5)
                            .setHMax(hmax)
                            .setHMin(hmin)
                            .setHausdorff(hausd)
                            .setGradation(hgrad)
                            .setAngleDetection(false)
                            .setBaseReferences(GammaD)
                            .setBoundaryReference(Gamma)
                            .discretize(dist)
                            // .split(Interior, {Interior, Exterior})
                            // .split(Exterior, {Interior, Exterior})
                            ;
    elementStep = 1.1 * elementStep > 1.0 ? 1.0 : elementStep * 1.1;
    hausd = 1.1 * hausd > 0.01 ? 0.01 : 1.1 * hausd;
    hmin = 1.1 * hmin > 0.5 * hmax ? 0.5 * hmax : 1.1 * hmin;
   }
   catch (Rodin::Alert::Exception& err)
   {
    Alert::Warning() << " | Failed, skipping iteration." << Alert::Raise;
    hmin /= 1.6;
    hausd /= 1.4;
    elementStep /= 1.6;
    continue;
   }

   Alert::Info() << "   | Saving results." << Alert::Raise;
   Omega.save("Omega.mesh", IO::FileFormat::MEDIT);
   Omega.save("out/Omega." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
   auto shape = Omega.keep(Interior);
   shape.save("out/trimmed." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
   shape.save("trimmed.mesh");

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


