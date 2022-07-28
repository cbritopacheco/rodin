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
using namespace Rodin::Variational;
using namespace Rodin::External;

// Parameters
static constexpr int Gamma = 6;
static constexpr int GammaD = 3;
static constexpr int GammaN = 2;

static constexpr int SigmaD = 1;
static constexpr int SigmaN = 2;

static constexpr size_t maxIt = 250;

static constexpr double hmax = 0.1;
static constexpr double alpha = 0.7;
static constexpr double epsilon = 0.1;
static constexpr double ell = 0.2;
static constexpr double tgv = std::numeric_limits<double>::max();

GridFunction<H1> getShapeGradient(
    FiniteElementSpace<H1>& vecFes, GridFunction<H1>& dist,
    const FunctionBase& expr, Solver::Solver& solver);

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/dirichlet-region-example.mesh";

  // Load and build finite element spaces on the volumetric domain
  Mesh Omega;
  Omega.load(meshFile);

  auto J = [&](GridFunction<H1>& u)
  {
    return Integral(u).compute() + ell * Omega.getPerimeter(GammaD);
  };

  std::ofstream fObj("obj.txt");
  for (size_t i = 0; i < maxIt; i++)
  {
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    // Skin the mesh, computing the borders of the new regions
    Alert::Info() << "    | Skinning mesh." << Alert::Raise;
    auto dOmega = Omega.skin();
    dOmega.trace({{{GammaD, Gamma}, SigmaD}, {{GammaN, Gamma}, SigmaN}});
    dOmega.save("trace.mesh", IO::FileFormat::MEDIT);

    // Build finite element spaces
    Alert::Info() << "    | Building finite element spaces." << Alert::Raise;
    FiniteElementSpace<H1> Vh(Omega);
    FiniteElementSpace<H1> Th(Omega, Omega.getSpaceDimension());

    FiniteElementSpace<H1> VhS(dOmega);
    FiniteElementSpace<H1> ThS(dOmega, dOmega.getSpaceDimension());

    Alert::Info() << "    | Distancing domain." << Alert::Raise;
    auto distS = MMG::Distancer(VhS).setInteriorDomain(GammaD)
                                    .distance(dOmega);

    GridFunction dist(Vh);
    distS.transfer(dist);

    Solver::CG solver;
    solver.setMaxIterations(1000);

    auto h = [](double r)
    {
      if (r < -1.0)
        return 1.0;
      else if (r > 1.0)
        return 0.0;
      else
        return 1.0 - 1.0 / (1.0 + std::exp(4 * r / (r * r - 1.0)));
    };
    auto he = compose(h, dist / epsilon) / epsilon;

    // State equation
    Alert::Info() << "    | Solving state equation." << Alert::Raise;
    ScalarFunction f = 1;
    ScalarFunction g = -2.0;

    TrialFunction u(Vh);
    TestFunction  v(Vh);
    Problem state(u, v);
    state = Integral(Grad(u), Grad(v))
          + BoundaryIntegral(he * u, v).over({Gamma, GammaD})
          - Integral(f, v)
          - BoundaryIntegral(g, v).over(GammaN);
    solver.solve(state);

    // Adjoint equation
    auto dj = -ScalarFunction(u.getGridFunction()) / Omega.getVolume();
    Alert::Info() << "    | Solving adjoint equation." << Alert::Raise;
    TrialFunction p(Vh);
    TestFunction  q(Vh);
    Problem adjoint(p, q);
    adjoint = Integral(Grad(p), Grad(q))
            + BoundaryIntegral(he * p, q).over({Gamma, GammaD})
            - Integral(dj, q);
    solver.solve(adjoint);

    double objective = J(u.getGridFunction());
    Alert::Info() << "    | Objective: " << objective
                  << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

    // Transfer the functions to the surfacic spaces
    GridFunction uS(VhS), pS(VhS);
    u.getGridFunction().transfer(uS);
    p.getGridFunction().transfer(pS);

    // Compute the shape gradient
    Alert::Info() << "    | Computing shape gradient." << Alert::Raise;
    auto expr = 1. / (epsilon * epsilon) * uS * pS + ell;
    auto gradS = getShapeGradient(ThS, distS, expr, solver);

    // Transfer back the vector field to the whole space
    GridFunction grad(Th);
    gradS.transfer(grad);

    grad *= -1.0;

    // Advect the distance function with the gradient
    Alert::Info() << "    | Advecting the distance function." << Alert::Raise;
    double gInf = std::max(gradS.max(), -gradS.min());
    double dt = hmax / gInf;
    MMG::Advect(dist, grad).avoidTimeTruncation().surface().step(dt);

    // Mesh only the surface part
    Alert::Info() << "    | Meshing the domain." << Alert::Raise;
    Omega = MMG::ImplicitDomainMesher().noSplit(GammaN)
                                       .split(GammaD, {GammaD, Gamma})
                                       .split(Gamma, {GammaD, Gamma})
                                       .setHMax(hmax)
                                       .surface()
                                       .discretize(dist);
    MMG::MeshOptimizer().setHMax(hmax).optimize(Omega);

    Omega.save("Omega.mesh", IO::FileFormat::MEDIT);

    // Omega.skin().save("out/dOmega." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
  }

  return 0;
}

GridFunction<H1> getShapeGradient(
    FiniteElementSpace<H1>& vecFes, GridFunction<H1>& dist,
    const FunctionBase& expr, Solver::Solver& solver)
{
  TrialFunction d(vecFes);
  TestFunction  v(vecFes);

  Problem conormalExt(d, v);
  conormalExt = Integral(alpha * Jacobian(d), Jacobian(v))
              + Integral(d, v)
              - BoundaryIntegral(Grad(dist).traceOf(GammaD), v).over(SigmaD);
  solver.solve(conormalExt);

  const auto& cnd = d.getGridFunction();
  const auto cn = cnd / Pow(cnd.x() * cnd.x() + cnd.y() * cnd.y() + cnd.z() * cnd.z(), 0.5);

  TrialFunction grad(vecFes);
  Problem velExt(grad, v);
  velExt = Integral(alpha * Jacobian(grad), Jacobian(v))
         + Integral(grad, v)
         + Integral(tgv * grad, v).over(GammaN)
         - BoundaryIntegral(expr * cn, v).over(SigmaD);
  solver.solve(velExt);

  return grad.getGridFunction();
}

