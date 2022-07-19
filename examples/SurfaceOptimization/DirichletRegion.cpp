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
static constexpr int Gamma = 2;
static constexpr int GammaD = 3;
static constexpr int GammaN = 6;

static constexpr int SigmaD = 1;
static constexpr int SigmaN = 2;

static constexpr size_t maxIt = 250;

static constexpr double hmax = 0.05;
static constexpr double ell = 0.1;
static constexpr double alpha = 0.1;
static constexpr double epsilon = 0.01;
static constexpr double tgv = std::numeric_limits<double>::max();

/**
 * @brief Computes the conormal field.
 * @param[in] Sh Scalar finite element space
 * @oaram[in] Vh Vectorial finite element space
 */
GridFunction<H1> getShapeGradient(
    FiniteElementSpace<H1>& scalarFes,
    FiniteElementSpace<H1>& vecFes,
    GridFunction<H1>& dist,
    GridFunction<H1>& state,
    GridFunction<H1>& adjoint,
    const ScalarFunctionBase& g,
    Solver::Solver& solver);

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/dirichlet-region-example.mesh";

  // Load and build finite element spaces on the volumetric domain
  Mesh Omega;
  Omega.load(meshFile);

  auto J = [&](GridFunction<H1>& u)
  {
    return Integral(u).compute() - ell * Omega.getPerimeter(GammaN);
  };

  for (size_t i = 0; i < maxIt; i++)
  {
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    // Skin the mesh, computing the borders of the new regions
    Alert::Info() << "    | Skinning mesh." << Alert::Raise;
    auto dOmega = Omega.skin({{GammaD, SigmaD}, {GammaN, SigmaN}});

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

    auto solver = Solver::CG();

    auto h = [](double r)
    {
      if (r < -1.0)
        return 1.0;
      else if (r > 1.0)
        return 0.0;
      else
        return 1.0 - 1.0 / (1.0 + std::exp(4 * r / (r * r - 1.0)));
    };
    auto he = compose(h, ScalarFunction(dist) / epsilon) / epsilon;

    // State equation
    Alert::Info() << "    | Solving state equation." << Alert::Raise;
    ScalarFunction f = 1;
    ScalarFunction g = 2.0;

    TrialFunction u(Vh);
    TestFunction  v(Vh);
    Problem state(u, v);
    state = Integral(Grad(u), Grad(v))
          + BoundaryIntegral(he * u, v).over({Gamma, GammaD})
          - Integral(f, v)
          - BoundaryIntegral(g, v).over(GammaN);
    solver.solve(state);

    // Adjoint equation
    Alert::Info() << "    | Solving adjoint equation." << Alert::Raise;
    TrialFunction p(Vh);
    TestFunction  q(Vh);
    Problem adjoint(p, q);
    adjoint = Integral(Grad(p), Grad(q))
            + BoundaryIntegral(he * p, q).over({Gamma, GammaD})
            + Integral(ScalarFunction(u.getGridFunction()), q);
    solver.solve(adjoint);

    Alert::Info() << "    | Objective: " << J(u.getGridFunction())
                  << Alert::Raise;

    // Transfer the functions to the surfacic spaces
    GridFunction uS(VhS), pS(VhS);
    u.getGridFunction().transfer(uS);
    p.getGridFunction().transfer(pS);

    // Compute the shape gradient
    Alert::Info() << "    | Computing shape gradient." << Alert::Raise;
    auto gradS = getShapeGradient(VhS, ThS, distS, uS, pS, g, solver);

    // Transfer back the vector field to the whole space
    GridFunction grad(Th);
    gradS.transfer(grad);

    // Advect the distance function with the gradient
    Alert::Info() << "    | Advecting the distance function." << Alert::Raise;
    double gInf = std::max(grad.max(), -grad.min());
    double dt = 4 * hmax / gInf;
    MMG::Advect(dist, grad).avoidTimeTruncation().surface().step(dt);

    // Mesh only the surface part
    Alert::Info() << "    | Meshing the domain." << Alert::Raise;
    Omega = MMG::ImplicitDomainMesher().noSplit(Gamma).surface().discretize(dist);
  }

  return 0;
}

GridFunction<H1> getShapeGradient(
    FiniteElementSpace<H1>& scalarFes,
    FiniteElementSpace<H1>& vecFes,
    GridFunction<H1>& dist,
    GridFunction<H1>& state,
    GridFunction<H1>& adjoint,
    const ScalarFunctionBase& g,
    Solver::Solver& solver)
{
  auto n0 = VectorFunction{Dx(dist), Dy(dist), Dz(dist)};

  TrialFunction nx(scalarFes);
  TrialFunction ny(scalarFes);
  TrialFunction nz(scalarFes);
  TestFunction  v(scalarFes);

  // Conormal calculation
  Problem conormalX(nx, v);
  conormalX = Integral(alpha * Grad(nx), Grad(v))
          + Integral(nx, v)
          - Integral(n0.x(), v).over(GammaD);
  solver.solve(conormalX);

  Problem conormalY(ny, v);
  conormalY = Integral(alpha * Grad(ny), Grad(v))
          + Integral(ny, v)
          - Integral(n0.y(), v).over(GammaD);
  solver.solve(conormalY);

  Problem conormalZ(nz, v);
  conormalZ = Integral(alpha * Grad(nz), Grad(v))
          + Integral(nz, v)
          - Integral(n0.z(), v).over(GammaD);
  solver.solve(conormalZ);

  auto d = VectorFunction{
    nx.getGridFunction(), ny.getGridFunction(), nz.getGridFunction()};
  auto conormal = d / Pow(d.x() * d.x() + d.y() * d.y() + d.z() * d.z(), 0.5);
  auto expr = -ScalarFunction(state) * ScalarFunction(adjoint) - ell;

  TrialFunction gx(scalarFes);
  TrialFunction gy(scalarFes);
  TrialFunction gz(scalarFes);
  Problem velextX(gx, v);
  velextX = Integral(alpha * Grad(gx), Grad(v))
          + Integral(gx, v)
          + Integral(tgv * gx, v).over(GammaN)
          - BoundaryIntegral(expr * conormal.x(), v).over(SigmaD);
  solver.solve(velextX);

  Problem velextY(gy, v);
  velextY = Integral(alpha * Grad(gy), Grad(v))
          + Integral(gy, v)
          + Integral(tgv * gy, v).over(GammaN)
          - BoundaryIntegral(expr * conormal.y(), v).over(SigmaD);
  solver.solve(velextY);

  Problem velextZ(gz, v);
  velextZ = Integral(alpha * Grad(gz), Grad(v))
          + Integral(gz, v)
          + Integral(tgv * gz, v).over(GammaN)
          - BoundaryIntegral(expr * conormal.z(), v).over(SigmaD);
  solver.solve(velextZ);

  GridFunction grad(vecFes);
  grad = VectorFunction{gx.getGridFunction(), gy.getGridFunction(), gz.getGridFunction()};

  return grad;
}

