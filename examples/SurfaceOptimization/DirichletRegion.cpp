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
static constexpr int Gamma = 1;
static constexpr int GammaD = 5;
static constexpr int SigmaD = 1;
static constexpr int GammaN = 6;
static constexpr int SigmaN = 2;
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
  const char* meshFile = "Omega.mesh";

  // Load and build finite element spaces on the volumetric domain
  Mesh Omega;
  Omega.load(meshFile);

  FiniteElementSpace<H1> Vh(Omega);
  FiniteElementSpace<H1> Th(Omega, 3);

  // Skin the mesh and build finite element spaces on the submesh
  auto dOmega = Omega.skin({{GammaD, SigmaD}, {GammaN, SigmaN}});

  FiniteElementSpace<H1> VhS(dOmega);
  FiniteElementSpace<H1> ThS(dOmega, 3);

  // Distance the surface
  auto distSurf = MMG::Distancer(VhS).setInteriorDomain(GammaD)
                                     .distance(dOmega);

  // Transfer the distance to the whole domain
  GridFunction dist(Vh);
  distSurf.transfer(dist);

  dist.save("faulty.sol", IO::GridFunctionFormat::MEDIT);
  Omega.save("faulty.mesh", IO::MeshFormat::MEDIT);

  // Mesh only the surface part of the mesh
  Omega = MMG::ImplicitDomainMesher().surface()
                                     .discretize(dist);

  std::exit(1);

  ScalarFunction g = 2.0;

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

  auto J = [&](GridFunction<H1>& u)
  {
    return Integral(u).compute() - ell * Omega.getPerimeter(GammaN);
  };

  auto solver = Solver::UMFPack();

  // State equation
  ScalarFunction f = 1;
  TrialFunction u(Vh);
  TestFunction  v(Vh);
  Problem state(u, v);
  state = Integral(Grad(u), Grad(v))
        + BoundaryIntegral(he * u, v).over({Gamma, GammaD})
        - Integral(f, v)
        - BoundaryIntegral(g, v).over(GammaN);
  solver.solve(state);

  // Adjoint equation
  TrialFunction p(Vh);
  TestFunction  q(Vh);
  Problem adjoint(p, q);
  adjoint = Integral(Grad(p), Grad(q))
          + BoundaryIntegral(he * p, q).over({Gamma, GammaD})
          + Integral(ScalarFunction(u.getGridFunction()), q);
  solver.solve(adjoint);

  // Shape gradient
  GridFunction uS(VhS), pS(VhS);
  u.getGridFunction().transfer(uS);
  p.getGridFunction().transfer(pS);

  GridFunction<H1> grad = getShapeGradient(VhS, ThS, distSurf, uS, pS, g, solver);

  grad.save("grad.gf");

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

