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

/**
 * @brief Computes the conormal field.
 * @param[in] Sh Scalar finite element space
 * @oaram[in] Vh Vectorial finite element space
 */
GridFunction<H1> getShapeGradient(
    FiniteElementSpace<H1>& scalarFes,
    FiniteElementSpace<H1>& vecFes,
    GridFunction<H1>& dist,
    GridFunction<H1>& adjoint,
    const ScalarFunctionBase& g,
    Solver::Solver& solver, double ell, double alpha=0.1);

int main(int, char**)
{
  const char* meshFile = "Omega.mesh";

  int Gamma = 1;
  int GammaD = 5;
  int GammaN = 6;

  Mesh Omega;
  Omega.load(meshFile);

  ScalarFunction g = -1.0;
  double ell = 0.01;

  auto J = [&](GridFunction<H1>& u)
  {
    return Integral(u).compute() - ell * Omega.getPerimeter(GammaN);
  };

  FiniteElementSpace<H1> Vh(Omega);
  FiniteElementSpace<H1> Th(Omega, 3);

  auto dOmega = Omega.skin();

  // Compute extended conormal
  auto solver = Solver::UMFPack();

  // State equation
  ScalarFunction f = 1;
  TrialFunction u(Vh);
  TestFunction  v(Vh);
  Problem state(u, v);
  state = Integral(Grad(u), Grad(v))
        - Integral(f, v)
        - BoundaryIntegral(g, v).over(GammaN)
        + DirichletBC(u, ScalarFunction(0.0)).on(GammaD);
  solver.solve(state);

  // Adjoint equation
  TrialFunction p(Vh);
  TestFunction  q(Vh);
  Problem adjoint(p, q);
  adjoint = Integral(Grad(p), Grad(q))
          - Integral(q)
          + DirichletBC(p, ScalarFunction(0.0)).on(GammaD);
  solver.solve(adjoint);

  // Shape gradient
  FiniteElementSpace<H1> VhS(dOmega);
  FiniteElementSpace<H1> ThS(dOmega, 3);

  GridFunction pS(VhS);
  p.getGridFunction().transfer(pS);

  pS.save("pS.gf");

  GridFunction uS(VhS);
  u.getGridFunction().transfer(uS);

  dOmega.save("dOmega.mesh");
  uS.save("uS.gf");


  double alpha = 0.1;

  auto mmgMesh = Cast(dOmega).to<MMG::MeshS>();
  auto mmgDist = MMG::DistancerS().setInteriorDomain(GammaN).distance(mmgMesh).setMesh(mmgMesh);
  auto dist = Cast(mmgDist).to<GridFunction<H1>>(VhS);

  GridFunction<H1> grad = getShapeGradient(VhS, ThS, dist, pS, g, solver, ell, alpha);
  grad *= -1.0;

  dist.save("dist.gf");
  grad.save("grad.gf");

  return 0;
}

GridFunction<H1> getShapeGradient(
    FiniteElementSpace<H1>& scalarFes,
    FiniteElementSpace<H1>& vecFes,
    GridFunction<H1>& dist,
    GridFunction<H1>& adjoint,
    const ScalarFunctionBase& g,
    Solver::Solver& solver, double ell, double alpha)
{
  auto d0 = VectorFunction{Dx(dist), Dy(dist), Dz(dist)};
  auto n0 = d0 / (Pow(d0.x() * d0.x() + d0.y() * d0.y() + d0.z() * d0.z(), 0.5) + ScalarFunction(0.01));

  TrialFunction nx(scalarFes);
  TrialFunction ny(scalarFes);
  TrialFunction nz(scalarFes);
  TestFunction  v(scalarFes);

  Problem velextX(nx, v);
  velextX = Integral(alpha * Grad(nx), Grad(v))
          + Integral(nx, v)
          - Integral((ScalarFunction(adjoint) * g - ScalarFunction(ell)) * n0.x(), v);
  solver.solve(velextX);

  Problem velextY(ny, v);
  velextY = Integral(alpha * Grad(ny), Grad(v))
          + Integral(ny, v)
          - Integral((ScalarFunction(adjoint) * g - ScalarFunction(ell)) * n0.y(), v);
  solver.solve(velextY);

  Problem velextZ(nz, v);
  velextZ = Integral(alpha * Grad(nz), Grad(v))
          + Integral(nz, v)
          - Integral((ScalarFunction(adjoint) * g - ScalarFunction(ell)) * n0.z(), v);
  solver.solve(velextZ);

  GridFunction grad(vecFes);
  grad = VectorFunction{nx.getGridFunction(), ny.getGridFunction(), nz.getGridFunction()};

  return grad;
}

