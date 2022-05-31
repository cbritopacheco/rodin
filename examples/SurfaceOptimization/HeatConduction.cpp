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
auto getConormal(
    FiniteElementSpace<H1>& scalarFes,
    FiniteElementSpace<H1>& vecFes,
    GridFunction<H1>& dist,
    Solver::Solver& solver, double alpha=0.1)
{
  auto n0 = VectorFunction{Dx(dist), Dy(dist), Dz(dist)};

  TrialFunction nx(scalarFes);
  TrialFunction ny(scalarFes);
  TrialFunction nz(scalarFes);
  TestFunction  v(scalarFes);

  Problem velextX(nx, v);
  velextX = Integral(alpha * Grad(nx), Grad(v))
          + Integral(nx, v)
          - Integral(n0.x(), v);
  solver.solve(velextX);

  Problem velextY(ny, v);
  velextY = Integral(alpha * Grad(ny), Grad(v))
          + Integral(ny, v)
          - Integral(n0.y(), v);
  solver.solve(velextY);

  Problem velextZ(nz, v);
  velextZ = Integral(alpha * Grad(nz), Grad(v))
          + Integral(nz, v)
          - Integral(n0.z(), v);
  solver.solve(velextZ);

  auto n = VectorFunction{nx.getGridFunction(), ny.getGridFunction(), nz.getGridFunction()};
  auto norm = Pow(n.x() * n.x() + n.y() * n.y() + n.z() * n.z(), 0.5);

  GridFunction conormal(vecFes);
  conormal = n / norm;

  return conormal;
}

int main(int, char**)
{
  const char* meshFile = "rodin.mesh";

  int Interior = 2, Exterior = 3;
  int Gamma = 4;
  int GammaD = 5;
  int GammaN = 6;

  Mesh Omega;
  Omega.load(meshFile);

  FiniteElementSpace<H1> Vh(Omega);
  FiniteElementSpace<H1> Th(Omega, 3);

  auto mmgMesh = Cast(Omega).to<MMG::MeshS>();
  auto mmgDist = MMG::DistancerS().distance(mmgMesh).setMesh(mmgMesh);
  auto dist = Cast(mmgDist).to<GridFunction<H1>>(Vh);

  // Compute extended conormal
  auto solver = Solver::UMFPack();

  auto conormal = getConormal(Vh, Th, dist, solver);

  conormal.save("conormal.gf");

  // double epsilon = 0.01;
  // ScalarFunction h = 1. / epsilon;
  // ScalarFunction f = 1;
  // ScalarFunction g = 1;
  // TrialFunction u(Vh);
  // TestFunction  v(Vh);
  // Problem state(u, v);
  // state = Integral(Grad(u), Grad(v))
  //       + BoundaryIntegral(h * u, v)
  //       - Integral(f, v).over({Gamma, GammaD})
  //       - Integral(g, v).over(GammaN);

  // ScalarFunction jp = 1;
  // TrialFunction p(Vh);
  // TestFunction  q(Vh);
  // Problem adjoint(p, q);
  // adjoint = Integral(Grad(p), Grad(q))
  //         + BoundaryIntegral(h * p, q)
  //         - Integral(jp, q).over({Gamma, GammaD});


  return 0;
}

