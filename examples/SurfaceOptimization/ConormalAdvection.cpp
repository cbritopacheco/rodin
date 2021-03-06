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
  const char* meshFile = "Omega.mesh";
  const double pi = std::atan(1) * 4;

  int Interior = 2, Exterior = 3;
  int Gamma = 4;

  Mesh Omega;
  Omega.load(meshFile);


  // GridFunction f(Vh);
  // f = VectorFunction{
  //   [](const double* x, int){ return 3 * x[0] * x[0] + sin(10 * x[1]); },
  //   [](const double* x, int){ return 3 * x[0] + cos(10 * x[1]); }};

  FiniteElementSpace<H1> Vh(Omega);

  // Sphere radius
  double r = 1;

  // Hole radius
  double rr = 0.2;

  // Hole centers on sphere
  std::vector<std::array<double, 3>> cs = {
    // { r * sin(0) * cos(0), r * sin(0) * sin(0), r * cos(0) },
    { r * sin(pi / 2) * cos(pi), r * sin(pi / 2) * sin(pi), r * cos(pi / 2) },
    { r * sin(pi / 2) * cos(pi / 2), r * sin(pi / 2) * sin(pi / 2), r * cos(pi / 2) }
  };

  // Geodesic distance
  auto gd = [&](const double* x, const double* c, int)
            {
              return std::acos((x[0] * c[0] + x[1] * c[1] + x[2] * c[2]));
            };

  // Function for generating holes
  auto f = ScalarFunction(
      [&](const double* x, int dim) -> double
      {
        double d = std::numeric_limits<double>::max();
        for (const auto& c : cs)
        {
          double dd = gd(x, c.data(), dim) - rr;
          d = std::min(d, dd);
        }
        if (d <= rr)
          return d;
        else
          return -d;
      });

  for (const auto& c : cs)
  {
    for (int i = 0; i < Omega.getHandle().GetNBE(); i++)
    {
      mfem::Array<int> vs;
      Omega.getHandle().GetBdrElement(i)->GetVertices(vs);
      double d = std::numeric_limits<double>::max();
      for (int v = 0; v < vs.Size(); v++)
      {
        double dd = gd(Omega.getHandle().GetVertex(vs[v]), c.data(), 3) - rr;
        d = std::min(d, dd);
        if (d <= rr)
          Omega.getHandle().SetBdrAttribute(i, 6);
        // else
        // Omega.getHandle().SetBdrAttribute(i, 1);
      }
    }
  }

  GridFunction dist(Vh);
  dist = f;

  Omega.save("Omega.mesh");
  dist.save("dist.gf");

  // Cast(Omega).to<MMG::Mesh3D>().save("mmg.mesh");
  // Cast(dist).to<MMG::ScalarSolution3D>(mmgMesh).save("mmg.sol");

  // mmgSol.save("mmg.sol");

  // MMG::Mesh3D mmgMesh;
  // mmgMesh.load("test.mesh");

  // MMG::ScalarSolution3D mmgDist(mmgMesh);
  // mmgDist.load("test.sol");

  // auto Omega = Cast(mmgMesh).to<SerialMesh>();

  // FiniteElementSpace<H1> Vh(Omega);
  // FiniteElementSpace<H1> Th(Omega, 3);

  // auto dist = Cast(mmgDist).to<GridFunction<H1>>(Vh);


  // auto skin = Omega.skin();
  // FiniteElementSpace<H1> VhS(skin);
  // FiniteElementSpace<H1> ThS(skin, 3);
  // GridFunction distS(VhS);
  // dist.transfer(distS);


  return 0;
}
