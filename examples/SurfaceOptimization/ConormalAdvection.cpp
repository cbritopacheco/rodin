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
  const char* meshFile = "miaow.o.mesh";
  const double pi = std::atan(1) * 4;


  Mesh Omega;
  Omega.load(meshFile, IO::FileFormat::MEDIT);

  // Omega.save("Omega.mesh", IO::MeshFormat::MEDIT);

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
        return d;
      });

  GridFunction dist(Vh);
  dist.projectOnBoundary(f);
  Omega.save("dist.mesh");
  dist.save("dist.gf");

  Omega = MMG::ImplicitDomainMesher().split(2, {2, 6})
                             .noSplit(3)
                             .setHMax(0.05)
                             .surface()
                             .discretize(dist);

  Omega.save("leeel.mesh", IO::FileFormat::MEDIT);

  return 0;
}
