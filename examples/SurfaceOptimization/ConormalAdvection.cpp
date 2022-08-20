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


int main(int, char**)
{
  // const char* meshFile = "miaow.o.mesh";
  // const double pi = std::atan(1) * 4;


  // Mesh Omega;
  // Omega.load(meshFile, IO::FileFormat::MEDIT);

  // // Omega.save("Omega.mesh", IO::MeshFormat::MEDIT);

  // H1 Vh(Omega);

  // // Sphere radius
  // double r = 1;

  // // Hole radius
  // double rr = 0.2;

  // // Hole centers on sphere
  // std::vector<std::array<double, 3>> cs = {
  //   // { r * sin(0) * cos(0), r * sin(0) * sin(0), r * cos(0) },
  //   { r * sin(pi / 2) * cos(pi), r * sin(pi / 2) * sin(pi), r * cos(pi / 2) },
  //   { r * sin(pi / 2) * cos(pi / 2), r * sin(pi / 2) * sin(pi / 2), r * cos(pi / 2) }
  // };

  // // Geodesic distance
  // auto gd = [&](const double* x, const double* c, int)
  //           {
  //             return std::acos((x[0] * c[0] + x[1] * c[1] + x[2] * c[2]));
  //           };

  // // Function for generating holes
  // auto f = ScalarFunction(
  //     [&](const double* x, int dim) -> double
  //     {
  //       double d = std::numeric_limits<double>::max();
  //       for (const auto& c : cs)
  //       {
  //         double dd = gd(x, c.data(), dim) - rr;
  //         d = std::min(d, dd);
  //       }
  //       return d;
  //     });

  // GridFunction dist(Vh);
  // dist.projectOnBoundary(f);
  // Omega.save("dist.mesh");
  // dist.save("dist.gf");

  // Omega = MMG::ImplicitDomainMesher().split(2, {2, 6})
  //                            .noSplit(3)
  //                            .setHMax(0.05)
  //                            .surface()
  //                            .discretize(dist);

  // Omega.save("leeel.mesh", IO::FileFormat::MEDIT);

  return 0;
}
