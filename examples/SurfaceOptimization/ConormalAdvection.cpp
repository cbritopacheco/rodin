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

int main(int, char**)
{
  const char* meshFile = "miaow.mesh";
  const double pi = std::atan(1) * 4;

  Mesh Omega;
  Omega.load(meshFile);

  H1 Vh(Omega);

  // Sphere radius
  double r = 1;

  // Hole radius
  double rr = 0.2;

  // Hole centers on sphere
  std::vector<std::array<double, 3>> cs = {
    { r * sin(0) * cos(0), r * sin(0) * sin(0), r * cos(0) }
    // { r * sin(pi / 2) * cos(pi), r * sin(pi / 2) * sin(pi), r * cos(pi / 2) },
    // { r * sin(pi / 2) * cos(pi / 2), r * sin(pi / 2) * sin(pi / 2), r * cos(pi / 2) }
  };

  // Geodesic distance
  auto gd = [&](const Point& c, std::array<double, 3> x)
            {
              return std::acos((x[0] * c(0) + x[1] * c(1) + x[2] * c(2)));
            };

  // Function for generating holes
  auto f = ScalarFunction(
      [&](const Point& v) -> double
      {
        double d = std::numeric_limits<double>::max();
        for (const auto& c : cs)
        {
          double dd = gd(v, c) - rr;
          d = std::min(d, dd);
        }
        return d;
      });

  GridFunction dist(Vh);
  dist.projectOnBoundary(f);

  Omega = MMG::ImplicitDomainMesher()
                             .split(6, {3, 6})
                             .noSplit(2)
                             .setHMax(0.05)
                             .surface()
                             .discretize(dist);

  Omega.save("leeel.mesh", IO::FileFormat::MEDIT);

  return 0;
}
