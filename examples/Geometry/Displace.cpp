/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  size_t n = 32;
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, n, n);
  mesh.scale(1. / (n - 1));

  auto displacement =
    VectorFunction{0, [&](const Geometry::Point& p) { return (0.5 - p.x()) * (0.5 - p.x()); }};

  mesh.save("Original.mesh");
  mesh.displace(displacement);
  mesh.save("Displaced.mesh");
}


