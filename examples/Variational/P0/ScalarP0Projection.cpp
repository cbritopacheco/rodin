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
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, 32, 32);

  P0 fes(mesh);
  GridFunction gf(fes);

  gf = [](const Point& p) { return p.x() * p.x() +  p.y() * p.y(); };

  mesh.save("Projection.mesh");
  gf.save("Projection.gf");
}



