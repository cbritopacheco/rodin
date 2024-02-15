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
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });

  P1 fes(mesh, 2);
  GridFunction gf(fes);
  gf = [](const Geometry::Point& p) { return Math::Vector{{p.x(), p.y()}}; };

  mesh.save("function.mesh");
  gf.save("function.gf");

  return 0;
}


