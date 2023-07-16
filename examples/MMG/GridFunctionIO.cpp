/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>

#include <Rodin/Geometry.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::External;
using namespace Rodin::Variational;

int main(int, char**)
{
  const size_t n = 16;
  MMG::Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Geometry::Triangle, n, n);

  P1 fes(mesh, 2);
  GridFunction gf(fes);
  gf = [](const Geometry::Point& p) { return Math::Vector{{p.x(), p.y()}}; };

  gf.save("MMG.medit.sol", IO::FileFormat::MEDIT);
  mesh.save("MMG.medit.mesh", IO::FileFormat::MEDIT);

  return 0;
}


