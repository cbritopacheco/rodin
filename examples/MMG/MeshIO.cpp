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
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, n, n);
  mesh.getConnectivity().compute(1, 2);

  mesh.setCorner(0);
  mesh.setCorner(15);
  mesh.setCorner(240);
  mesh.setCorner(255);

  for (auto it = mesh.getBoundary(); !it.end(); ++it)
    mesh.setRidge(it->getIndex());

  mesh.save("MMG.medit.mesh", IO::FileFormat::MEDIT);

  MMG::Mesh mesh2;
  mesh2.load("MMG.medit.mesh", IO::FileFormat::MEDIT);
  mesh2.save("Woof.medit.mesh", IO::FileFormat::MEDIT);

  return 0;
}

