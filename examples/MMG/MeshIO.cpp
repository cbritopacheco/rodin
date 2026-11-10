/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [mmg-io-source]
#include <cstdlib>
#include <iostream>
#include <Rodin/Geometry.h>
#include <Rodin/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main()
{
  //! [mmg-io-create]
  constexpr size_t nodes = 8;
  MMG::Mesh mesh;
  mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  mesh.scale(1.0 / (nodes - 1));
  mesh.getConnectivity().compute(1, 2);
  mesh.setCorner(0);
  mesh.setCorner(nodes - 1);
  mesh.setCorner(nodes * (nodes - 1));
  mesh.setCorner(nodes * nodes - 1);
  mesh.setRequiredVertex(0);
  //! [mmg-io-create]

  //! [mmg-io-roundtrip]
  mesh.save("MMG.medit.mesh", IO::FileFormat::MEDIT);
  MMG::Mesh loaded;
  loaded.load("MMG.medit.mesh", IO::FileFormat::MEDIT);
  loaded.save("MMG.copy.medit.mesh", IO::FileFormat::MEDIT);
  //! [mmg-io-roundtrip]

  const bool valid = loaded.getVertexCount() == mesh.getVertexCount() &&
    loaded.getCellCount() == mesh.getCellCount() &&
    loaded.getCorners() == mesh.getCorners() &&
    loaded.getRequiredVertices() == mesh.getRequiredVertices();
  std::cout << "Vertices: " << loaded.getVertexCount()
            << ", triangles: " << loaded.getCellCount() << std::endl;
  return valid ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [mmg-io-source]
