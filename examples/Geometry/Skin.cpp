/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>
#include <Rodin/Alert.h>

using namespace Rodin;
using namespace Geometry;

const char* meshFile = "../resources/examples/Geometry/Skinning.mesh";

int main(int, char**)
{
  Mesh mesh;
  mesh.load(meshFile);

  // Needed for computing the boundary
  mesh.getConnectivity().compute(2, 3);

  // Optional: for computing the edges on the boundary
  mesh.getConnectivity().compute(2, 1);

  std::cout << mesh.getCount(1) << std::endl;

  for (auto it = mesh.getFace(); !it.end(); ++it)
    mesh.setAttribute({ it->getDimension(), it->getIndex() }, it->getIndex() + 1);

  const size_t edgeDimension = 1;
  for (auto it = mesh.getPolytope(edgeDimension); !it.end(); ++it)
    mesh.setAttribute({ it->getDimension(), it->getIndex() }, it->getIndex() + 1);

  auto skin = mesh.skin();

  skin.save("Skin.mesh", IO::FileFormat::MEDIT);
}