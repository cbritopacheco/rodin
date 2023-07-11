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
  mesh.getConnectivity().compute(2, 3);
  mesh.getConnectivity().compute(2, 1);
  std::cout << mesh.getCount(1) << std::endl;
  // mesh.getConnectivity().compute(1, 3);

  for (auto it = mesh.getFace(); !it.end(); ++it)
    mesh.setAttribute({ it->getDimension(), it->getIndex() }, it->getIndex() + 1);

  for (auto it = mesh.getPolytope(1); !it.end(); ++it)
    mesh.setAttribute({ it->getDimension(), it->getIndex() }, it->getIndex() + 1);

  auto skin = mesh.skin();

  skin.save("Skin.mesh", IO::FileFormat::MEDIT);
}
