/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Alert.h>
#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main(int, char**)
{
  const char* meshFile = "../resources/examples/Geometry/Skinning.mesh";
  Mesh mesh;
  mesh.load(meshFile);
  mesh.getConnectivity().compute(2, 3);
  // mesh.getConnectivity().compute(2, 3);
  auto d = mesh.skin();
  d.getConnectivity().compute(1, 2);
  d.trace(
      {
        { { 2, 6 }, 9 },
        { { 3, 6 }, 3 }
      });

  d.save("miaow.mesh", IO::FileFormat::MEDIT);

  return 0;
}

