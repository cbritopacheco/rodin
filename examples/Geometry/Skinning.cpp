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

  Alert::Info() << "Skinning mesh..." << Alert::Raise;

  auto skin = mesh.skin();

  // skin.trace({
  //    {{3, 6}, 666},
  //    {{2, 6}, 777}
  //    });
  //
  skin.save("skin.mesh", IO::FileFormat::MEDIT);
  mesh.save("volume.mesh", IO::FileFormat::MEDIT);

  Alert::Info() << "Saved mesh to skin.mesh" << Alert::Raise;

  // skin.save("skin.mesh", IO::FileFormat::MEDIT);
}
