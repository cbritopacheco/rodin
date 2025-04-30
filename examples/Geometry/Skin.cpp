/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>
#include <Rodin/Alert.h>

using namespace Rodin;
using namespace Rodin::Geometry;

static constexpr Attribute trimAttribute = 2;

int main(int, char**)
{
  size_t n = 8;
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.getConnectivity().compute(1, 0);
  for (auto it = mesh.getFace(); !it.end(); ++it)
    mesh.setAttribute({ it->getDimension(), it->getIndex() }, it->getIndex() + 1);
  mesh.save("Grid.mesh", IO::FileFormat::MEDIT);
  auto skin = mesh.skin();
  skin.save("Skin.mesh", IO::FileFormat::MEDIT);
}
