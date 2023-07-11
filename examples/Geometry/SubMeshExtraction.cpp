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
  size_t n = 6;
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Geometry::Triangle, n, n);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  for (auto it = mesh.getElement(); !it.end(); ++it)
  {
    for (auto vit = it->getVertex(); !vit.end(); ++vit)
    {
      const auto v = *vit;
      if (2 <= v.x() && v.x() <= 3 && 2 <= v.y() && v.y() <= 3)
      {
        mesh.setAttribute({ it->getDimension(), it->getIndex() }, trimAttribute);
        break;
      }
    }
  }

  for (auto it = mesh.getFace(); !it.end(); ++it)
    mesh.setAttribute({ it->getDimension(), it->getIndex() }, it->getIndex() + 1);

  auto trimmed = mesh.trim(trimAttribute);
  auto skin = mesh.skin();

  skin.save("Skin.mesh", IO::FileFormat::MEDIT);
  trimmed.save("Trimmed.mesh", IO::FileFormat::MEDIT);
}

