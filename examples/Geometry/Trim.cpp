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

static constexpr size_t n = 6;
static constexpr Attribute trimAttribute = 2;

std::array rainbowAttribute = {1, 3, 4, 5, 6, 7, 8, 9, 10};

int main(int, char**)
{
  // Build the mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Geometry::Triangle, n, n);

  // We're gonna trim the center of the mesh.
  for (auto it = mesh.getElement(); !it.end(); ++it)
  {
    for (auto vit = it->getVertex(); !vit.end(); ++vit)
    {
      const auto v = *vit;
      if (2 <= v.x() && v.x() <= 3 && 2 <= v.y() && v.y() <= 3)
      {
        // Set the attribute to the trim attribute
        mesh.setAttribute({ it->getDimension(), it->getIndex() }, trimAttribute);
        break;
      }
      else
      {
        // Set other attributes to produce colors
        mesh.setAttribute(
            { it->getDimension(), it->getIndex() }, rainbowAttribute[it->getIndex() % rainbowAttribute.size()]);
      }
    }
  }

  // Perform the actual trimming of the attribute
  auto trimmed = mesh.trim(trimAttribute);
  trimmed.save("Trimmed.mfem.mesh", IO::FileFormat::MFEM);
  Alert::Info() << "Saved trimmed mesh to: Trimmed.mfem.mesh" << Alert::Raise;
}
