/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Geometry/Polytope.h"
#include <Rodin/Geometry.h>
#include <Rodin/Alert.h>

using namespace Rodin;
using namespace Geometry;

int main(int, char**)
{
  constexpr size_t n = 4;
  Mesh mesh;

  mesh = mesh.Box(Polytope::Type::Point, { n });
  mesh.save("Box_Point.mesh", IO::FileFormat::MEDIT);

  mesh = mesh.Box(Polytope::Type::Segment, { n, n });
  mesh.save("Box_Segment.mesh", IO::FileFormat::MEDIT);

  mesh = mesh.Box(Polytope::Type::Triangle, { n, n, n });
  mesh.save("Box_Triangle.mesh", IO::FileFormat::MEDIT);

  mesh = mesh.Box(Polytope::Type::Quadrilateral, { n, n, n });
  mesh.save("Box_Quadrilateral.mesh", IO::FileFormat::MEDIT);

  return 0;
}




