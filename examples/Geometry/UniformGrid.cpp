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

int main(int, char**)
{
  constexpr size_t n = 4;
  constexpr Geometry::Polytope::Type g = Geometry::Polytope::Type::TriangularPrism;
  Mesh mesh;
  mesh = LocalMesh::UniformGrid(g, { n, n, n });
  mesh.getConnectivity().compute(3, 2);
  mesh.save("UniformGrid.mesh", IO::FileFormat::MEDIT);
  return 0;
}



