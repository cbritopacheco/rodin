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
  const size_t dim = 2;

  Mesh mesh;
  mesh.load("../resources/mfem/StarSquare.mfem.mesh", IO::FileFormat::MFEM);
  // mesh.initialize(dim, dim)
  //     .vertex({0, 0})
  //     .vertex({1, 0})
  //     .vertex({0, 1})
  //     .vertex({1, 1})
  //     .element(Geometry::Type::Triangle, {0, 1, 2})
  //     .element(Geometry::Type::Triangle, {2, 1, 3})
  //     .finalize();

  mesh.save("miaow.mesh", IO::FileFormat::MEDIT);

  return 0;
}


