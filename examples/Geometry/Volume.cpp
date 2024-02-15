/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  size_t n = 32;
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.getConnectivity().compute(1, 2);

  std::cout << "Volume: " << mesh.getVolume() << std::endl;
  std::cout << "Perimeter: " << mesh.getPerimeter() << std::endl;
}




