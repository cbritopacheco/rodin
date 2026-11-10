/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main()
{
  constexpr size_t nodes = 8;
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  mesh.scale(1.0 / (nodes - 1));
  mesh.getConnectivity().compute(1, 2);
  const Real area = mesh.getArea();
  const Real perimeter = mesh.getPerimeter();
  std::cout << "Area: " << area << '\n' << "Perimeter: " << perimeter << std::endl;
  constexpr Real tolerance = 1e-12;
  return std::abs(area - 1) < tolerance && std::abs(perimeter - 4) < tolerance
    ? EXIT_SUCCESS
    : EXIT_FAILURE;
}
