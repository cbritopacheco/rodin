/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [interface-source]
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  //! [interface-setup]
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(1, 2);
  mesh.setAttribute({2, 0}, 1);
  mesh.setAttribute({2, 1}, 2);
  P0 vh(mesh);
  GridFunction field(vh);
  field = RealFunction([](const Point& point) {
    return point.getPolytope().getAttribute() == 1 ? 1.0 : 3.0;
  });
  //! [interface-setup]

  //! [interface-evaluate]
  const auto jump = Jump(field);
  const auto average = Average(field);

  constexpr Real tolerance = 1e-12;
  bool valid = true;
  size_t interfaces = 0;
  std::ofstream samples("InterfaceTraces.csv");
  samples << "x,y,first,second,jump,average\n";
  for (auto face = mesh.getInterface(); !face.end(); ++face)
  {
    const Point point(*face, Math::SpatialPoint{0.5});
    const auto firstCell = mesh.getCell(0);
    const auto secondCell = mesh.getCell(1);
    Math::SpatialPoint firstReference, secondReference;
    firstCell->getTransformation().inverse(firstReference, point.getPhysicalCoordinates());
    secondCell->getTransformation().inverse(secondReference, point.getPhysicalCoordinates());
    const Real first = field(Point(*firstCell, firstReference));
    const Real second = field(Point(*secondCell, secondReference));
    const Real value = jump(point);
    const Real mean = average(point);
    samples << point.x() << ',' << point.y() << ',' << first << ','
            << second << ',' << value << ',' << mean << '\n';
    std::cout << "Traces: " << first << ", " << second
              << "; jump: " << value << "; average: " << mean << std::endl;
    valid = valid && std::abs(first - 1) < tolerance
      && std::abs(second - 3) < tolerance
      && std::abs(std::abs(value) - 2) < tolerance && std::abs(mean - 2) < tolerance;
    ++interfaces;
  }
  //! [interface-evaluate]
  samples.flush();
  return valid && interfaces == 1 && samples ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [interface-source]
