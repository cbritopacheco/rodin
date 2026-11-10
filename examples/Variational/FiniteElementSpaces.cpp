/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [spaces-source]
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
  //! [spaces-setup]
  constexpr size_t nodes = 8;
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  mesh.scale(1.0 / (nodes - 1));
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 1);
  P0 cellSpace(mesh);
  P0g constantSpace(mesh);
  P1 linearSpace(mesh);
  H1 quadraticSpace(std::integral_constant<size_t, 2>{}, mesh);
  P1 vectorSpace(mesh, 2);
  P1<Complex> complexSpace(mesh);
  //! [spaces-setup]

  //! [spaces-counts]
  std::cout << "Vertices: " << mesh.getVertexCount() << '\n'
            << "Triangles: " << mesh.getCellCount() << '\n'
            << "P0: " << cellSpace.getSize() << '\n'
            << "P0g: " << constantSpace.getSize() << '\n'
            << "P1: " << linearSpace.getSize() << '\n'
            << "H1 degree 2: " << quadraticSpace.getSize() << '\n'
            << "Vector P1: " << vectorSpace.getSize() << '\n'
            << "Complex P1: " << complexSpace.getSize() << std::endl;
  //! [spaces-counts]

  //! [spaces-fields]
  RealFunction exact = [](const Point& point) {
    return point.x() * point.x() + point.y() * point.y();
  };
  GridFunction constant(cellSpace);
  GridFunction linear(linearSpace);
  GridFunction quadratic(quadraticSpace);
  constant = exact;
  linear = exact;
  quadratic = exact;
  //! [spaces-fields]

  //! [spaces-sampling]
  // Refine only the sampling grid in each reference triangle, not the FE mesh.
  constexpr size_t subdivisions = 6;
  std::ofstream samples("FiniteElementSpaces.csv");
  samples.precision(17);
  samples << "cell,i,j,x,y,exact,p0,p1,p2\n";
  for (auto cell = mesh.getCell(); !cell.end(); ++cell)
  {
    for (size_t i = 0; i <= subdivisions; ++i)
    {
      for (size_t j = 0; j <= subdivisions - i; ++j)
      {
        const Point point(
          *cell, Math::SpatialPoint{Real(i) / subdivisions, Real(j) / subdivisions});
        samples << cell->getIndex() << ',' << i << ',' << j << ',' << point.x() << ','
                << point.y() << ',' << exact(point) << ',' << constant(point) << ','
                << linear(point) << ',' << quadratic(point) << '\n';
      }
    }
  }
  //! [spaces-sampling]

  //! [spaces-reference]
  P1Element<Real> element(Polytope::Type::Triangle);
  std::ofstream basis("ReferenceP1.csv");
  basis.precision(17);
  basis << "i,j,r,s,phi0,phi1,phi2\n";
  for (size_t i = 0; i <= subdivisions; ++i)
  {
    for (size_t j = 0; j <= subdivisions - i; ++j)
    {
      const Math::SpatialPoint point{Real(i) / subdivisions, Real(j) / subdivisions};
      basis << i << ',' << j << ',' << point.x() << ',' << point.y();
      for (size_t k = 0; k < element.getCount(); ++k)
        basis << ',' << element.getBasis(k)(point);
      basis << '\n';
    }
  }
  //! [spaces-reference]
  samples.flush();
  basis.flush();
  return samples && basis ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [spaces-source]
