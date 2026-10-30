/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [mesh-basics-source]
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main()
{
  //! [mesh-basics-build]
  Mesh mesh = Mesh<Context::Local>::Builder()
                .initialize(2)
                .nodes(6)
                .vertex({0, 0})
                .vertex({1, 0})
                .vertex({2, 0})
                .vertex({0, 1})
                .vertex({1, 1})
                .vertex({2, 1})
                .polytope(Polytope::Type::Triangle, {0, 1, 4})
                .attribute({2, 0}, 1)
                .polytope(Polytope::Type::Triangle, {0, 4, 3})
                .attribute({2, 1}, 1)
                .polytope(Polytope::Type::Quadrilateral, {1, 2, 5, 4})
                .attribute({2, 2}, 2)
                .finalize();
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 2);
  //! [mesh-basics-build]

  //! [mesh-basics-queries]
  const Real area = mesh.getArea();
  const Real perimeter = mesh.getPerimeter();
  std::cout << "Vertices: " << mesh.getVertexCount() << '\n'
            << "Edges: " << mesh.getFaceCount() << '\n'
            << "Cells: " << mesh.getCellCount() << '\n'
            << "Area: " << area << ", perimeter: " << perimeter << '\n'
            << "Region 1 area: " << mesh.getArea(1) << '\n'
            << "Region 2 area: " << mesh.getArea(2) << std::endl;
  //! [mesh-basics-queries]

  //! [mesh-basics-iterate]
  for (auto cell = mesh.getCell(); !cell.end(); ++cell)
  {
    std::cout << "Cell " << cell->getIndex();
    if (const auto attribute = cell->getAttribute())
      std::cout << ": attribute " << *attribute;
    std::cout << ", measure " << cell->getMeasure() << '\n';
  }
  size_t boundaryCount = 0;
  for (auto face = mesh.getBoundary(); !face.end(); ++face)
    ++boundaryCount;
  //! [mesh-basics-iterate]

  //! [mesh-basics-extract]
  auto region = mesh.keep(1);
  auto boundary = mesh.skin();
  const auto components = mesh.ccl([](const Polytope& first, const Polytope& second) {
    return first.getAttribute() == second.getAttribute();
  });
  std::cout << "Region cells: " << region.getCellCount() << '\n'
            << "Boundary cells: " << boundary.getCellCount() << '\n'
            << "Attribute components: " << components.getComponents().size() << std::endl;
  //! [mesh-basics-extract]

  //! [mesh-basics-partition]
  BalancedCompactPartitioner partitioner(mesh);
  partitioner.partition(2);
  for (auto cell = mesh.getCell(); !cell.end(); ++cell)
  {
    std::cout << "Cell " << cell->getIndex() << ": partition "
              << partitioner.getPartition(cell->getIndex()) << '\n';
  }
  //! [mesh-basics-partition]

  //! [mesh-basics-box]
  Mesh surface = Mesh().Box(Polytope::Type::Triangle, {3, 3, 3});
  surface.scale(1.0 / 3.0);
  std::cout << "Box boundary: dimension " << surface.getDimension()
            << ", embedding dimension " << surface.getSpaceDimension() << ", area "
            << surface.getArea() << std::endl;
  //! [mesh-basics-box]

  std::ofstream edges("MeshBasics.csv");
  edges << "edge,x0,y0,x1,y1,boundary\n";
  for (auto edge = mesh.getFace(); !edge.end(); ++edge)
  {
    edges << edge->getIndex();
    for (auto vertex = edge->getVertex(); !vertex.end(); ++vertex)
      edges << ',' << vertex->getCoordinates().x() << ',' << vertex->getCoordinates().y();
    edges << ',' << mesh.isBoundary(edge->getIndex()) << '\n';
  }
  mesh.save("MeshBasics.mesh", IO::FileFormat::MFEM);
  constexpr Real tolerance = 1e-12;
  const bool valid = mesh.getVertexCount() == 6 && mesh.getFaceCount() == 8 &&
    mesh.getCellCount() == 3 && boundaryCount == 6 && region.getCellCount() == 2 &&
    boundary.getCellCount() == 6 && components.getComponents().size() == 2 &&
    std::abs(area - 2) < tolerance && std::abs(perimeter - 6) < tolerance &&
    std::abs(mesh.getArea(1) - 1) < tolerance &&
    std::abs(mesh.getArea(2) - 1) < tolerance && surface.getDimension() == 2 &&
    surface.getSpaceDimension() == 3 && std::abs(surface.getArea() - 6) < tolerance;
  edges.flush();
  return valid && edges ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [mesh-basics-source]
