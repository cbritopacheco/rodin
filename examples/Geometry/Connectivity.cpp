/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [connectivity-source]
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main()
{
  //! [connectivity-mesh]
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto& connectivity = mesh.getConnectivity();
  connectivity.compute(1, 2);
  connectivity.compute(2, 1);
  //! [connectivity-mesh]

  //! [connectivity-incidence]
  size_t boundaryEdges = 0;
  for (auto edge = mesh.getFace(); !edge.end(); ++edge)
  {
    const auto& cells = connectivity.getIncidence({1, 2}, edge->getIndex());
    std::cout << "Edge " << edge->getIndex() << ": " << cells.size()
              << " incident cell(s)";
    if (edge->isBoundary())
    {
      ++boundaryEdges;
      std::cout << " [boundary]";
    }
    std::cout << std::endl;
  }
  std::cout << "Vertices: " << mesh.getVertexCount()
            << ", edges: " << connectivity.getCount(1)
            << ", cells: " << mesh.getCellCount() << ", boundary edges: " << boundaryEdges
            << std::endl;
  //! [connectivity-incidence]

  // Export the actual topology for the tutorial's labeled diagram.
  std::ofstream vertices("Connectivity.vertices.csv");
  vertices.precision(17);
  vertices << "id,x,y\n";
  for (auto vertex = mesh.getVertex(); !vertex.end(); ++vertex)
  {
    vertices << vertex->getIndex() << ',' << vertex->getCoordinates().x() << ','
             << vertex->getCoordinates().y() << '\n';
  }
  std::ofstream edges("Connectivity.edges.csv");
  edges << "id,v0,v1,cells\n";
  for (auto edge = mesh.getFace(); !edge.end(); ++edge)
  {
    const auto& key = connectivity.getPolytope(1, edge->getIndex());
    edges << edge->getIndex() << ',' << key[0] << ',' << key[1] << ','
          << connectivity.getIncidence({1, 2}, edge->getIndex()).size() << '\n';
  }
  vertices.flush();
  edges.flush();
  return vertices && edges ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [connectivity-source]
