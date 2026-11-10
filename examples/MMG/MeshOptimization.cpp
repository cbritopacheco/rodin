/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [mmg-optimizer-source]
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>

#include <Rodin/Geometry.h>
#include <Rodin/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main()
{
  //! [mmg-optimizer-setup]
  constexpr size_t nodes = 9;
  MMG::Mesh mesh;
  mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  mesh.scale(1.0 / (nodes - 1));
  mesh.getConnectivity().compute(1, 2);
  for (const Index corner :
    {size_t{0}, nodes - 1, nodes * (nodes - 1), nodes * nodes - 1})
    mesh.setCorner(corner);
  for (auto vertex = mesh.getVertex(); !vertex.end(); ++vertex)
  {
    auto coordinates = vertex->getCoordinates();
    if (coordinates.x() > 0 && coordinates.x() < 1 && coordinates.y() > 0 &&
      coordinates.y() < 1)
    {
      // A quarter-grid displacement creates a visible quality-improvement task.
      coordinates.x() += 0.03 * std::sin(17 * coordinates.y());
      coordinates.y() += 0.03 * std::sin(13 * coordinates.x());
      mesh.setVertexCoordinates(vertex->getIndex(), coordinates);
    }
  }
  mesh.flush();
  //! [mmg-optimizer-setup]

  bool valid = true;
  for (size_t stage = 0; stage < 2; ++stage)
  {
    if (stage == 1)
    {
      //! [mmg-optimizer-call]
      MMG::Optimizer().setHMin(0.05).setHMax(0.2).optimize(mesh);
      //! [mmg-optimizer-call]
    }
    const std::string stem = stage == 0 ? "Before" : "After";
    std::ofstream triangles(stem + ".csv");
    triangles.precision(17);
    triangles << "cell,x0,y0,x1,y1,x2,y2,quality\n";
    Real minimumQuality = std::numeric_limits<Real>::max();
    for (auto cell = mesh.getCell(); !cell.end(); ++cell)
    {
      Math::SpatialPoint points[3];
      size_t i = 0;
      for (auto vertex = cell->getVertex(); !vertex.end(); ++vertex)
        points[i++] = vertex->getCoordinates();
      Real squaredEdges = 0;
      for (size_t edge = 0; edge < 3; ++edge)
        squaredEdges += (points[edge] - points[(edge + 1) % 3]).squaredNorm();
      const Real quality = 4 * std::sqrt(3.0) * cell->getMeasure() / squaredEdges;
      minimumQuality = std::min(minimumQuality, quality);
      triangles << cell->getIndex();
      for (const auto& point : points)
        triangles << ',' << point.x() << ',' << point.y();
      triangles << ',' << quality << '\n';
    }
    constexpr Real areaTolerance = 1e-10;
    valid = valid && !mesh.isEmpty() && minimumQuality > 0 &&
      std::abs(mesh.getArea() - 1) < areaTolerance;
    triangles.flush();
    valid = valid && bool(triangles);
    mesh.save(stem + ".medit.mesh", IO::FileFormat::MEDIT);
    std::cout << stem << ": " << mesh.getCellCount() << " triangles, area "
              << mesh.getArea() << ", minimum quality " << minimumQuality << std::endl;
  }
  return valid ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [mmg-optimizer-source]
