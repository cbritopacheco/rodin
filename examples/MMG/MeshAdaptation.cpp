/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [mmg-adapt-source]
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <Rodin/Geometry.h>
#include <Rodin/MMG.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  //! [mmg-adapt-setup]
  MMG::Mesh mesh;
  mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {9, 9});
  mesh.scale(1.0 / 8.0);
  mesh.getConnectivity().compute(1, 2);
  for (const Index corner : {0, 8, 72, 80})
    mesh.setCorner(corner);
  //! [mmg-adapt-setup]

  bool valid = true;
  for (size_t stage = 0; stage < 2; ++stage)
  {
    if (stage == 1)
    {
      //! [mmg-adapt-call]
      P1 vh(mesh);
      MMG::RealGridFunction sizeMap(vh);
      sizeMap = [](const Point& point) { return 0.05 + 0.15 * point.x(); };
      MMG::Adapt().setHMin(0.04).setHMax(0.21).adapt(mesh, sizeMap);
      //! [mmg-adapt-call]
    }
    const std::string stem = stage == 0 ? "Before" : "After";
    std::ofstream triangles(stem + ".csv");
    triangles.precision(17);
    triangles << "cell,x0,y0,x1,y1,x2,y2\n";
    for (auto cell = mesh.getCell(); !cell.end(); ++cell)
    {
      triangles << cell->getIndex();
      for (auto vertex = cell->getVertex(); !vertex.end(); ++vertex)
        triangles << ',' << vertex->getCoordinates().x() << ',' << vertex->getCoordinates().y();
      triangles << '\n';
      valid = valid && cell->getMeasure() > 0;
    }
    triangles.flush();
    constexpr Real areaTolerance = 1e-10;
    valid = valid && bool(triangles) && !mesh.isEmpty()
      && std::abs(mesh.getArea() - 1) < areaTolerance;
    mesh.save(stem + ".medit.mesh", IO::FileFormat::MEDIT);
    std::cout << stem << ": " << mesh.getVertexCount() << " vertices, "
              << mesh.getCellCount() << " triangles, area " << mesh.getArea() << std::endl;
  }
  return valid ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [mmg-adapt-source]
