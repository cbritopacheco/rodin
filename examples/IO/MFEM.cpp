/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [mfem-source]
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  //! [mfem-mesh]
  Mesh mesh = Mesh<Context::Local>::Builder()
    .initialize(2)
    .nodes(4)
    .vertex({0, 0})
    .vertex({1, 0})
    .vertex({0, 1})
    .vertex({1, 1})
    .polytope(Polytope::Type::Triangle, {0, 1, 2})
    .polytope(Polytope::Type::Triangle, {1, 3, 2})
    .finalize();
  //! [mfem-mesh]

  //! [mfem-save]
  P1 vh(mesh);
  GridFunction u(vh);
  u = RealFunction([](const Point& point) { return point.x() + point.y(); });
  mesh.save("Square.mfem.mesh", IO::FileFormat::MFEM);
  u.save("Square.mfem.gf", IO::FileFormat::MFEM);
  //! [mfem-save]

  //! [mfem-load]
  Mesh loaded;
  loaded.load("Square.mfem.mesh", IO::FileFormat::MFEM);
  P1 loadedSpace(loaded);
  GridFunction restored(loadedSpace);
  restored.load("Square.mfem.gf", IO::FileFormat::MFEM);
  //! [mfem-load]

  //! [mfem-check]
  Real error = 0;
  for (auto cell = loaded.getCell(); !cell.end(); ++cell)
  {
    const Point point(*cell, Math::SpatialPoint{1.0 / 3.0, 1.0 / 3.0});
    error = std::max(error, std::abs(restored(point) - point.x() - point.y()));
  }
  constexpr Real tolerance = 1e-12;
  const bool valid =
    loaded.getVertexCount() == 4 && loaded.getCellCount() == 2 && error < tolerance;
  std::cout << "Maximum cell-centroid error: " << error << std::endl;
  //! [mfem-check]

  IO::XDMF xdmf("Square");
  xdmf.grid().setMesh(loaded).add("u", restored);
  xdmf.write().flush();
  return valid ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [mfem-source]
