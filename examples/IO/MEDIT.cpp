/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [medit-source]
#include <cstdlib>
#include <iostream>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main()
{
  //! [medit-write]
  Mesh square = Mesh().UniformGrid(Polytope::Type::Triangle, {5, 5});
  square.scale(0.25);
  square.save("Square.medit.mesh", IO::FileFormat::MEDIT);
  //! [medit-write]

  //! [medit-convert]
  Mesh loaded;
  loaded.load("Square.medit.mesh", IO::FileFormat::MEDIT);
  loaded.save("Square.mfem.mesh", IO::FileFormat::MFEM);

  Mesh converted;
  converted.load("Square.mfem.mesh", IO::FileFormat::MFEM);
  //! [medit-convert]

  const bool valid = loaded.getVertexCount() == square.getVertexCount() &&
    loaded.getCellCount() == square.getCellCount() &&
    converted.getVertexCount() == square.getVertexCount() &&
    converted.getCellCount() == square.getCellCount();
  std::cout << "Vertices: " << converted.getVertexCount()
            << ", triangles: " << converted.getCellCount() << '\n'
            << "Mesh count round trip: " << (valid ? "passed" : "failed") << std::endl;
  return valid ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [medit-source]
