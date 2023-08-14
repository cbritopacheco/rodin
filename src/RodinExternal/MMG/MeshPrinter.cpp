/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"
#include "Rodin/IO/MeshPrinter.h"
#include "Rodin/IO/Helpers.h"

#include "MeshPrinter.h"

namespace Rodin::External::MMG
{
  void MeshPrinter::print(std::ostream& os)
  {
    const auto& mesh = getObject();

    IO::MeshPrinter<IO::FileFormat::MEDIT, Context::Serial> printer(mesh);
    printer.footer(false).print(os);

    // Print edges
    if (mesh.getDimension() == 3)
    {
      os << '\n'
        << IO::Medit::Keyword::Edges
        << '\n'
        << mesh.getEdges().size()
        << '\n';

      for (const auto& e : mesh.getEdges())
        os << e.endpoints.first + 1 << " " << e.endpoints.second + 1 << " " << e.ref << '\n';
    }

    // Print corners
    os << '\n'
      << IO::Medit::Keyword::Corners
      << '\n'
      << mesh.getCorners().size()
      << '\n';
    for (const auto& c : mesh.getCorners())
      os << c + 1 << '\n';

    // Print ridges
    os << '\n'
      << IO::Medit::Keyword::Ridges
      << '\n'
      << mesh.getRidges().size()
      << '\n';
    for (const auto& r : mesh.getRidges())
      os << r + 1 << '\n';

    // Print footer
    os << "\n\n" << IO::Medit::Keyword::End;
  }
}

