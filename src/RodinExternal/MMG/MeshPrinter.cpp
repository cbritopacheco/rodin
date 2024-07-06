/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"
#include "Rodin/IO/MeshPrinter.h"
#include "Rodin/IO/MEDIT.h"

#include "MeshPrinter.h"

namespace Rodin::External::MMG
{
  void MeshPrinter::print(std::ostream& os)
  {
    m_printer.print(os, false);
    printCorners(os);
    printRidges(os);
    printRequiredVertices(os);
    printRequiredEdges(os);
    m_printer.printEnd(os);
  }

  void MeshPrinter::printCorners(std::ostream& os)
  {
    const auto& mesh = getObject();
    os << IO::MEDIT::Keyword::Corners
       << '\n'
       << mesh.getCorners().size()
       << '\n';
    for (const auto& c : mesh.getCorners())
      os << c + 1 << '\n';
    os << '\n';
  }

  void MeshPrinter::printRidges(std::ostream& os)
  {
    const auto& mesh = getObject();
    os << IO::MEDIT::Keyword::Ridges
       << '\n'
       << mesh.getRidges().size()
       << '\n';
    for (const auto& r : mesh.getRidges())
      os << r + 1 << '\n';
    os << '\n';
  }

  void MeshPrinter::printRequiredVertices(std::ostream& os)
  {
    const auto& mesh = getObject();
    os << IO::MEDIT::Keyword::RequiredVertices
       << '\n'
       << mesh.getRequiredVertices().size()
       << '\n';
    for (const auto& r : mesh.getRequiredVertices())
      os << r + 1 << '\n';
    os << '\n';
  }

  void MeshPrinter::printRequiredEdges(std::ostream& os)
  {
    const auto& mesh = getObject();
    os << IO::MEDIT::Keyword::RequiredEdges
       << '\n'
       << mesh.getRequiredEdges().size()
       << '\n';
    for (const auto& r : mesh.getRequiredEdges())
      os << r + 1 << '\n';
    os << '\n';
  }
}

