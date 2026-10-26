/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iomanip>
#include <sstream>
#include <string>

#include "Rodin/Alert.h"
#include "Rodin/Geometry.h"
#include "Rodin/IO/MEDIT.h"

#include "MeshLoader.h"

namespace Rodin::MMG
{
  Index MeshLoader::readIndex(std::istream& is, IO::MEDIT::Keyword keyword)
  {
    std::string line;
    std::getline(is, line);
    std::istringstream entry(line);
    Index index = 0;
    entry >> std::ws;
    if (entry.peek() == '-' || !(entry >> index) || index == 0 ||
      !(entry >> std::ws).eof())
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Invalid one-based " << std::quoted(IO::MEDIT::toCharString(keyword))
        << " index " << std::quoted(line) << "." << Alert::Raise;
    }
    return index - 1;
  }

  void MeshLoader::load(std::istream& is)
  {
    Parent::load(is);
    is.clear();

    // At this point, the mesh topology has been built. It only remains to
    // track the geometric information (ie. edges, corners, ridges, etc.) that
    // MEDIT/MMG provides with its mesh format.
    auto& mesh = this->getObject();
    mesh.getCorners().clear();
    mesh.getRidges().clear();
    mesh.getRequiredVertices().clear();
    mesh.getRequiredEdges().clear();
    mesh.getRequiredTriangles().clear();
    mesh.getRequiredTetrahedra().clear();

    // We perform another pass on the file, seeking the geometric entities.
    const auto& pos = getPositionMap();
    const auto& count = getCountMap();
    auto find = pos.find(IO::MEDIT::Keyword::Corners);
    if (find != pos.end())
    {
      is.seekg(find->second);
      for (size_t i = 0; i < count.at(IO::MEDIT::Keyword::Corners); i++)
      {
        mesh.setCorner(readIndex(is, IO::MEDIT::Keyword::Corners));
      }
    }

    find = pos.find(IO::MEDIT::Keyword::Ridges);
    if (find != pos.end())
    {
      is.seekg(find->second);
      for (size_t i = 0; i < count.at(IO::MEDIT::Keyword::Ridges); i++)
      {
        mesh.setRidge(readIndex(is, IO::MEDIT::Keyword::Ridges));
      }
    }

    find = pos.find(IO::MEDIT::Keyword::RequiredEdges);
    if (find != pos.end())
    {
      is.seekg(find->second);
      for (size_t i = 0; i < count.at(IO::MEDIT::Keyword::RequiredEdges); i++)
      {
        mesh.setRequiredEdge(readIndex(is, IO::MEDIT::Keyword::RequiredEdges));
      }
    }

    find = pos.find(IO::MEDIT::Keyword::RequiredVertices);
    if (find != pos.end())
    {
      is.seekg(find->second);
      for (size_t i = 0; i < count.at(IO::MEDIT::Keyword::RequiredVertices); i++)
      {
        mesh.setRequiredVertex(readIndex(is, IO::MEDIT::Keyword::RequiredVertices));
      }
    }

    find = pos.find(IO::MEDIT::Keyword::RequiredTriangles);
    if (find != pos.end())
    {
      is.seekg(find->second);
      for (size_t i = 0; i < count.at(IO::MEDIT::Keyword::RequiredTriangles); i++)
      {
        mesh.setRequiredTriangle(readIndex(is, IO::MEDIT::Keyword::RequiredTriangles));
      }
    }

    find = pos.find(IO::MEDIT::Keyword::RequiredTetrahedra);
    if (find != pos.end())
    {
      is.seekg(find->second);
      for (size_t i = 0; i < count.at(IO::MEDIT::Keyword::RequiredTetrahedra); i++)
      {
        mesh.setRequiredTetrahedron(
          readIndex(is, IO::MEDIT::Keyword::RequiredTetrahedra));
      }
    }
  }
}
