/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <boost/algorithm/string.hpp>

#include "Rodin/Alert.h"
#include "Rodin/Geometry.h"
#include "Rodin/IO/MEDIT.h"

#include "MeshLoader.h"

namespace Rodin::External::MMG
{
  void MeshLoader::load(std::istream& is)
  {
    Parent::load(is);
    is.clear();

    // At this point, the mesh topology has been built. It only remains to
    // track the geometric information (ie. edges, corners, ridges, etc.) that
    // MEDIT/MMG provides with its mesh format.
    auto& mesh = this->getObject();

    // We perform another pass on the file, seeking the geometric entities.
    std::string line;
    const auto& pos = getPositionMap();
    const auto& count = getCountMap();
    auto find = pos.find(IO::MEDIT::Keyword::Corners);
    if (find != pos.end())
    {
      is.seekg(find->second);
      for (size_t i = 0; i < count.at(IO::MEDIT::Keyword::Corners); i++)
      {
        std::getline(is, line);
        mesh.setCorner(std::stoul(line) - 1);
      }
    }

    find = pos.find(IO::MEDIT::Keyword::Ridges);
    if (find != pos.end())
    {
      is.seekg(find->second);
      for (size_t i = 0; i < count.at(IO::MEDIT::Keyword::Ridges); i++)
      {
        std::getline(is, line);
        mesh.setRidge(std::stoul(line) - 1);
      }
    }
  }
}

