/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <boost/algorithm/string.hpp>

#include "Rodin/Alert.h"
#include "Rodin/Geometry.h"
#include "Rodin/IO/MeshLoader.h"

#include "MeshLoader.h"

namespace Rodin::External::MMG
{
  void MeshLoader::load(std::istream& is)
  {
    IO::MeshLoader<IO::FileFormat::MEDIT, Context::Serial> loader(getObject());
    loader.load(is);
    // At this point, the mesh topology has been built. It only remains to
    // track the geometric information (ie. edges, corners, ridges, etc.) that
    // MEDIT/MMG provides with its mesh format.

    // We perform another pass on the file, seeking the geometric entities.
    is.clear();
    std::string line;
    for (const auto& [_, kw] : IO::Medit::KeywordMap)
    {
      if (!loader.getContext().pos.count(kw))
        continue;
      auto g = loader.getContext().pos.at(kw);
      if (!g.has_value())
        continue;
      is.seekg(*g);

      if (!loader.getContext().count.count(kw))
      {
        Alert::Exception()
          << "Failed to parse \"" << kw << "\" count."
          << Alert::Raise;
      }

      switch (kw)
      {
        case IO::Medit::Keyword::Edges:
        {
          if (getObject().getDimension() >= 3)
          {
            for (size_t i = 0; i < loader.getContext().count.at(kw); i++)
            {
              if(!std::getline(is, line))
                Alert::Exception("Bad mesh format.").raise();
              std::istringstream lss(line);
              int v1, v2, ref;
              lss >> v1 >> v2 >> ref;
              getObject().edge({v1 - 1, v2 - 1}, ref);
            }
          }
          break;
        }
        case IO::Medit::Keyword::Corners:
        {
          for (size_t i = 0; i < loader.getContext().count.at(kw); i++)
          {
            if(!std::getline(is, line))
              Alert::Exception("Bad mesh format.").raise();
            std::istringstream lss(line);
            int vertexIdx;
            lss >> vertexIdx;
            getObject().corner(vertexIdx - 1);
          }
          break;
        }
        case IO::Medit::Keyword::Ridges:
        {
          for (size_t i = 0; i < loader.getContext().count.at(kw); i++)
          {
            if(!std::getline(is, line))
              Alert::Exception("Bad mesh format.").raise();
            std::istringstream lss(line);
            int edgeIdx;
            lss >> edgeIdx;
            getObject().ridge(edgeIdx - 1);
          }
          break;
        }
        default:
        {
          // No action since we should have already parsed the syntax and
          // verified it is correct.
          break;
        }
      }
    }
  }
}

