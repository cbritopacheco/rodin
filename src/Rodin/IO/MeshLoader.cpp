/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <boost/algorithm/string.hpp>

#include "Rodin/Alert.h"

#include "MeshLoader.h"

namespace Rodin::IO
{
  void MeshLoader<FileFormat::MFEM, Context::Serial>::load(std::istream& is)
  {
    assert(is);
    auto fmt = getMeshFormat(is);
    is.clear();
    is.seekg(0, std::ios::beg);
    if (!fmt.has_value())
      Alert::Exception("Unrecognized mesh format.").raise();
    if (fmt == FileFormat::MFEM)
    {
      getObject() =
        Rodin::Geometry::Mesh<Context::Serial>(
            mfem::Mesh(is, 0, 1, getFixOrientation()));
    }
    else
    {
      Alert::Exception("Cannot determine MFEM format version.").raise();
    }
  }

  void MeshLoader<FileFormat::GMSH, Context::Serial>::load(std::istream& is)
  {
    assert(is);
    auto fmt = getMeshFormat(is);
    is.clear();
    is.seekg(0, std::ios::beg);
    if (!fmt.has_value())
      Alert::Exception("Unrecognized mesh format.").raise();
    if (fmt == FileFormat::GMSH)
    {
      getObject() = Rodin::Geometry::Mesh<Context::Serial>(mfem::Mesh(is, 0, 1, getFixOrientation()));
    }
    else
    {
      Alert::Exception("Cannot determine GMSH format version.").raise();
    }
  }

  void MeshLoader<FileFormat::MEDIT, Context::Serial>::load(std::istream& is)
  {
    std::string line;
    while (std::getline(is, line))
    {
      boost::algorithm::trim(line);
      if (!line.empty())
        break;
    }

    std::istringstream iss(line);
    std::string kw;
    iss >> kw;

    int version = 0;
    if (Medit::KeywordMap.left.count(kw))
    {
      switch (Medit::KeywordMap.left.at(kw))
      {
        case Medit::Keyword::MeshVersionFormatted:
        {
          if (!(iss >> version)) // Version is not on this line
          {
            // Try next line
            std::getline(is, line);
            if (!(std::istringstream(line) >> version))
            {
              Alert::Exception()
                << "Bad mesh format. Encountered \""
                << line << "\" after keyword \""
                << Medit::Keyword::MeshVersionFormatted
                << "\"."
                << Alert::Raise;
            }
          }
          break;
        }
        default:
        {
          Alert::Exception()
            << "Bad mesh format. Encountered keyword \""
            << kw
            << "\" before keyword \""
            << Medit::Keyword::MeshVersionFormatted
            << "\"."
            << Alert::Raise;
          break;
        }
      }
    }
    else
    {
      Alert::Exception()
        << "Bad mesh format. Encountered \""
        << line << "\" before keyword \""
        << Medit::Keyword::MeshVersionFormatted
        << "\"."
        << Alert::Raise;
    }

    assert(version >= 2 && version <= 3);

    Rodin::Geometry::Mesh<Context::Serial> mesh;
    Rodin::Geometry::Mesh<Context::Serial>::Builder build;

    int spaceDim = 0;
    while (std::getline(is, line))
    {
      boost::algorithm::trim(line);
      if (line.empty())
        continue;

      std::istringstream iss(line);
      std::string kw;
      iss >> kw;
      if (kw == Medit::KeywordMap.right.at(Medit::Keyword::Dimension))
      {
        if (!(iss >> spaceDim)) // Dimension not on this line
        {
          // Try next line
          std::getline(is, line);
          if (!(std::istringstream(line) >> spaceDim))
          {
            Alert::Exception()
              << "Bad solution format. Encountered \""
              << line << "\" after keyword \""
              << Medit::Keyword::Dimension
              << "\"."
              << Alert::Raise;
          }
        }
        if (spaceDim < 2 || spaceDim > 3)
          Alert::Exception() << "Invalid dimension " << spaceDim << Alert::Raise;
        break;
      }
      else
      {
        Alert::Exception()
          << "Bad solution format. Encountered keyword \"" << kw
          << "\" before keyword \""
          << Medit::Keyword::Dimension
          << "\"."
          << Alert::Raise;
        return;
      }
    }

    while (std::getline(is, line))
    {
      boost::algorithm::trim(line);
      if (line.empty())
        continue;

      std::string kw;
      std::istringstream(line) >> kw;

      if (!std::isalpha(kw[0]))
        continue;

      if (Medit::KeywordMap.left.count(kw))
      {
        switch (Medit::KeywordMap.left.at(kw))
        {
          case Medit::Keyword::Corners:
          case Medit::Keyword::Ridges:
          case Medit::Keyword::Tetrahedra:
          case Medit::Keyword::Triangles:
          case Medit::Keyword::Edges:
          case Medit::Keyword::Vertices:
          {
            auto ent = Medit::KeywordMap.left.at(kw);
            std::getline(is, line);
            m_context.pos[ent] = is.tellg();
            std::istringstream lss(line);
            lss >> m_context.count[ent];
            continue; // This is to break out of the for loop
          }
          case Medit::Keyword::End:
          {
            break;
          }
          default:
          {
            Alert::Warning()
              << "Ignoring unrecognized keyword \"" << kw << "\"."
              << Alert::Raise;
            break;
          }
        }
        break; // If we reach here we should break out of the loop
      }
      else
      {
        Alert::Warning()
          << "Ignoring unrecognized keyword \"" << kw << "\"."
          << Alert::Raise;
        continue;
      }
    }

    // Set entities with zero count
    for (const auto& [_, kw] : Medit::KeywordMap)
    {
      switch (kw)
      {
        case Medit::Keyword::Corners:
        case Medit::Keyword::Ridges:
        case Medit::Keyword::Tetrahedra:
        case Medit::Keyword::Triangles:
        case Medit::Keyword::Edges:
        case Medit::Keyword::Vertices:
        {
          if (!m_context.count.count(kw))
            m_context.count[kw] = 0;
          break;
        }
        default:
        {
          // No-op
          break;
        }
      }
    }

    // Infer type of mesh
    bool isSurfaceMesh;
    if (spaceDim == 3 && m_context.count.at(Medit::Keyword::Tetrahedra) == 0)
    {
      build = mesh.initialize(spaceDim - 1, spaceDim);
      isSurfaceMesh = true;
    }
    else if (spaceDim == 2 && m_context.count.at(Medit::Keyword::Triangles) == 0)
    {
      build = mesh.initialize(spaceDim - 1, spaceDim);
      isSurfaceMesh = true;
    }
    else if (spaceDim == 3)
    {
      build = mesh.initialize(spaceDim, spaceDim);
      isSurfaceMesh = false;
    }
    else if (spaceDim == 2)
    {
      build = mesh.initialize(spaceDim, spaceDim);
      isSurfaceMesh = false;
    }
    else
    {
      Alert::Exception("Unhandled case.").raise();
    }

    is.clear();
    for (const auto& [_, kw] : Medit::KeywordMap)
    {
      if (!m_context.pos.count(kw))
        continue;
      auto g = m_context.pos.at(kw);
      if (!g.has_value())
        continue;
      is.seekg(*g);

      if (!m_context.count.count(kw))
      {
        Alert::Exception()
          << "Failed to parse \"" << kw << "\" count." << Alert::Raise;
      }

      switch (kw)
      {
        case Medit::Keyword::Vertices:
        {
          // Read all vertices
          if (spaceDim == 2)
          {
            for (size_t i = 0; i < m_context.count.at(kw); i++)
            {
              if(!std::getline(is, line))
                Alert::Exception("Bad mesh format.").raise();
              std::istringstream lss(line);
              double x, y;
              lss >> x >> y;
              // We ignore the reference
              build.vertex({x, y});
            }
          }
          else if (spaceDim == 3)
          {
            for (size_t i = 0; i < m_context.count.at(kw); i++)
            {
              if(!std::getline(is, line))
                Alert::Exception("Bad mesh format.").raise();
              std::istringstream lss(line);
              double x, y, z;
              lss >> x >> y >> z;
              // We ignore the reference
              build.vertex({x, y, z});
            }
          }
          break;
        }
        case Medit::Keyword::Edges:
        {
          assert(spaceDim >= 2);
          // Read all edges
          for (size_t i = 0; i < m_context.count.at(kw); i++)
          {
            if(!std::getline(is, line))
              Alert::Exception("Bad mesh format.").raise();
            std::istringstream lss(line);
            size_t v1, v2, ref;
            lss >> v1 >> v2 >> ref;
            switch (spaceDim)
            {
              case 2: // Planar mesh
              {
                build.face(Geometry::Type::Segment, {v1 - 1, v2 - 1});
                break;
              }
              case 3: // Surface mesh
              {
                if (isSurfaceMesh)
                {
                  build.face(Geometry::Type::Segment, {v1 - 1, v2 - 1});
                }
                break;
              }
              default:
              {
                Alert::Exception("Unhandled case.").raise();
                break;
              }
            }
          }
          break;
        }
        case Medit::Keyword::Triangles:
        {
          assert(spaceDim >= 2);
          // Read all triangles
          for (size_t i = 0; i < m_context.count.at(kw); i++)
          {
            if(!std::getline(is, line))
              Alert::Exception("Bad mesh format.").raise();
            std::istringstream lss(line);
            size_t v1, v2, v3, ref;
            lss >> v1 >> v2 >> v3 >> ref;
            switch (spaceDim)
            {
              case 2:
              {
                build.element(Geometry::Type::Triangle, {v1 - 1, v2 - 1, v3 - 1});
                break;
              }
              case 3:
              {
                if (isSurfaceMesh)
                {
                  build.element(Geometry::Type::Triangle, {v1 - 1, v2 - 1, v3 - 1});
                }
                else
                {
                  build.face(Geometry::Type::Triangle, {v1 - 1, v2 - 1, v3 - 1});
                }
                break;
              }
              default:
              {
                Alert::Exception("Unhandled case.").raise();
                break;
              }
            }
          }
          break;
        }
        case Medit::Keyword::Tetrahedra:
        {
          assert(spaceDim >= 3);
          // Read all tetrahedra
          for (size_t i = 0; i < m_context.count.at(kw); i++)
          {
            if(!std::getline(is, line))
              Alert::Exception("Bad mesh format.").raise();
            std::istringstream lss(line);
            size_t v1, v2, v3, v4, ref;
            lss >> v1 >> v2 >> v3 >> v4 >> ref;
            build.element(Geometry::Type::Tetrahedron, {v1 - 1, v2 - 1, v3 - 1, v4 - 1});
          }
          break;
        }
        default:
        {
          break;
        }
      }
    }
    build.finalize();
    getObject() = std::move(mesh);
  }
}
