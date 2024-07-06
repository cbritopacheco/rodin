/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <boost/algorithm/string.hpp>

#include "MEDIT.h"

namespace Rodin::IO
{
  std::istream& MeshLoader<FileFormat::MEDIT, Context::Sequential>::getline(std::istream& is, std::string& line)
  {
    m_currentLineNumber++;
    return std::getline(is, line);
  }

  std::string MeshLoader<FileFormat::MEDIT, Context::Sequential>::skipEmptyLines(std::istream& is)
  {
    std::string line;
    while (getline(is, line))
    {
      if (!MEDIT::ParseEmptyLine()(line.begin(), line.end()))
        break;
    }
    return line;
  }

  void MeshLoader<FileFormat::MEDIT, Context::Sequential>::readVersion(std::istream& is)
  {
    auto line = skipEmptyLines(is);
    std::optional<unsigned int> version = MEDIT::ParseMeshVersionFormatted()(line.begin(), line.end());
    if (version) // Version was on the same line
    {
      m_version = *version;
    }
    else // Version is not on the same line
    {
      auto line = skipEmptyLines(is);
      version = MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
      if (version)
        m_version = *version;
      else
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Failed to parse version number of mesh." << Alert::Raise;
      }
    }
  }

  void MeshLoader<FileFormat::MEDIT, Context::Sequential>::readDimension(std::istream& is)
  {
    auto line = skipEmptyLines(is);
    std::optional<unsigned int> dimension = MEDIT::ParseDimension()(line.begin(), line.end());
    if (dimension) // Version was on the same line
      m_spaceDimension = *dimension;
    else // Version is not on the same line
    {
      auto line = skipEmptyLines(is);
      dimension = MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
      if (dimension)
        m_spaceDimension = *dimension;
      else
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Failed to parse dimension of mesh." << Alert::Raise;
      }
    }
  }

  void MeshLoader<FileFormat::MEDIT, Context::Sequential>::readEntities(std::istream& is)
  {
    std::string line;
    while (getline(is, line))
    {
      if (MEDIT::ParseEmptyLine()(line.begin(), line.end()))
        continue;
      assert(line.size() > 0);
      if (!(std::isalpha(line[0]) || std::isspace(line[0])))
        continue;
      auto kw = MEDIT::ParseKeyword()(line.begin(), line.end());
      if (!kw)
        continue;
      auto entity = MEDIT::toKeyword(kw->c_str());
      if (!entity)
      {
        Alert::Warning() << "Ignoring unrecognized keyword: " << *kw << Alert::Raise;
        continue;
      }
      if (*entity == MEDIT::Keyword::End)
        break;

      line = skipEmptyLines(is);
      auto count = MEDIT::ParseUnsignedInteger()(line.begin(), line.end());

      if (!count)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Failed to determine number of "
          << std::quoted(*kw) << "." << Alert::Raise;
      }

      m_pos[*entity] = is.tellg();
      m_count[*entity] = *count;

      switch (*entity)
      {
        case MEDIT::Keyword::Vertices:
        {
          m_build.nodes(*count);
          for (size_t i = 0; i < *count; i++)
          {
            getline(is, line);
            auto data = MEDIT::ParseVertex(m_spaceDimension)(line.begin(), line.end());
            if (!data)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to parse Vertex on line "
                << std::to_string(m_currentLineNumber)
                << "."
                << Alert::Raise;
            }
            m_build.vertex(std::move(data->vertex));
            if (data->attribute != RODIN_DEFAULT_POLYTOPE_ATTRIBUTE)
              m_build.attribute({ 0, i }, data->attribute);
          }
          continue; // Continue the while loop
        }
        case MEDIT::Keyword::Edges:
        {
          m_build.reserve(1, *count);
          for (size_t i = 0; i < *count; i++)
          {
            getline(is, line);
            auto data = MEDIT::ParseEntity(2)(line.begin(), line.end());
            if (!data)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to parse Edge on line "
                << std::to_string(m_currentLineNumber)
                << "."
                << Alert::Raise;
            }
            data->vertices -= 1;
            m_build.polytope(Geometry::Polytope::Type::Segment, std::move(data->vertices));
            if (data->attribute != RODIN_DEFAULT_POLYTOPE_ATTRIBUTE)
              m_build.attribute({ 1, i }, data->attribute);
          }
          continue; // Continue the while loop
        }
        case MEDIT::Keyword::Triangles:
        {
          m_build.reserve(2, *count);
          for (size_t i = 0; i < *count; i++)
          {
            getline(is, line);
            auto data = MEDIT::ParseEntity(3)(line.begin(), line.end());
            if (!data)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to parse Triangle on line "
                << std::to_string(m_currentLineNumber)
                << "."
                << Alert::Raise;
            }
            data->vertices -= 1;
            m_build.polytope(Geometry::Polytope::Type::Triangle, std::move(data->vertices));
            if (data->attribute != RODIN_DEFAULT_POLYTOPE_ATTRIBUTE)
              m_build.attribute({ 2, i }, data->attribute);
          }
          continue; // Continue the while loop
        }
        case MEDIT::Keyword::Quadrilaterals:
        {
          m_build.reserve(2, *count);
          for (size_t i = 0; i < *count; i++)
          {
            getline(is, line);
            auto data = MEDIT::ParseEntity(4)(line.begin(), line.end());
            if (!data)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to parse Quadrilateral on line "
                << std::to_string(m_currentLineNumber)
                << "."
                << Alert::Raise;
            }
            data->vertices -= 1;
            std::swap(data->vertices(2), data->vertices(3));
            m_build.polytope(Geometry::Polytope::Type::Quadrilateral, std::move(data->vertices));
            if (data->attribute != RODIN_DEFAULT_POLYTOPE_ATTRIBUTE)
              m_build.attribute({ 2, i }, data->attribute);
          }
          continue; // Continue the while loop
        }
        case MEDIT::Keyword::Tetrahedra:
        {
          m_build.reserve(3, *count);
          for (size_t i = 0; i < *count; i++)
          {
            getline(is, line);
            auto data = MEDIT::ParseEntity(4)(line.begin(), line.end());
            if (!data)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to parse Tetrahedron on line "
                << std::to_string(m_currentLineNumber)
                << "."
                << Alert::Raise;
            }
            data->vertices -= 1;
            m_build.polytope(Geometry::Polytope::Type::Tetrahedron, std::move(data->vertices));
            if (data->attribute != RODIN_DEFAULT_POLYTOPE_ATTRIBUTE)
              m_build.attribute({ 3, i }, data->attribute);
          }
          continue; // Continue the while loop
        }
        default:
          continue; // Continue the while loop
      }
      break; // If we reach here we should break out of the while loop
    }
  }

  void MeshLoader<FileFormat::MEDIT, Context::Sequential>::load(std::istream& is)
  {
    readVersion(is);
    readDimension(is);
    m_build.initialize(m_spaceDimension);
    readEntities(is);
    getObject() = m_build.finalize();
  }

  void MeshPrinter<FileFormat::MEDIT, Context::Sequential>::printDimension(std::ostream& os)
  {
    const auto& mesh = getObject();
    os << MEDIT::Keyword::Dimension << '\n' << mesh.getSpaceDimension() << "\n\n";
  }

  void MeshPrinter<FileFormat::MEDIT, Context::Sequential>::printVersion(std::ostream& os)
  {
    os << MEDIT::Keyword::MeshVersionFormatted << "\n2" << "\n\n";
  }

  void MeshPrinter<FileFormat::MEDIT, Context::Sequential>::printEntities(std::ostream& os)
  {
    const auto& mesh = getObject();

    for (auto g : Geometry::Polytope::Types)
    {
      switch (g)
      {
        case Geometry::Polytope::Type::Point:
        {
          os << MEDIT::Keyword::Vertices << '\n' << mesh.getVertexCount() << '\n';
          for (auto it = mesh.getVertex(); !it.end(); ++it)
          {
            for (const auto& x : it->getCoordinates())
              os << x << " ";
            os << it->getAttribute() << '\n';
          }

          os << '\n';
          break;
        }
        default:
        {
          switch (g)
          {
            case Geometry::Polytope::Type::Point:
            {
              assert(false);
              break;
            }
            case Geometry::Polytope::Type::Segment:
            {
              os << MEDIT::Keyword::Edges << '\n';
              break;
            }
            case Geometry::Polytope::Type::Triangle:
            {
              os << MEDIT::Keyword::Triangles << '\n';
              break;
            }
            case Geometry::Polytope::Type::Quadrilateral:
            {
              os << MEDIT::Keyword::Quadrilaterals << '\n';
              break;
            }
            case Geometry::Polytope::Type::Tetrahedron:
            {
              os << MEDIT::Keyword::Tetrahedra << '\n';
              break;
            }
          }
          const size_t d = Geometry::Polytope::getGeometryDimension(g);
          if (d <= mesh.getDimension())
          {
            os << mesh.getPolytopeCount(g) << '\n';
            for (auto it = mesh.getPolytope(d); !it.end(); ++it)
            {
              if (it->getGeometry() == g)
              {
                const auto& vertices = it->getVertices();
                switch (it->getGeometry())
                {
                  case Geometry::Polytope::Type::Point:
                  {
                    os << vertices(0) + 1;
                    break;
                  }
                  case Geometry::Polytope::Type::Triangle:
                  {
                    os << vertices(0) + 1 << ' ' << vertices(1) + 1 << ' ' << vertices(2) + 1;
                    break;
                  }
                  case Geometry::Polytope::Type::Segment:
                  {
                    os << vertices(0) + 1 << ' ' << vertices(1) + 1;
                    break;
                  }
                  case Geometry::Polytope::Type::Tetrahedron:
                  {
                    os << vertices(0) + 1 << ' ' << vertices(1) + 1 << ' '
                       << vertices(2) + 1 << ' ' << vertices(3) + 1;
                    break;
                  }
                  case Geometry::Polytope::Type::Quadrilateral:
                  {
                    os << vertices(0) + 1 << ' ' << vertices(1) + 1 << ' '
                       << vertices(3) + 1 << ' ' << vertices(2) + 1;
                    break;
                  }
                }
                os << ' ' << it->getAttribute() << '\n';
              }
            }
          }
          else
          {
            os << 0 << '\n';
          }
          os << '\n';
          break;
        }
      }
    }
  }

  void MeshPrinter<FileFormat::MEDIT, Context::Sequential>::printEnd(std::ostream& os)
  {
    os << '\n' << IO::MEDIT::Keyword::End;
  }

  void MeshPrinter<FileFormat::MEDIT, Context::Sequential>::print(std::ostream& os, bool pEnd)
  {
    printVersion(os);
    printDimension(os);
    printEntities(os);
    if (pEnd) printEnd(os);
  }
}
