/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MFEM.h"

namespace Rodin::IO::MFEM
{
  std::istream& getline(std::istream& is, std::string& line, size_t& currentLineNumber)
  {
    currentLineNumber++;
    return std::getline(is, line);
  }

  std::string skipEmptyLinesAndComments(std::istream& is, size_t& currentLineNumber)
  {
    std::string line;
    while (getline(is, line))
    {
      if (!MFEM::ParseEmptyLineOrComment()(line.begin(), line.end()))
        break;
      currentLineNumber++;
    }
    return line;
  }
}

namespace Rodin::IO
{
  void MeshLoader<FileFormat::MFEM, Context::Sequential>::readHeader(std::istream& is)
  {
    auto line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
    auto header = MFEM::ParseMeshHeader()(line.begin(), line.end());
    if (!header)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to determine MFEM mesh type and version." << Alert::Raise;
    }
    m_header = *header;
  }

  void MeshLoader<FileFormat::MFEM, Context::Sequential>::readDimension(std::istream& is)
  {
    auto line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
    auto kw = MFEM::ParseKeyword()(line.begin(), line.end());
    if (!kw && *kw != MFEM::Keyword::dimension)
    {
      Alert::MemberFunctionException(*this, __func__)
         << "Expected "
         << std::quoted(toCharString(MFEM::Keyword::dimension))
         << " on line " << m_currentLineNumber << "."
         << Alert::Raise;
    }
    line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
    auto dimension = MFEM::ParseUnsignedInteger()(line.begin(), line.end());
    if (!dimension)
    {
      Alert::MemberFunctionException(*this, __func__) << "Failed to determine mesh dimension on line "
                         << m_currentLineNumber << "."
                         << Alert::Raise;
    }
    m_dimension = *dimension;
  }

  void MeshLoader<FileFormat::MFEM, Context::Sequential>::readMesh(std::istream& is)
  {
    Geometry::Connectivity<Context::Sequential> connectivity;
    connectivity.initialize(m_dimension);

    Geometry::PolytopeIndexed<Geometry::Attribute> attrs;
    attrs.initialize(m_dimension);

    std::string line;
    while (MFEM::getline(is, line, m_currentLineNumber))
    {
      if (MFEM::ParseEmptyLineOrComment()(line.begin(), line.end()))
        continue;
      assert(line.size() > 0);
      auto kw = MFEM::ParseKeyword()(line.begin(), line.end());
      if (!kw)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Expected keyword on line "
          << m_currentLineNumber
          << "." << line << Alert::Raise;
      }
      auto keyword = MFEM::toKeyword(kw->c_str());
      if (!keyword)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Unrecognized keyword " << std::quoted(*kw)
          << " on line " << m_currentLineNumber << "."
          << Alert::Raise;
      }

      switch (*keyword)
      {
        case MFEM::Keyword::boundary:
        {
          MFEM::getline(is, line, m_currentLineNumber);
          auto count = MFEM::ParseUnsignedInteger()(line.begin(), line.end());
          if (!count)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to determine " << std::quoted(*kw)
              << " count on line " << m_currentLineNumber << "."
              << Alert::Raise;
          }
          attrs.reserve(m_dimension - 1, *count);
          for (size_t i = 0; i < *count; i++)
          {
            MFEM::getline(is, line, m_currentLineNumber);
            auto g = MFEM::ParseGeometry()(line.begin(), line.end());
            if (!g)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to parse geometry on line "
                << m_currentLineNumber << "."
                << Alert::Raise;
            }
            if (g->geometry == Geometry::Polytope::Type::Quadrilateral)
              std::swap(g->vertices(2), g->vertices(3));
            connectivity.polytope(g->geometry, std::move(g->vertices));
            attrs.track({ m_dimension - 1, i }, g->attribute);
          }
          continue;
        }
        case MFEM::Keyword::elements:
        {
          MFEM::getline(is, line, m_currentLineNumber);
          auto count = MFEM::ParseUnsignedInteger()(line.begin(), line.end());
          if (!count)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to determine " << std::quoted(*kw)
              << " count on line " << m_currentLineNumber << "."
              << Alert::Raise;
          }
          attrs.reserve(m_dimension, *count);
          for (size_t i = 0; i < *count; i++)
          {
            MFEM::getline(is, line, m_currentLineNumber);
            auto g = MFEM::ParseGeometry()(line.begin(), line.end());
            if (!g)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to parse geometry on line "
                << m_currentLineNumber << "."
                << Alert::Raise;
            }
            if (g->geometry == Geometry::Polytope::Type::Quadrilateral)
              std::swap(g->vertices(2), g->vertices(3));
            connectivity.polytope(g->geometry, std::move(g->vertices));
            attrs.track({ m_dimension, i }, g->attribute);
          }
          continue;
        }
        case MFEM::Keyword::vertices:
        {
          MFEM::getline(is, line, m_currentLineNumber);
          auto count = MFEM::ParseUnsignedInteger()(line.begin(), line.end());
          if (!count)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to determine " << std::quoted(*kw)
              << " count on line " << m_currentLineNumber << "."
              << Alert::Raise;
          }
          connectivity.nodes(*count);
          MFEM::getline(is, line, m_currentLineNumber);
          auto vdim = MFEM::ParseUnsignedInteger()(line.begin(), line.end());
          if (vdim)
          {
            m_spaceDimension = *vdim;
            m_build.initialize(m_spaceDimension).nodes(*count);
            for (size_t i = 0; i < *count; i++)
            {
              MFEM::getline(is, line, m_currentLineNumber);
              auto x = MFEM::ParseVertex(m_spaceDimension)(line.begin(), line.end());
              if (!x)
              {
                Alert::MemberFunctionException(*this, __func__)
                  << "Failed to parse vertex on line " << m_currentLineNumber
                  << Alert::Raise;
              }
              m_build.vertex(std::move(*x));
            }
          }
          else
          {
            m_spaceDimension = m_dimension;
            m_build.initialize(m_spaceDimension).nodes(m_spaceDimension);
            auto x = MFEM::ParseVertex(m_spaceDimension)(line.begin(), line.end());
            if (!x)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to parse vertex on line " << m_currentLineNumber
                << Alert::Raise;
            }
            m_build.vertex(std::move(*x));
            for (size_t i = 1; i < *count; i++)
            {
              auto x = MFEM::ParseVertex(m_spaceDimension)(line.begin(), line.end());
              if (!x)
              {
                Alert::MemberFunctionException(*this, __func__)
                  << "Failed to parse vertex on line " << m_currentLineNumber
                  << Alert::Raise;
              }
              m_build.vertex(std::move(*x));
            }
          }
          continue;
        }
        default:
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Unexpected keyword " << std::quoted(*kw)
            << " on line " << m_currentLineNumber << "."
            << Alert::Raise;
        }
      }
    }
    m_build.setConnectivity(std::move(connectivity));
    m_build.setAttributeIndex(std::move(attrs));
  }

  void MeshLoader<IO::FileFormat::MFEM, Context::Sequential>::load(std::istream &is)
  {
    readHeader(is);
    readDimension(is);
    readMesh(is);
    getObject() = m_build.finalize();
  }

  void MeshPrinter<FileFormat::MFEM, Context::Sequential>::print(std::ostream &os)
  {
    printHeader(os);
    printDimension(os);
    printMesh(os);
  }

  void MeshPrinter<FileFormat::MFEM, Context::Sequential>::printHeader(std::ostream &os)
  {
    os << "MFEM mesh v1.0\n\n";
  }

  void MeshPrinter<FileFormat::MFEM, Context::Sequential>::printDimension(std::ostream &os)
  {
    const auto& mesh = getObject();
    os << "dimension\n" << mesh.getDimension() << "\n\n";
  }

  void MeshPrinter<FileFormat::MFEM, Context::Sequential>::printMesh(std::ostream &os)
  {
    const auto& mesh = getObject();
    os << MFEM::Keyword::elements << '\n' << mesh.getCellCount() << '\n';
    for (auto it = mesh.getCell(); !it.end(); ++it)
    {
      auto g = MFEM::getGeometry(it->getGeometry());
      if (!g)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "MFEM format does not support geometry: "
          << it->getGeometry() << "."
          << Alert::Raise;
      }
      os << it->getAttribute() << ' ' << *g << ' ';

      const auto& vertices = it->getVertices();
      switch (it->getGeometry())
      {
        case Geometry::Polytope::Type::Point:
        {
          os << vertices(0);
          break;
        }
        case Geometry::Polytope::Type::Triangle:
        {
          os << vertices(0) << ' ' << vertices(1) << ' ' << vertices(2);
          break;
        }
        case Geometry::Polytope::Type::Segment:
        {
          os << vertices(0) << ' ' << vertices(1);
          break;
        }
        case Geometry::Polytope::Type::Tetrahedron:
        {
          os << vertices(0) << ' ' << vertices(1) << ' ' << vertices(2) << ' ' << vertices(3);
          break;
        }
        case Geometry::Polytope::Type::Quadrilateral:
        {
          os << vertices(0) << ' ' << vertices(1) << ' ' << vertices(3) << ' ' << vertices(2);
          break;
        }
      }
      os << '\n';
    }
    os << '\n';

    os << MFEM::Keyword::boundary << '\n' << mesh.getFaceCount() << '\n';
    for (auto it = mesh.getFace(); !it.end(); ++it)
    {
      auto g = MFEM::getGeometry(it->getGeometry());
      if (!g)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "MFEM format does not support geometry: "
          << it->getGeometry() << "."
          << Alert::Raise;
      }
      os << it->getAttribute() << ' ' << *g << ' ';

      const auto& vertices = it->getVertices();
      switch (it->getGeometry())
      {
        case Geometry::Polytope::Type::Point:
        {
          os << vertices(0);
          break;
        }
        case Geometry::Polytope::Type::Triangle:
        {
          os << vertices(0) << ' ' << vertices(1) << ' ' << vertices(2);
          break;
        }
        case Geometry::Polytope::Type::Segment:
        {
          os << vertices(0) << ' ' << vertices(1);
          break;
        }
        case Geometry::Polytope::Type::Tetrahedron:
        {
          os << vertices(0) << ' ' << vertices(1) << ' ' << vertices(2) << ' ' << vertices(3);
          break;
        }
        case Geometry::Polytope::Type::Quadrilateral:
        {
          os << vertices(0) << ' ' << vertices(1) << ' ' << vertices(3) << ' ' << vertices(2);
          break;
        }
      }
      os << '\n';
    }
    os << '\n';

    os << MFEM::Keyword::vertices << '\n'
       << mesh.getVertexCount() << '\n'
       << mesh.getSpaceDimension() << '\n';

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto& x = it->getCoordinates();
      for (int i = 0; i < x.size() - 1; i++)
        os << x(i) << ' ';
      os << x(x.size() - 1) << '\n';
    }
  }
}


