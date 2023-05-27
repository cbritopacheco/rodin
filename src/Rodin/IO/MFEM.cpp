/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MFEM.h"

namespace Rodin::IO::MFEM
{}

namespace Rodin::IO
{
  void MeshLoader<FileFormat::MFEM, Context::Serial>::readHeader(std::istream& is)
  {
    auto line = skipEmptyLinesAndComments(is);
    auto header = MFEM::ParseHeader()(line.begin(), line.end());
    if (!header)
      Alert::Exception("Failed to determine MFEM mesh type and version.").raise();
    m_header = *header;
  }

  void MeshLoader<FileFormat::MFEM, Context::Serial>::readDimension(std::istream& is)
  {
    auto line = skipEmptyLinesAndComments(is);
    auto kw = MFEM::ParseKeyword()(line.begin(), line.end());
    if (!kw && *kw != MFEM::Keyword::dimension)
    {
      Alert::Exception() << "Expected "
                         << std::quoted(toCharString(MFEM::Keyword::dimension))
                         << " on line " << m_currentLineNumber << "."
                         << Alert::Raise;
    }
    line = skipEmptyLinesAndComments(is);
    auto dimension = MFEM::ParseUnsignedInteger()(line.begin(), line.end());
    if (!dimension)
    {
      Alert::Exception() << "Failed to determine mesh dimension on line "
                         << m_currentLineNumber << "."
                         << Alert::Raise;
    }
    m_dimension = *dimension;
  }

  void MeshLoader<FileFormat::MFEM, Context::Serial>::readMesh(std::istream& is)
  {
    Geometry::MeshConnectivity connectivity;
    connectivity.initialize(m_dimension);

    Geometry::PolytopeIndexed<Geometry::Attribute> attrs;
    attrs.initialize(m_dimension);

    std::string line;
    while (getline(is, line))
    {
      if (MFEM::ParseEmptyLineOrComment()(line.begin(), line.end()))
        continue;
      assert(line.size() > 0);
      auto kw = MFEM::ParseKeyword()(line.begin(), line.end());
      if (!kw)
      {
        Alert::Exception() << "Expected keyword on line "
                           << m_currentLineNumber
                           << "." << Alert::Raise;
      }
      auto keyword = MFEM::toKeyword(kw->c_str());
      if (!keyword)
      {
        Alert::Exception() << "Unrecognized keyword " << std::quoted(*kw)
                           << " on line " << m_currentLineNumber << "."
                           << Alert::Raise;
      }

      switch (*keyword)
      {
        case MFEM::Keyword::boundary:
        {
          getline(is, line);
          auto count = MFEM::ParseUnsignedInteger()(line.begin(), line.end());
          if (!count)
          {
            Alert::Exception() << "Failed to determine " << std::quoted(*kw)
                               << " count on line " << m_currentLineNumber << "."
                               << Alert::Raise;
          }
          attrs.reserve(m_dimension - 1, *count);
          for (size_t i = 0; i < *count; i++)
          {
            getline(is, line);
            auto g = MFEM::ParseGeometry()(line.begin(), line.end());
            if (!g)
            {
              Alert::Exception() << "Failed to parse geometry on line "
                                 << m_currentLineNumber << "."
                                 << Alert::Raise;
            }
            connectivity.polytope(g->geometry, std::move(g->vertices));
            attrs.track(m_dimension - 1, i, g->attribute);
          }
          continue;
        }
        case MFEM::Keyword::elements:
        {
          getline(is, line);
          auto count = MFEM::ParseUnsignedInteger()(line.begin(), line.end());
          if (!count)
          {
            Alert::Exception() << "Failed to determine " << std::quoted(*kw)
                               << " count on line " << m_currentLineNumber << "."
                               << Alert::Raise;
          }
          attrs.reserve(m_dimension, *count);
          for (size_t i = 0; i < *count; i++)
          {
            getline(is, line);
            auto g = MFEM::ParseGeometry()(line.begin(), line.end());
            if (!g)
            {
              Alert::Exception() << "Failed to parse geometry on line "
                                 << m_currentLineNumber << "."
                                 << Alert::Raise;
            }
            connectivity.polytope(g->geometry, std::move(g->vertices));
            attrs.track(m_dimension, i, g->attribute);
          }
          continue;
        }
        case MFEM::Keyword::vertices:
        {
          getline(is, line);
          auto count = MFEM::ParseUnsignedInteger()(line.begin(), line.end());
          if (!count)
          {
            Alert::Exception() << "Failed to determine " << std::quoted(*kw)
                               << " count on line " << m_currentLineNumber << "."
                               << Alert::Raise;
          }
          connectivity.nodes(*count);
          getline(is, line);
          auto vdim = MFEM::ParseUnsignedInteger()(line.begin(), line.end());
          if (vdim)
          {
            m_spaceDimension = *vdim;
            m_build.initialize(m_spaceDimension).nodes(*count);
            for (size_t i = 0; i < *count; i++)
            {
              getline(is, line);
              auto x = MFEM::ParseVertex(m_spaceDimension)(line.begin(), line.end());
              if (!x)
              {
                Alert::Exception() << "Failed to parse vertex on line " << m_currentLineNumber
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
              Alert::Exception() << "Failed to parse vertex on line " << m_currentLineNumber
                                 << Alert::Raise;
            }
            m_build.vertex(std::move(*x));
            for (size_t i = 1; i < *count; i++)
            {
              auto x = MFEM::ParseVertex(m_spaceDimension)(line.begin(), line.end());
              if (!x)
              {
                Alert::Exception() << "Failed to parse vertex on line " << m_currentLineNumber
                                   << Alert::Raise;
              }
              m_build.vertex(std::move(*x));
            }
          }
          continue;
        }
        default:
        {
          Alert::Exception() << "Unexpected keyword " << std::quoted(*kw)
                             << " on line " << m_currentLineNumber << "."
                             << Alert::Raise;
        }
      }
    }
    m_build.setConnectivity(std::move(connectivity));
    m_build.setAttributeIndex(std::move(attrs));
  }

  std::string
  MeshLoader<FileFormat::MFEM, Context::Serial>::skipEmptyLinesAndComments(std::istream& is)
  {
    std::string line;
    while (getline(is, line))
    {
      if (!MFEM::ParseEmptyLineOrComment()(line.begin(), line.end()))
        break;
    }
    return line;
  }

  std::istream&
  MeshLoader<FileFormat::MFEM, Context::Serial>::getline(std::istream& is, std::string& line)
  {
    m_currentLineNumber++;
    return std::getline(is, line);
  }

  void MeshLoader<IO::FileFormat::MFEM, Context::Serial>::load(std::istream &is)
  {
    readHeader(is);
    readDimension(is);
    readMesh(is);
    getObject() = m_build.finalize();
  }

  void MeshPrinter<FileFormat::MFEM, Context::Serial>::print(std::ostream &os)
  {
    printHeader(os);
    printDimension(os);
    printMesh(os);
  }

  void MeshPrinter<FileFormat::MFEM, Context::Serial>::printHeader(std::ostream &os)
  {
    os << "MFEM mesh v1.0\n\n";
  }

  void MeshPrinter<FileFormat::MFEM, Context::Serial>::printDimension(std::ostream &os)
  {
    const auto& mesh = getObject();
    os << "dimension\n" << mesh.getDimension() << "\n\n";
  }

  void MeshPrinter<FileFormat::MFEM, Context::Serial>::printMesh(std::ostream &os)
  {
    const auto& mesh = getObject();
    os << MFEM::Keyword::elements << '\n' << mesh.getElementCount() << '\n';
    for (auto it = mesh.getElement(); !it.end(); ++it)
    {
      auto g = MFEM::getGeometry(it->getGeometry());
      if (!g)
      {
        Alert::Exception() << "MFEM format does not support geometry: "
                           << it->getGeometry() << "."
                           << Alert::Raise;
      }
      os << it->getAttribute() << ' ' << *g << ' ';

      const auto& vertices = it->getVertices();
      for (int i = 0; i < vertices.size() - 1; i++)
        os << vertices(i) << ' ';
      os << vertices(vertices.size() - 1) << '\n';
    }
    os << '\n';

    const size_t fdim = mesh.getDimension() - 1;
    os << MFEM::Keyword::boundary << '\n' << mesh.getAttributeIndex().size(fdim) << '\n';
    const auto& attrs = mesh.getAttributeIndex();
    for (auto it = attrs.begin(mesh.getDimension() - 1); it != attrs.end(mesh.getDimension() - 1); ++it)
    {
      const auto& [i, attr] = *it;
      auto g = MFEM::getGeometry(mesh.getConnectivity().getGeometry(fdim, i));
      if (!g)
      {
        Alert::Exception() << "MFEM format does not support geometry: "
                           << mesh.getConnectivity().getGeometry(fdim, i) << "."
                           << Alert::Raise;
      }
      os << attr << ' ' << *g << ' ';

      const auto& vertices = mesh.getConnectivity().getPolytope(fdim, i);
      for (int i = 0; i < vertices.size() - 1; i++)
        os << vertices(i) << ' ';
      os << vertices(vertices.size() - 1) << '\n';
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


