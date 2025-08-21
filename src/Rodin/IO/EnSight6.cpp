#include "EnSight6.h"

#include "Rodin/Configure.h"

namespace Rodin::IO
{
  MeshLoader<FileFormat::ENSIGHT6, Context::Local>::MeshLoader(ObjectType& mesh)
    : MeshLoaderBase(mesh),
      m_build(mesh),
      m_currentLineNumber(0),
      m_spaceDimension(3)
  {}

  std::istream& MeshLoader<FileFormat::ENSIGHT6, Context::Local>::getline(std::istream& is, std::string& line)
  {
    m_currentLineNumber++;
    return std::getline(is, line);
  }

  std::string MeshLoader<FileFormat::ENSIGHT6, Context::Local>::skipEmptyLines(std::istream& is)
  {
    std::string line;
    while (getline(is, line))
    {
      if (!EnSight6::ParseEmptyLine()(line.begin(), line.end()))
        break;
    }
    return line;
  }

  void MeshLoader<FileFormat::ENSIGHT6, Context::Local>::load(std::istream& is)
  {
    readHeader(is);
    readCoordinates(is);
    readParts(is);
    getObject().finalize();
  }

  void MeshLoader<FileFormat::ENSIGHT6, Context::Local>::readHeader(std::istream& is)
  {
    // Read the first description line
    auto line = skipEmptyLines(is);
    m_descriptionLine1 = line;

    // Read the second description line
    line = skipEmptyLines(is);
    m_descriptionLine2 = line;

    // Read node id configuration
    line = skipEmptyLines(is);
    // Expected format: "node id given" or "node id off"

    // Read element id configuration  
    line = skipEmptyLines(is);
    // Expected format: "element id given" or "element id off"
  }

  void MeshLoader<FileFormat::ENSIGHT6, Context::Local>::readCoordinates(std::istream& is)
  {
    auto line = skipEmptyLines(is);
    auto keyword = EnSight6::ParseKeyword()(line.begin(), line.end());
    if (!keyword || *keyword != "coordinates")
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Expected 'coordinates' keyword but found: " << (keyword ? *keyword : "invalid")
        << Alert::Raise;
    }

    // Read number of vertices
    line = skipEmptyLines(is);
    auto vertexCount = EnSight6::ParseUnsignedInteger()(line.begin(), line.end());
    if (!vertexCount)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to parse vertex count at line " << m_currentLineNumber
        << Alert::Raise;
    }

    // Reserve space for vertices
    m_build.reserve(0, *vertexCount);

    // Read vertex coordinates
    for (size_t i = 0; i < *vertexCount; i++)
    {
      line = skipEmptyLines(is);
      std::istringstream iss(line);
      
      Index vertexId;
      Real x, y, z;
      iss >> vertexId >> x >> y >> z;
      
      if (iss.fail())
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Failed to parse vertex coordinates at line " << m_currentLineNumber
          << Alert::Raise;
      }

      Math::Vector<Real> coords(m_spaceDimension);
      coords(0) = x;
      if (m_spaceDimension > 1) coords(1) = y;
      if (m_spaceDimension > 2) coords(2) = z;

      m_build.vertex(coords);
    }
  }

  void MeshLoader<FileFormat::ENSIGHT6, Context::Local>::readParts(std::istream& is)
  {
    std::string line;
    while (getline(is, line))
    {
      if (EnSight6::ParseEmptyLine()(line.begin(), line.end()))
        continue;

      auto keyword = EnSight6::ParseKeyword()(line.begin(), line.end());
      if (!keyword)
        continue;

      if (*keyword == "part")
      {
        // Read part number
        std::istringstream iss(line);
        std::string partKeyword;
        Geometry::Attribute partId;
        iss >> partKeyword >> partId;

        // Read part description/name
        line = skipEmptyLines(is);
        std::string partName = line;

        // Read elements in this part
        while (getline(is, line))
        {
          if (EnSight6::ParseEmptyLine()(line.begin(), line.end()))
            continue;

          auto elementKeyword = EnSight6::ParseKeyword()(line.begin(), line.end());
          if (!elementKeyword)
            continue;

          auto elementType = EnSight6::toElementType(elementKeyword->c_str());
          if (!elementType)
          {
            // Might be another part or end of file
            is.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            m_currentLineNumber--;
            break;
          }

          auto polytopeType = EnSight6::toPolytopeType(*elementType);
          if (!polytopeType)
          {
            Alert::Warning() << "Unsupported element type: " << *elementKeyword << Alert::Raise;
            continue;
          }

          // Read element count
          line = skipEmptyLines(is);
          auto elementCount = EnSight6::ParseUnsignedInteger()(line.begin(), line.end());
          if (!elementCount)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to parse element count for " << *elementKeyword
              << " at line " << m_currentLineNumber
              << Alert::Raise;
          }

          // Determine dimension for reserve
          size_t dim = 0;
          switch (*polytopeType)
          {
            case Geometry::Polytope::Type::Point:
              dim = 0;
              break;
            case Geometry::Polytope::Type::Segment:
              dim = 1;
              break;
            case Geometry::Polytope::Type::Triangle:
            case Geometry::Polytope::Type::Quadrilateral:
              dim = 2;
              break;
            case Geometry::Polytope::Type::Tetrahedron:
            case Geometry::Polytope::Type::Wedge:
              dim = 3;
              break;
          }

          // Reserve space for elements
          m_build.reserve(dim, *elementCount);

          // Read elements
          for (size_t i = 0; i < *elementCount; i++)
          {
            line = skipEmptyLines(is);
            std::istringstream iss(line);
            
            Index elementId;
            iss >> elementId;

            Array<Index> vertices;
            int numVertices = 0;
            
            // Determine number of vertices based on element type
            switch (*elementType)
            {
              case EnSight6::ElementType::point:
                numVertices = 1;
                break;
              case EnSight6::ElementType::bar2:
                numVertices = 2;
                break;
              case EnSight6::ElementType::tria3:
                numVertices = 3;
                break;
              case EnSight6::ElementType::quad4:
                numVertices = 4;
                break;
              case EnSight6::ElementType::tetra4:
                numVertices = 4;
                break;
              case EnSight6::ElementType::penta6:
                numVertices = 6;
                break;
              case EnSight6::ElementType::hexa8:
                numVertices = 8;
                break;
              default:
                Alert::Warning() << "Unsupported element type for vertex count: " << *elementKeyword << Alert::Raise;
                continue;
            }

            vertices.resize(numVertices);
            for (int j = 0; j < numVertices; j++)
            {
              iss >> vertices[j];
              // Convert from 1-based to 0-based indexing if needed
              // EnSight6 typically uses 0-based indexing, but check the data format
            }

            if (iss.fail())
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to parse element vertices at line " << m_currentLineNumber
                << Alert::Raise;
            }

            m_build.polytope(*polytopeType, std::move(vertices));
            if (partId != 0) // 0 is typically the default attribute
              m_build.attribute({ dim, i }, partId);
          }
        }
      }
    }
  }

  MeshPrinter<FileFormat::ENSIGHT6, Context::Local>::MeshPrinter(const ObjectType& mesh)
    : MeshPrinterBase(mesh),
      m_descriptionLine1("EnSight6 Geometry File Format"),
      m_descriptionLine2("Rodin v" RODIN_VERSION)
  {}

  void MeshPrinter<FileFormat::ENSIGHT6, Context::Local>::print(std::ostream& os)
  {
    printHeader(os);
    printCoordinates(os);
    printParts(os);
  }

  void MeshPrinter<FileFormat::ENSIGHT6, Context::Local>::printHeader(std::ostream& os)
  {
    os << m_descriptionLine1 << '\n'
       << m_descriptionLine2 << '\n'
       << EnSight6::Keyword::node << ' ' << EnSight6::Keyword::id << ' ' << EnSight6::Keyword::given << '\n'
       << EnSight6::Keyword::element << ' ' << EnSight6::Keyword::id << ' ' << EnSight6::Keyword::given << '\n';
  }

  void MeshPrinter<FileFormat::ENSIGHT6, Context::Local>::printCoordinates(std::ostream& os)
  {
    os << EnSight6::Keyword::coordinates << '\n';
    const auto& mesh = getObject();
    os << mesh.getVertexCount() << '\n';
    os << std::setprecision(5) << std::scientific;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      // Print the vertex index
      os << it->getIndex() << ' ';

      // Retrieve dynamic coordinate vector and default to zero for missing dimensions
      const auto& coords = it->getCoordinates();
      Real x0 = 0.0, x1 = 0.0, x2 = 0.0;
      if (coords.size() > 0) x0 = coords(0);
      if (coords.size() > 1) x1 = coords(1);
      if (coords.size() > 2) x2 = coords(2);

      // Always write three components: X, Y, Z
      os << std::setw(12) << x0
         << std::setw(12) << x1
         << std::setw(12) << x2
         << '\n';
    }
  }

  void MeshPrinter<FileFormat::ENSIGHT6, Context::Local>::printParts(std::ostream& os)
  {
    os.unsetf(std::ios::scientific | std::ios::fixed); // removes scientific and fixed
    os << std::setprecision(6);                        // reset to default precision
    os << std::defaultfloat;                           // optional, returns to normal float format

    const auto& mesh = getObject();
    // The attributes are the part names.
    IndexSet attributes;
    for (size_t d = 0; d <= mesh.getDimension(); d++)
    {
      const auto& as = mesh.getAttributes(d);
      attributes.insert(as.begin(), as.end());
    }

    FlatMap<Geometry::Attribute, Geometry::GeometryIndexed<size_t>> count;
    FlatMap<Geometry::Attribute, Geometry::GeometryIndexed<std::ostringstream>> ess;
    for (size_t d = 0; d <= mesh.getDimension(); d++)
    {
      for (auto it = mesh.getPolytope(d); it; ++it)
      {
        const auto& geometry = it->getGeometry();
        const auto& vertices = it->getVertices();
        if (mesh.getPolytopeCount(geometry) == 0)
          continue;
        switch (geometry)
        {
          case Geometry::Polytope::Type::Point:
          {
            ess[it->getAttribute()][geometry]
              << it->getIndex()
              << ' ' << it->getIndex() << '\n';
            break;
          }
          case Geometry::Polytope::Type::Segment:
          {
            ess[it->getAttribute()][geometry]
              << it->getIndex() << ' ' << vertices(0) << ' ' << vertices(1) << '\n';
            break;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            ess[it->getAttribute()][geometry]
              << it->getIndex() << ' ' << vertices(0) << ' ' << vertices(1) << ' ' << vertices(2) << '\n';
            break;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            ess[it->getAttribute()][geometry]
              << it->getIndex() << ' ' << vertices(0) << ' ' << vertices(1) << ' '
              << vertices(3) << ' ' << vertices(2) << '\n';
            break;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            ess[it->getAttribute()][geometry]
              << it->getIndex() << ' ' << vertices(0) << ' ' << vertices(1) << ' '
              << vertices(2) << ' ' << vertices(3) << '\n';
            break;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            ess[it->getAttribute()][geometry]
              << it->getIndex() << ' ' << vertices(0) << ' ' << vertices(1) << ' '
              << vertices(2) << ' ' << vertices(3) << ' ' << vertices(4) << ' ' << vertices(5) << '\n';
            break;
          }
        }
        count[it->getAttribute()][geometry] += 1;
      }
    }

    for (const auto& attr : attributes)
    {
      os << EnSight6::Keyword::part << ' ' << attr << '\n' << "Attribute_" << attr << '\n';
      for (const auto& geometry : Geometry::Polytope::Types)
      {
        const size_t cnt = count[attr][geometry];
        if (cnt == 0)
          continue;
        switch (geometry)
        {
          case Geometry::Polytope::Type::Point:
          {
            os << EnSight6::ElementType::point << '\n';
            break;
          }
          case Geometry::Polytope::Type::Segment:
          {
            os << EnSight6::ElementType::bar2 << '\n';
            break;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            os << EnSight6::ElementType::tria3 << '\n';
            break;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            os << EnSight6::ElementType::quad4 << '\n';
            break;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            os << EnSight6::ElementType::tetra4 << '\n';
            break;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            os << EnSight6::ElementType::penta6 << '\n';
            break;
          }
        }
        os << cnt << '\n' << ess[attr][geometry].str();
      }
    }
  }
}
