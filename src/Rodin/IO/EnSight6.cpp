#include "EnSight6.h"

#include "Rodin/Configure.h"

namespace Rodin::IO
{
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
      const auto as = mesh.getAttributeIndex().getAttributes(d);
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
