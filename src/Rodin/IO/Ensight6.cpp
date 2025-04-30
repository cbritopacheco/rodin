#include "Ensight6.h"

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
      os << it->getIndex() << ' ';
      const auto& x = it->getCoordinates();
      for (int i = 0; i < x.size(); i++)
        os << std::setw(12) << x(i);
      os << '\n';
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
              << it->getIndex() << ' ' << it->getIndex() << '\n';
            break;
          }
          case Geometry::Polytope::Type::Segment:
          {
            ess[it->getAttribute()][geometry]
              << ' ' << it->getIndex() << ' ' << vertices(0) << ' ' << vertices(1) << '\n';
            break;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            ess[it->getAttribute()][geometry]
              << ' ' << it->getIndex() << ' ' << vertices(0) << ' ' << vertices(1) << ' ' << vertices(2) << '\n';
            break;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            ess[it->getAttribute()][geometry]
              << ' ' << it->getIndex() << ' ' << vertices(0) << ' ' << vertices(1) << ' '
              << vertices(3) << ' ' << vertices(2) << '\n';
            break;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            ess[it->getAttribute()][geometry]
              << ' ' << it->getIndex() << ' ' << vertices(0) << ' ' << vertices(1) << ' '
              << vertices(2) << ' ' << vertices(3) << '\n';
            break;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            ess[it->getAttribute()][geometry]
              << ' ' << it->getIndex() << ' ' << vertices(0) << ' ' << vertices(1) << ' '
              << vertices(2) << ' ' << vertices(3) << ' ' << vertices(4) << ' ' << vertices(5) << '\n';
            break;
          }
        }
      }
    }

    for (const auto& attr : attributes)
    {
      os << EnSight6::Keyword::part << ' ' << attr << '\n' << "Rodin::Geometry::Attribute" << '\n';
      for (const auto& geometry : Geometry::Polytope::Types)
      {
        const size_t count = mesh.getPolytopeCount(geometry);
        if (count == 0)
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
        os << count << '\n' << ess[attr][geometry].str();
      }
    }
  }
}
