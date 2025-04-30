#ifndef RODIN_IO_ENSIGHT6_H
#define RODIN_IO_ENSIGHT6_H

#include <iomanip>
#include <boost/bimap.hpp>
#include <boost/spirit/home/x3.hpp>

#include "Rodin/Types.h"
#include "Rodin/Alert.h"
#include "Rodin/Context.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Geometry/Types.h"

#include "ForwardDecls.h"
#include "MeshLoader.h"
#include "MeshPrinter.h"
#include "GridFunctionLoader.h"
#include "GridFunctionPrinter.h"

namespace Rodin::IO::EnSight6
{
  enum class Keyword
  {
    node,
    id,
    off,
    given,
    assign,
    ignore,
    element,
    coordinates,
    part,
    block,
    iblanked
  };

  inline
  constexpr
  const char* toCharString(Keyword kw)
  {
    switch (kw)
    {
      case Keyword::node:
        return "node";
      case Keyword::id:
        return "id";
      case Keyword::off:
        return "off";
      case Keyword::given:
        return "given";
      case Keyword::assign:
        return "assign";
      case Keyword::ignore:
        return "ignore";
      case Keyword::element:
        return "element";
      case Keyword::coordinates:
        return "coordinates";
      case Keyword::part:
        return "part";
      case Keyword::block:
        return "block";
      case Keyword::iblanked:
        return "iblanked";
    }
    return nullptr;
  }

  inline
  std::ostream& operator<<(std::ostream& os, Keyword kw)
  {
    os << toCharString(kw);
    return os;
  }

  enum class ElementType
  {
    point,
    bar2,
    bar3,
    tria3,
    tria6,
    quad4,
    quad8,
    tetra4,
    tetra10,
    pyramid5,
    pyramid13,
    hexa8,
    hexa20,
    penta6,
    penta15,
  };

  inline
  constexpr
  const char* toCharString(ElementType kw)
  {
    switch (kw)
    {
      case ElementType::point:
        return "point";
      case ElementType::bar2:
        return "bar2";
      case ElementType::bar3:
        return "bar3";
      case ElementType::tria3:
        return "tria3";
      case ElementType::tria6:
        return "tria6";
      case ElementType::quad4:
        return "quad4";
      case ElementType::quad8:
        return "quad8";
      case ElementType::tetra4:
        return "tetra4";
      case ElementType::tetra10:
        return "tetra10";
      case ElementType::pyramid5:
        return "pyramid5";
      case ElementType::pyramid13:
        return "pyramid13";
      case ElementType::hexa8:
        return "hexa8";
      case ElementType::hexa20:
        return "hexa20";
      case ElementType::penta6:
        return "penta6";
      case ElementType::penta15:
        return "penta15";
    }
    return nullptr;
  }

  inline
  std::ostream& operator<<(std::ostream& os, ElementType kw)
  {
    os << toCharString(kw);
    return os;
  }

  inline
  constexpr
  std::optional<ElementType> getGeometry(Geometry::Polytope::Type t)
  {
    switch (t)
    {
      case Geometry::Polytope::Type::Point:
        return ElementType::point;
      case Geometry::Polytope::Type::Segment:
        return ElementType::bar2;
      case Geometry::Polytope::Type::Triangle:
        return ElementType::tria3;
      case Geometry::Polytope::Type::Quadrilateral:
        return ElementType::quad4;
      case Geometry::Polytope::Type::Tetrahedron:
        return ElementType::tetra4;
      case Geometry::Polytope::Type::Wedge:
        return ElementType::penta6;
      default:
        return {};
    }
    assert(false);
    return {};
  }

  enum VariableType
  {
    scalar,
    vector
  };

  enum class Location
  {
    node,
    element
  };
}

namespace Rodin::IO
{
  template <>
  class MeshPrinter<FileFormat::ENSIGHT6, Context::Local>
    : public MeshPrinterBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = MeshPrinterBase<ContextType>;

      MeshPrinter(const ObjectType& mesh);

      void print(std::ostream& os) override;

      void printHeader(std::ostream& os);
      void printCoordinates(std::ostream& os);
      void printParts(std::ostream& os);

    private:
      std::string m_descriptionLine1;
      std::string m_descriptionLine2;
  };

  template <class FES>
  class GridFunctionPrinter<FileFormat::ENSIGHT6, FES>
    : public GridFunctionPrinterBase<FES>
  {
    public:
      using FESType = FES;

      using ObjectType = Variational::GridFunction<FESType>;

      using Parent = GridFunctionPrinterBase<FESType>;

      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      void print(std::ostream& os) override
      {
        printHeader(os);
        printData(os);
      }

      void printHeader(std::ostream& os)
      {}

      void printData(std::ostream& os)
      {}
  };

}

#endif
