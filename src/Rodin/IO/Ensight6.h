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

namespace Rodin::IO::Ensight6
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
      case Keyword::point:
        return "point";
      case Keyword::bar2:
        return "bar2";
      case Keyword::bar3:
        return "bar3";
      case Keyword::tria3:
        return "tria3";
      case Keyword::tria6:
        return "tria6";
      case Keyword::quad4:
        return "quad4";
      case Keyword::quad8:
        return "quad8";
      case Keyword::tetra4:
        return "tetra4";
      case Keyword::tetra10:
        return "tetra10";
      case Keyword::pyramid5:
        return "pyramid5";
      case Keyword::pyramid13:
        return "pyramid13";
      case Keyword::hexa8:
        return "hexa8";
      case Keyword::hexa20:
        return "hexa20";
      case Keyword::penta6:
        return "penta6";
      case Keyword::penta15:
        return "penta15";
      case Keyword::block:
        return "block";
      case Keyword::iblanked:
        return "iblanked";
    }
    return nullptr;
  }

  inline
  bool operator==(const std::string& str, Keyword kw)
  {
    return str == toCharString(kw);
  }

  inline
  bool operator==(Keyword kw, const std::string& str)
  {
    return str == toCharString(kw);
  }

  inline
  bool operator==(Keyword kw, const char* str)
  {
    return strcmp(toCharString(kw), str) == 0;
  }

  inline
  bool operator==(const char* str, Keyword kw)
  {
    return strcmp(toCharString(kw), str) == 0;
  }

  inline
  bool operator!=(const char* str, Keyword kw)
  {
    return !operator==(str, kw);
  }

  inline
  bool operator!=(const std::string& str, Keyword kw)
  {
    return !operator==(str, kw);
  }

  inline
  bool operator!=(Keyword kw, const std::string& str)
  {
    return !operator==(str, kw);
  }

  inline
  std::ostream& operator<<(std::ostream& os, Keyword kw)
  {
    os << toCharString(kw);
    return os;
  }

  inline
  std::optional<Keyword> toKeyword(const char* str)
  {
    Keyword res;
    if (str == Keyword::node)
      res = Keyword::node;
    else if (str == Keyword::id)
      res = Keyword::id;
    else if (str == Keyword::off)
      res = Keyword::off;
    else if (str == Keyword::given)
      res = Keyword::given;
    else if (str == Keyword::assign)
      res = Keyword::assign;
    else if (str == Keyword::ignore)
      res = Keyword::ignore;
    else if (str == Keyword::element)
      res = Keyword::element;
    else if (str == Keyword::coordinates)
      res = Keyword::coordinates;
    else if (str == Keyword::part)
      res = Keyword::part;
    else if (str == Keyword::point)
      res = Keyword::point;
    else if (str == Keyword::bar2)
      res = Keyword::bar2;
    else if (str == Keyword::bar3)
      res = Keyword::bar3;
    else if (str == Keyword::tria3)
      res = Keyword::tria3;
    else if (str == Keyword::tria6)
      res = Keyword::tria6;
    else if (str == Keyword::quad4)
      res = Keyword::quad4;
    else if (str == Keyword::quad8)
      res = Keyword::quad8;
    else if (str == Keyword::tetra4)
      res = Keyword::tetra4;
    else if (str == Keyword::tetra10)
      res = Keyword::tetra10;
    else if (str == Keyword::pyramid5)
      res = Keyword::pyramid5;
    else if (str == Keyword::pyramid13)
      res = Keyword::pyramid13;
    else if (str == Keyword::hexa8)
      res = Keyword::hexa8;
    else if (str == Keyword::hexa20)
      res = Keyword::hexa20;
    else if (str == Keyword::penta6)
      res = Keyword::penta6;
    else if (str == Keyword::penta15)
      res = Keyword::penta15;
    else if (str == Keyword::block)
      res = Keyword::block;
    else if (str == Keyword::iblanked)
      res = Keyword::iblanked;
    else
      return {};
    assert(res == str);
    return res;
  }
}

#endif
