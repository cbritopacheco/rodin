#ifndef RODIN_GEOMETRY_REGION_H
#define RODIN_GEOMETRY_REGION_H

namespace Rodin::Geometry
{
  enum class Region
  {
    Cells,
    Faces,
    Boundary,
    Interface
  };
}

#endif
