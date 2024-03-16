#ifndef RODIN_GEOMETRY_TYPES_H
#define RODIN_GEOMETRY_TYPES_H

#include <vector>
#include <boost/range/adaptors.hpp>

#include "Rodin/Types.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /// Standard type for representing material attributes in a mesh.
  using Attribute = std::size_t;

  /// Represents the incidence of a polytope.
  using Incidence = std::vector<IndexSet>;
}

#endif
