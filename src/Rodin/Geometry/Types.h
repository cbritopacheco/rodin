/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_TYPES_H
#define RODIN_GEOMETRY_TYPES_H

/**
 * @file
 * @brief Type definitions for the Rodin::Geometry module.
 */

#include <vector>
#include <boost/range/adaptors.hpp>

#include "Rodin/Types.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Standard type for representing material attributes in a mesh.
   *
   * Attributes are used to mark different regions or materials in a mesh.
   * For example, in a multi-material problem, each material region would
   * have a distinct attribute value.
   */
  using Attribute = std::size_t;

  /**
   * @brief Represents the incidence relation between polytopes.
   *
   * The incidence relation @f$ d \longrightarrow d' @f$ stores for each
   * polytope of dimension @f$ d @f$ the indices of all incident polytopes
   * of dimension @f$ d' @f$.
   *
   * @see Connectivity
   */
  using Incidence = std::vector<std::vector<Index>>;
}

#endif
