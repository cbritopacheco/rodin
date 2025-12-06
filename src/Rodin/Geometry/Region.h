/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_REGION_H
#define RODIN_GEOMETRY_REGION_H

/**
 * @file
 * @brief Region enumeration for identifying mesh subregions.
 */

namespace Rodin::Geometry
{
  /**
   * @brief Enumeration of standard mesh regions.
   *
   * This enumeration defines commonly used subregions of a mesh that
   * can be referenced in variational formulations and boundary conditions.
   */
  enum class Region
  {
    Cells,      ///< Interior cells (elements) of the mesh
    Faces,      ///< All faces (d-1 dimensional polytopes) in the mesh
    Boundary,   ///< Boundary faces (faces incident to only one cell)
    Interface   ///< Interior faces (faces incident to two cells)
  };
}

#endif
