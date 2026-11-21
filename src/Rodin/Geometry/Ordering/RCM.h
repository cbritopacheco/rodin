/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ORDERING_RCM_H
#define RODIN_GEOMETRY_ORDERING_RCM_H

#include "Rodin/Geometry/Mesh.h"

namespace Rodin::Geometry::Ordering
{
  /**
   * @brief Deterministic reverse Cuthill-McKee ordering of local mesh vertices.
   *
   * Vertices sharing a maximal-dimensional cell form a clique in the graph.
   * Each connected component starts from a pseudo-peripheral vertex. Initial
   * seeds, last-level candidates and BFS neighbors are selected by increasing
   * degree, with original vertex index breaking ties. The BFS ordering is
   * reversed separately in each component; components are visited by their
   * minimum (degree, original index) seed. Isolated vertices are included.
   *
   * This is a bandwidth-reduction heuristic for the vertex graph, not an
   * ordering of high-order or mixed-space degrees of freedom. No derived
   * incidence is required or computed, and the input mesh is unchanged.
   */
  class RCM
  {
    public:
      /**
       * @brief Computes the reverse Cuthill-McKee vertex ordering.
       * @param[in] mesh Local mesh supplying cell-to-vertex lists.
       * @returns Permutation with newToOld[newIndex] equal to the old vertex index.
       */
      IndexVector operator()(const Mesh<Context::Local>& mesh) const;
  };
}

#endif
