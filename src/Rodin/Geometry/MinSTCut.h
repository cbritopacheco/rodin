/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MINSTCUT_H
#define RODIN_GEOMETRY_MINSTCUT_H

#include <limits>
#include <vector>

#include "Rodin/Types.h"

namespace Rodin::Geometry
{
  /**
   * @brief Serial binary Potts classifier solved by an s-t min cut.
   *
   * Labels follow the level-set convention used by topology reconstruction:
   * - -1 is inside;
   * - +1 is outside.
   */
  class MinSTCut
  {
    public:
      static constexpr int Inside = -1;
      static constexpr int Outside = 1;
      static constexpr Index InvalidIndex =
        std::numeric_limits<Index>::max();

      struct Edge
      {
        Index first;
        Index second;
        Real capacity;
        Index index = InvalidIndex;
      };

      struct Result
      {
        std::vector<int> labels;
        std::vector<Index> insideCells;
        std::vector<Index> outsideCells;
        std::vector<Edge> cutEdges;
        Real energy = 0;
      };

      static Real getInsideCost(Real volume, Real moment) noexcept;

      static Real getOutsideCost(Real volume, Real moment) noexcept;

      Result classify(
          const std::vector<Real>& volumes,
          const std::vector<Real>& moments,
          const std::vector<Edge>& edges) const;

      Result solve(
          const std::vector<Real>& insideCosts,
          const std::vector<Real>& outsideCosts,
          const std::vector<Edge>& edges) const;
  };
}

#endif
