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
      /// @brief Label value for cells inside the level set.
      static constexpr int Inside = -1;
      /// @brief Label value for cells outside the level set.
      static constexpr int Outside = 1;
      /// @brief Sentinel index denoting an absent edge/index.
      static constexpr Index InvalidIndex = std::numeric_limits<Index>::max();

      /// @brief A weighted graph edge between two cells.
      struct Edge
      {
          /// @brief Index of the first incident cell.
          Index first;
          /// @brief Index of the second incident cell.
          Index second;
          /// @brief Edge capacity (pairwise smoothing weight).
          Real capacity;
          /// @brief Optional index of this edge in the caller's edge list.
          Index index = InvalidIndex;
      };

      /// @brief Result of a classification: per-cell labels and cut data.
      struct Result
      {
          /// @brief Per-cell label (Inside or Outside).
          std::vector<int> labels;
          /// @brief Indices of the cells labeled Inside.
          std::vector<Index> insideCells;
          /// @brief Indices of the cells labeled Outside.
          std::vector<Index> outsideCells;
          /// @brief Edges crossing the cut (the interface skeleton).
          std::vector<Edge> cutEdges;
          /// @brief Total energy (cut cost) of the solution.
          Real energy = 0;
      };

      /**
       * @brief Optional classification policies layered on top of the bare
       * Potts graph cut.
       *
       * These options let callers control checkerboard suppression,
       * narrow-band restriction, and zero-level-aware facet weighting
       * without changing the underlying max-flow engine.
       */
      struct Options
      {
          /**
           * Multiplicative scale applied to every pairwise capacity
           * (Potts smoothing strength lambda).
           */
          Real lambdaScale = 1;

          /// Multiplicative scale applied to every unary cost.
          Real unaryScale = 1;

          /**
           * If `|moment[i]| >= farFieldThreshold` the cell `i` is pinned
           * to `sign(-moment[i])` by injecting a large unary on the
           * opposite terminal. Disabled when negative.
           */
          Real farFieldThreshold = -1;

          /**
           * Per-edge multiplier replacing `lambdaScale` on the i-th edge.
           * Use for zero-level-aware facet weighting (e.g. down-weight
           * pairwise term where |phi| is small so the cut prefers to
           * align with the actual interface).
           * Size must be `edges.size()` or `0` (disabled).
           */
          std::vector<Real> perEdgeLambda;

          /**
           * Per-cell free/fixed mask. `cellInBand[i] == false` pins cell
           * `i` to `sign(-moment[i])` via a large unary. Size must be
           * `volumes.size()` or `0` (disabled, all cells free).
           */
          std::vector<Boolean> cellInBand;
      };

      /**
       * @brief Returns the unary cost of labeling a cell Inside, given its
       * volume and signed moment.
       * @returns The unary cost of labeling a cell Inside, given its volume and signed moment.
       * @param volume Cell volume.
       * @param moment Integral of the field over the cell.
       */
      static Real getInsideCost(Real volume, Real moment) noexcept;

      /**
       * @brief Returns the unary cost of labeling a cell Outside, given its
       * volume and signed moment.
       * @returns The unary cost of labeling a cell Outside, given its volume and signed moment.
       * @param volume Cell volume.
       * @param moment Integral of the field over the cell.
       */
      static Real getOutsideCost(Real volume, Real moment) noexcept;

      /**
       * @brief Classifies cells into Inside/Outside via a Potts min s-t cut
       * built from per-cell volumes and moments.
       * @param volumes Cell volumes.
       * @param moments Cell integrals of the field.
       * @param edges Edges carrying the interface penalty.
       * @returns Cell labels and the value of the minimum-cut objective.
       */
      Result classify(const std::vector<Real>& volumes, const std::vector<Real>& moments,
        const std::vector<Edge>& edges) const;

      /**
       * @brief Classifies cells with additional @ref Options (narrow-band
       * restriction, far-field pinning, per-edge weighting).
       * @param volumes Cell volumes.
       * @param moments Cell integrals of the field.
       * @param edges Edges carrying the interface penalty.
       * @param options Options controlling the classification.
       * @returns Cell labels and the value of the minimum-cut objective.
       */
      Result classify(const std::vector<Real>& volumes, const std::vector<Real>& moments,
        const std::vector<Edge>& edges, const Options& options) const;

      /**
       * @brief Solves the s-t min cut directly from precomputed unary costs
       * and edges.
       * @param insideCosts Costs of assigning each cell to the inside region.
       * @param outsideCosts Costs of assigning each cell to the outside region.
       * @param edges Edges carrying the interface penalty.
       * @returns Cell labels and the value of the minimum-cut objective.
       */
      Result solve(const std::vector<Real>& insideCosts,
        const std::vector<Real>& outsideCosts, const std::vector<Edge>& edges) const;
  };
}

#endif
