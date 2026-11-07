/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_LAPLACIAN_SMOOTHING_H
#define KELVIN_BALL_LAPLACIAN_SMOOTHING_H

#include "Common.h"

namespace KelvinBall
{
  /**
   * @brief Constrained Laplacian smoothing of the affine tetrahedral chamber.
   *
   * Architecture: unique edge neighbours define a Jacobi displacement toward
   * their arithmetic mean. Gamma and Outer are pinned; cut displacements are
   * restricted to the intersection of their incident planes. Two relaxed passes
   * are attempted, each with bounded backtracking. A candidate is accepted only
   * when cell orientations and minimum and mean mean-ratio quality are preserved.
   * Thus smoothing changes neither the fitted interface nor the material topology.
   */
  class LaplacianSmoothing
  {
    public:
      /**
       * @brief Smooths a reconstructed mesh before rebuilding its spaces.
       * @param[in,out] mesh Affine tetrahedral chamber with labelled cut faces.
       * @returns Number of accepted smoothing passes (zero leaves geometry unchanged).
       */
      size_t execute(Mesh& mesh) const;
  };
}

#endif
