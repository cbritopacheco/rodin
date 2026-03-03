/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Rvachev.h
 * @brief Rvachev normalization for level set functions.
 *
 * This file provides the Rvachev class, which normalizes level set functions
 * to approximate signed distance functions using Rvachev's method.
 */
#ifndef RODIN_MODELS_DISTANCE_RVACHEV_H
#define RODIN_MODELS_DISTANCE_RVACHEV_H

#include <utility>

#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/GridFunction.h"

namespace Rodin::Distance
{
  /**
   * @brief Rvachev normalization for a level set function.
   *
   * This class implements the Rvachev normalization method, which transforms
   * an arbitrary level set function into an approximate signed distance function.
   * The transformation is given by:
   * @f[
   *   d(x) = \frac{u(x)}{\sqrt{u(x)^2 + |\nabla u(x)|^2}}
   * @f]
   * where @f$ u @f$ is the input level set function and @f$ d @f$ is the
   * normalized approximate signed distance.
   *
   * ## Mathematical Properties
   * - Preserves the zero level set: @f$ d(x) = 0 \iff u(x) = 0 @f$
   * - Preserves the sign: @f$ \text{sign}(d) = \text{sign}(u) @f$
   * - Gradient magnitude approaches 1 near the interface
   *
   * ## References
   * @cite rvachev1974methods @cite rvachev2001transfinite
   *
   * ## Usage Example
   * ```cpp
   * GridFunction u(fes);  // Some level set function
   * Rvachev rvachev;
   * auto dist = rvachev(u);  // Normalized distance approximation
   * ```
   */
  class Rvachev
  {
    public:
      /**
       * @brief Applies Rvachev normalization to a grid function.
       *
       * @tparam FES Finite element space type
       * @tparam Data Data storage type
       * @param[in] gf Input level set grid function
       * @return Normalized grid function approximating signed distance
       */
      template <class FES, class Data>
      auto operator()(const Variational::GridFunction<FES, Data>& gf)
      {
        Variational::GridFunction dist(gf.getFiniteElementSpace());
        Math::SpatialVector<Real> gu;
        dist =
          [&](const Geometry::Point& p)
          {
            const Real u = gf(p);
            Variational::Grad grad(gf);
            grad.getValue(gu, p);
            return u / Math::sqrt(u * u + gu.squaredNorm());
          };
        return dist;
      }
  };
}

#endif


