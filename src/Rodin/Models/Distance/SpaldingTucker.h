/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SpaldingTucker.h
 * @brief Spalding-Tucker normalization for level set functions.
 *
 * This file provides the SpaldingTucker class, which normalizes level set
 * functions to approximate signed distance functions using the Spalding-Tucker
 * method.
 */
#ifndef RODIN_MODELS_DISTANCE_SPALDINGTUCKER_H
#define RODIN_MODELS_DISTANCE_SPALDINGTUCKER_H

#include "Rodin/Variational/GridFunction.h"

namespace Rodin::Models::Distance
{
  /**
   * @brief Spalding-Tucker normalization for a level set function.
   *
   * This class implements the Spalding-Tucker normalization method, which
   * transforms an arbitrary level set function into an approximate signed
   * distance function. The transformation is given by:
   * @f[
   *   d(x) = \frac{2u(x)}{|\nabla u(x)| + \sqrt{|\nabla u(x)|^2 + 2|u(x)|}}
   * @f]
   * where @f$ u @f$ is the input level set function and @f$ d @f$ is the
   * normalized approximate signed distance.
   *
   * ## Mathematical Properties
   * - Preserves the zero level set: @f$ d(x) = 0 \iff u(x) = 0 @f$
   * - Preserves the sign: @f$ \text{sign}(d) = \text{sign}(u) @f$
   * - Better gradient normalization than simple Rvachev method
   * - More accurate near the interface
   *
   * ## References
   * @cite belyaev2015variational @cite spalding1994calculation @cite tucker1998assessment
   *
   * ## Usage Example
   * ```cpp
   * GridFunction u(fes);  // Some level set function
   * SpaldingTucker st;
   * auto dist = st(u);  // Normalized distance approximation
   * ```
   */
  class SpaldingTucker
  {
    public:
      /**
       * @brief Applies Spalding-Tucker normalization to a grid function.
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
        dist =
          [&](const Geometry::Point& p)
          {
            const Real u = gf(p);
            Variational::Grad grad(gf);
            const Real norm = grad.getValue(p).norm();
            return (2 * u) / (norm + Math::sqrt(norm * norm + 2 * Math::abs(u)));
          };
        return dist;
      }
  };
}

#endif

