/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_DISTANCE_RVACHEV_H
#define RODIN_MODELS_DISTANCE_RVACHEV_H

#include <utility>

#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/GridFunction.h"

namespace Rodin::Models::Distance
{
  /**
   * @brief Rvachev normalization for a level set function.
   * @cite rvachev1974methods @cite rvachev2001transfinite
   */
  class Rvachev
  {
    public:
      template <class FES>
      auto operator()(const Variational::GridFunction<FES>& gf)
      {
        Variational::GridFunction dist(gf.getFiniteElementSpace());
        Math::SpatialVector<Scalar> gu;
        dist =
          [&](const Geometry::Point& p)
          {
            const Scalar u = gf(p);
            Variational::Grad grad(gf);
            grad.getValue(gu, p);
            return u / Math::sqrt(u * u + gu.squaredNorm());
          };
        return dist;
      }
  };
}

#endif


