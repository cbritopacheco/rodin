/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_DISTANCE_SPALDINGTUCKER_H
#define RODIN_MODELS_DISTANCE_SPALDINGTUCKER_H

#include <utility>

#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/GridFunction.h"

namespace Rodin::Models::Distance
{
  /**
   * @brief Spalding-Tucker normalization for a level set function.
   * @cite belyaev2015variational @cite spalding1994calculation, @cite tucker1998assessment.
   */
  class SpaldingTucker
  {
    public:
      template <class FES>
      auto operator()(const Variational::GridFunction<FES>& gf)
      {
        Variational::GridFunction dist(gf.getFiniteElementSpace());
        Math::SpatialVector<Real> gu;
        dist =
          [&](const Geometry::Point& p)
          {
            const Real u = gf(p);
            Variational::Grad grad(gf);
            grad.getValue(gu, p);
            const Real norm = gu.norm();
            return (2 * u) / (norm + Math::sqrt(norm * norm + 2 * Math::abs(u)));
          };
        return dist;
      }
  };
}

#endif

