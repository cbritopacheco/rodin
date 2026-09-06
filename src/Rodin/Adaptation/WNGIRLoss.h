/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRLOSS_H
#define RODIN_ADAPTATION_WNGIRLOSS_H

#include <cassert>
#include <cmath>

#include "Rodin/Types.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Fixed-scale Welsch loss and its consistent IRLS weight.
   *
   * For a residual @f$r@f$, the influence and weight satisfy
   * @f$\rho'(r)=w(r)r@f$. The same loss is used by the objective, first
   * variation, and observation metric throughout one nonlinear solve.
   */
  class WNGIRLoss
  {
    public:
      /// @brief Constructs a Welsch loss with positive scale.
      explicit WNGIRLoss(Real scale)
        : m_scale2(scale * scale)
      {
        assert(scale > Real(0));
      }

      /// @brief Evaluates @f$\rho(r)@f$.
      Real getValue(Real residual) const
      {
        const Real s2 = residual * residual / m_scale2;
        return Real(0.5) * m_scale2 * (Real(1) - std::exp(-s2));
      }

      /// @brief Evaluates the Welsch weight @f$w(r)=\rho'(r)/r@f$.
      Real getWeight(Real residual) const
      {
        const Real s2 = residual * residual / m_scale2;
        return std::exp(-s2);
      }

      /// @brief Evaluates the influence @f$\rho'(r)@f$.
      Real getInfluence(Real residual) const
      {
        return getWeight(residual) * residual;
      }

    private:
      Real m_scale2;
  };
}

#endif
