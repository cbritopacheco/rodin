/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_SWIFT_LOSS_H
#define RODIN_ADAPTATION_SWIFT_LOSS_H

#include <cassert>
#include <cmath>

#include "Rodin/Types.h"

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Fixed-scale Welsch loss with consistent influence and curvature.
   *
   * For a residual @f$r@f$, the influence and weight satisfy
   * @f$\rho'(r)=w(r)r@f$. The same loss is used by the objective, first
   * variation, and directional curvature throughout one nonlinear solve.
   * The canonical observation metric is unweighted squared-fit curvature.
   */
  class Loss
  {
    public:
      /**
       * @brief Constructs a Welsch loss with positive scale.
       * @param scale Scale controlling the loss function.
       */
      explicit Loss(Real scale)
        : m_scale2(scale * scale)
      {
        assert(scale > Real(0));
      }

      /**
       * @brief Returns the squared robust residual scale.
       * @returns The fixed squared scale.
       */
      Real getScaleSquared() const noexcept
      {
        return m_scale2;
      }

      /**
       * @brief Evaluates @f$\rho(r)@f$.
       * @param residual Signed level-set residual.
       * @returns The Welsch loss value.
       */
      Real getValue(Real residual) const
      {
        const Real s2 = residual * residual / m_scale2;
        return -Real(0.5) * m_scale2 * std::expm1(-s2);
      }

      /**
       * @brief Evaluates the weight @f$w(r)=\rho'(r)/r@f$.
       * @param residual Signed level-set residual.
       * @returns The robust influence weight, including its limit at zero.
       */
      Real getWeight(Real residual) const
      {
        const Real s2 = residual * residual / m_scale2;
        return std::exp(-s2);
      }

      /**
       * @brief Evaluates the influence @f$\rho'(r)@f$.
       * @param residual Residual to evaluate.
       * @returns The influence.
       */
      Real getInfluence(Real residual) const
      {
        return getWeight(residual) * residual;
      }

      /**
       * @brief Second derivative in the scalar residual, without level-set Hessian.
       * @param residual Signed level-set residual.
       * @returns The scalar Welsch curvature, which may be negative.
       */
      Real getCurvature(Real residual) const
      {
        return getWeight(residual) * (Real(1) - Real(2) * residual * residual / m_scale2);
      }

    private:
      Real m_scale2;
  };
}

#endif
