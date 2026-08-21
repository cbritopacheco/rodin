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
  /// @brief Robust loss used by the interface registration objective.
  enum class WNGIRLossType
  {
    Welsch,
    Cauchy,
    PseudoHuber
  };

  /**
   * @brief Fixed-scale robust loss and its consistent IRLS weight.
   *
   * For a residual @f$r@f$, the influence and weight satisfy
   * @f$\rho'(r)=w(r)r@f$. The same loss is used by the objective, first
   * variation, and observation metric throughout one nonlinear solve.
   */
  class WNGIRLoss
  {
    public:
      WNGIRLoss(WNGIRLossType type, Real scale)
        : m_type(type),
          m_scale2(scale * scale)
      {
        assert(scale > Real(0));
      }

      /// @brief Evaluates @f$\rho(r)@f$.
      Real getValue(Real residual) const
      {
        const Real s2 = residual * residual / m_scale2;
        switch (m_type)
        {
          case WNGIRLossType::Welsch:
            return Real(0.5) * m_scale2 * (Real(1) - std::exp(-s2));
          case WNGIRLossType::Cauchy:
            return Real(0.5) * m_scale2 * std::log1p(s2);
          case WNGIRLossType::PseudoHuber:
            return m_scale2 * (std::sqrt(Real(1) + s2) - Real(1));
        }
        return Real(0);
      }

      /// @brief Evaluates the IRLS weight @f$w(r)=\rho'(r)/r@f$.
      Real getWeight(Real residual) const
      {
        const Real s2 = residual * residual / m_scale2;
        switch (m_type)
        {
          case WNGIRLossType::Welsch:
            return std::exp(-s2);
          case WNGIRLossType::Cauchy:
            return Real(1) / (Real(1) + s2);
          case WNGIRLossType::PseudoHuber:
            return Real(1) / std::sqrt(Real(1) + s2);
        }
        return Real(0);
      }

      /// @brief Evaluates the influence @f$\rho'(r)@f$.
      Real getInfluence(Real residual) const
      {
        return getWeight(residual) * residual;
      }

      /// @brief Returns the selected loss family.
      WNGIRLossType getType() const noexcept
      {
        return m_type;
      }

      /// @brief Returns @f$\sigma^2@f$.
      Real getScaleSquared() const noexcept
      {
        return m_scale2;
      }

    private:
      WNGIRLossType m_type;
      Real m_scale2;
  };
}

#endif
