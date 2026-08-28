/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRNORMALOFFSET_H
#define RODIN_ADAPTATION_WNGIRNORMALOFFSET_H

#include "WNGIRLoss.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Capped normal offset from an interface point to the target level set.
  class WNGIRNormalOffset
  {
    public:
      template <class PhiType, class GradType>
      /// @brief Constructs the WNGIR normal offset.
      WNGIRNormalOffset(const PhiType& phi, const GradType& grad,
        const Geometry::Point& p, const WNGIRParameters& params, Real sigma2)
      {
        constexpr Real epsG = Real(1e-12);
        const Real r = phi.getValue(p);
        const auto g = grad.getValue(p);
        const Real gNorm = std::sqrt(std::max(g.dot(g), Real(0)));
        m_degenerate = gNorm <= Real(1e-14);
        if (m_degenerate)
        {
          m_offset = 0;
          m_attenuation = 0;
          return;
        }
        m_normal = g / (gNorm + epsG);
        m_offset = -r / (gNorm + epsG);
        if (params.initialGuessCapH > Real(0) && params.h > Real(0))
        {
          const Real cap = params.initialGuessCapH * params.h;
          m_offset = std::clamp(m_offset, -cap, cap);
        }
        m_attenuation = WNGIRLoss(std::sqrt(sigma2)).getWeight(r);
      }

      /// @brief Whether degenerate.
      bool isDegenerate() const
      {
        return m_degenerate;
      }

      /// @brief The normal.
      const Math::SpatialVector<Real>& getNormal() const
      {
        return m_normal;
      }

      /// @brief The offset.
      Real getOffset() const
      {
        return m_offset;
      }

      /// @brief The attenuation.
      Real getAttenuation() const
      {
        return m_attenuation;
      }

    private:
      bool m_degenerate = false;
      Math::SpatialVector<Real> m_normal;
      Real m_offset = 0;
      Real m_attenuation = 0;
  };
}

#endif
