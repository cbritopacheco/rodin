/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_HOLZAPFELREDUCEDLAW_H
#define RODIN_HEART_CCMLC2014_HOLZAPFELREDUCEDLAW_H

#include <cmath>

#include "Rodin/Heart/CCMLC2014/PassiveEnergy.h"

namespace Rodin::Heart
{
  template <class Scalar>
  class HolzapfelReducedLaw
  {
    public:
      struct Parameters
      {
        Scalar mu1 = Scalar(0);
        Scalar mu2 = Scalar(0);
        Scalar C0 = Scalar(0);
        Scalar C1 = Scalar(0);
        Scalar C2 = Scalar(0);
        Scalar C3 = Scalar(0);
      };

      explicit HolzapfelReducedLaw(const Parameters& params = Parameters())
        : m_params(params)
      {}

      PassiveEnergyDerivatives<Scalar> evaluate(const ReducedInvariants<Scalar>& I) const
      {
        PassiveEnergyDerivatives<Scalar> out;

        const Scalar x1 = I.J1 - Scalar(3);
        const Scalar x4 = I.J4 - Scalar(1);

        const Scalar e1 = std::exp(m_params.C1 * x1 * x1);
        const Scalar e4 = std::exp(m_params.C3 * x4 * x4);

        out.grad.dW_dJ1 =
            m_params.mu1
          + Scalar(2) * m_params.C0 * m_params.C1 * x1 * e1;
        out.grad.dW_dJ2 = m_params.mu2;
        out.grad.dW_dJ4 =
          Scalar(2) * m_params.C2 * m_params.C3 * x4 * e4;

        out.hess.d2W_dJ1dJ1 =
            Scalar(2) * m_params.C0 * m_params.C1 * e1
          + Scalar(4) * m_params.C0 * m_params.C1 * m_params.C1 * x1 * x1 * e1;
        out.hess.d2W_dJ1dJ2 = Scalar(0);
        out.hess.d2W_dJ1dJ4 = Scalar(0);

        out.hess.d2W_dJ2dJ2 = Scalar(0);
        out.hess.d2W_dJ2dJ4 = Scalar(0);

        out.hess.d2W_dJ4dJ4 =
            Scalar(2) * m_params.C2 * m_params.C3 * e4
          + Scalar(4) * m_params.C2 * m_params.C3 * m_params.C3 * x4 * x4 * e4;

        return out;
      }

      const Parameters& getParameters() const noexcept
      {
        return m_params;
      }

    private:
      Parameters m_params;
  };
}

#endif
