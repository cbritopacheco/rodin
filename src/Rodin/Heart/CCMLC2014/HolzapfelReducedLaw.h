/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file HolzapfelReducedLaw.h
 * @brief Reduced Holzapfel-like passive-energy law used in CCMLC2014.
 */
#ifndef RODIN_HEART_CCMLC2014_HOLZAPFELREDUCEDLAW_H
#define RODIN_HEART_CCMLC2014_HOLZAPFELREDUCEDLAW_H

#include <cmath>

#include "Rodin/Heart/CCMLC2014/PassiveEnergy.h"

namespace Rodin::Heart
{
  /**
   * @brief Passive-energy evaluator based on reduced invariants.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  class HolzapfelReducedLaw
  {
    public:
      /**
       * @brief Parameters of the reduced passive-energy law.
       */
      struct Parameters
      {
        Scalar mu1 = Scalar(0); ///< Linear coefficient for @f$ J_1 @f$.
        Scalar mu2 = Scalar(0); ///< Linear coefficient for @f$ J_2 @f$.
        Scalar C0 = Scalar(0);  ///< Exponential scale for @f$ J_1 @f$.
        Scalar C1 = Scalar(0);  ///< Exponential exponent factor for @f$ J_1 @f$.
        Scalar C2 = Scalar(0);  ///< Exponential scale for @f$ J_4 @f$.
        Scalar C3 = Scalar(0);  ///< Exponential exponent factor for @f$ J_4 @f$.
      };

      /**
       * @brief Constructs the law with a parameter set.
       * @param params Reduced passive-energy parameters.
       */
      explicit HolzapfelReducedLaw(const Parameters& params = Parameters())
        : m_params(params)
      {}

      /**
       * @brief Evaluate first- and second-order derivatives at reduced invariants.
       *
       * @param[in] I Reduced invariants.
       * @returns Gradient and Hessian of passive energy.
       */
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

      /**
       * @brief Returns the parameter set used by this law.
       * @return Const parameter reference.
       */
      const Parameters& getParameters() const noexcept
      {
        return m_params;
      }

    private:
      Parameters m_params;
  };
}

#endif
