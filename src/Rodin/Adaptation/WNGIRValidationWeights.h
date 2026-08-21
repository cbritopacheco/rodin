/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRVALIDATIONWEIGHTS_H
#define RODIN_ADAPTATION_WNGIRVALIDATIONWEIGHTS_H

#include <cmath>
#include <functional>

#include "Rodin/QF/QuadratureFormula.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Positive weights for sampled nonlinear geometry functionals.
  class WNGIRValidationWeights
  {
    public:
      /// @brief Constructs the w n g i r validation weights.
      explicit WNGIRValidationWeights(const QF::QuadratureFormulaBase& qf)
        : m_qf(qf)
      {
        Real signedWeight = 0;
        Real absoluteWeight = 0;
        for (std::size_t q = 0; q < qf.getSize(); ++q)
        {
          signedWeight += qf.getWeight(q);
          absoluteWeight += std::abs(qf.getWeight(q));
        }
        m_scale = absoluteWeight > Real(0) ? signedWeight / absoluteWeight : Real(0);
      }

      /// @brief The weight.
      Real getWeight(std::size_t q) const
      {
        return m_scale * std::abs(m_qf.get().getWeight(q));
      }

      /**
       * @brief Correction that converts the signed quadrature weight at an
       * integration point into its positive, measure-preserving counterpart.
       *
       * Multiplying an integrand by this value before standard Rodin assembly
       * gives
       * @f[
       *   w_q\,c_q=\frac{\sum_i w_i}{\sum_i|w_i|}|w_q|.
       * @f]
       * This is required for sampled non-negative functionals, such as the
       * robust energy and its metric, for which cancellation by a signed
       * polynomial cubature rule has no variational interpretation.
       */
      static Real getCorrection(const Variational::IntegrationPoint& ip)
      {
        const auto* qf = ip.getQuadratureFormula();
        assert(qf);
        return getCorrection(*qf, ip.getIndex());
      }

      /// @brief Correction for one entry of a signed quadrature formula.
      static Real getCorrection(
        const QF::QuadratureFormulaBase& qf, std::size_t quadraturePoint)
      {
        thread_local const QF::QuadratureFormulaBase* cachedFormula = nullptr;
        thread_local Real cachedScale = Real(0);
        if (cachedFormula != &qf)
        {
          Real signedWeight = 0;
          Real absoluteWeight = 0;
          for (std::size_t q = 0; q < qf.getSize(); ++q)
          {
            signedWeight += qf.getWeight(q);
            absoluteWeight += std::abs(qf.getWeight(q));
          }
          cachedFormula = &qf;
          cachedScale =
            absoluteWeight > Real(0) ? signedWeight / absoluteWeight : Real(0);
        }
        const Real weight = qf.getWeight(quadraturePoint);
        if (weight == Real(0))
          return Real(0);
        return cachedScale * std::abs(weight) / weight;
      }

    private:
      std::reference_wrapper<const QF::QuadratureFormulaBase> m_qf;
      Real m_scale;
  };
}

#endif
