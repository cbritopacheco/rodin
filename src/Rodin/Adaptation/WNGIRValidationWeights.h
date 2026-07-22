/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRVALIDATIONWEIGHTS_H
#define RODIN_ADAPTATION_WNGIRVALIDATIONWEIGHTS_H

#include <cmath>
#include <functional>

#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Positive weights for sampled nonlinear geometry functionals.
  class WNGIRValidationWeights
  {
    public:
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
        m_scale = absoluteWeight > Real(0)
          ? signedWeight / absoluteWeight
          : Real(0);
      }

      Real getWeight(std::size_t q) const
      {
        return m_scale * std::abs(m_qf.get().getWeight(q));
      }

    private:
      std::reference_wrapper<const QF::QuadratureFormulaBase> m_qf;
      Real m_scale;
  };
}

#endif
