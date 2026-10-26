/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_QF_REPEATEDFORMULA_H
#define RODIN_TESTS_QF_REPEATEDFORMULA_H
#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Tests
{
  /**
   * @brief Publishes distinct formula keys with identical numerical samples.
   *
   * A repetition index replaces the former use of cloned formula identities in
   * cache-eviction tests and working-set benchmarks. The source is itself
   * persistent; generation and publication are outside the measured region.
   */
  class RepeatedQuadratureFormula final : public ::Rodin::QF::QuadratureFormulaBase
  {
    public:
      /**
       * @param formula Persistent source samples.
       * @param repetition Distinguishes otherwise identical test or benchmark inputs.
       */
      RepeatedQuadratureFormula(
        const ::Rodin::QF::QuadratureFormula& formula, [[maybe_unused]] size_t repetition)
        : m_formula(formula)
      {}

      /// @returns Borrowed arrays from the persistent source formula.
      Optional<Samples> getSamples() const override
      {
        return m_formula.get().getSamples();
      }

      /// @returns Source sample count.
      size_t getSize() const override
      {
        return m_formula.get().getSize();
      }

      /**
       * @param index Sample index.
       * @returns Source weight.
       */
      Real getWeight(size_t index) const override
      {
        return m_formula.get().getWeight(index);
      }

      /**
       * @param index Sample index.
       * @returns Source reference coordinates.
       */
      const Math::SpatialPoint& getPoint(size_t index) const override
      {
        return m_formula.get().getPoint(index);
      }

    private:
      std::reference_wrapper<const ::Rodin::QF::QuadratureFormula> m_formula;
  };
}
#endif
