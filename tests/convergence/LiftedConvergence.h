/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_CONVERGENCE_LIFTED_CONVERGENCE_H
#define RODIN_TESTS_CONVERGENCE_LIFTED_CONVERGENCE_H

#include <gtest/gtest.h>
#include "LiftedErrorNorm.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Shared assertions for represented and lifted scalar field studies.
   * @par Contract
   * Separate field and geometry defects must add pointwise. The triangle and
   * reverse-triangle inequalities are checked for their norms, not equality.
   * Every adjacent interval is checked. The caller supplies independent
   * physical data and a regular geometry map; finite-resolution rate windows
   * are acceptance policies, not asymptotic theorems.
   */
  class LiftedConvergence
  {
    public:
      using Components = std::array<ErrorNorms, 4>;
      static constexpr Real L2Margin = 0.55, H1Margin = 0.45;
      static constexpr Real RoundoffTolerance = 1e-11, SensitivityTolerance = 1e-6;

      static Components components(
        const ErrorNorms& represented, const LiftedErrorNorm::Result& lifted)
      {
        return {represented, lifted.field, lifted.geometry, lifted.total};
      }

      static void expectDecomposition(const LiftedErrorNorm::Result& error)
      {
        for (const auto& norm : {error.field, error.geometry, error.total})
          EXPECT_TRUE(norm.isFinite());
        for (const auto& values :
          {std::array{error.field.getL2(), error.geometry.getL2(), error.total.getL2()},
            std::array{error.field.getH1Seminorm(), error.geometry.getH1Seminorm(),
              error.total.getH1Seminorm()}})
        {
          EXPECT_LE(values[2], values[0] + values[1] + RoundoffTolerance);
          EXPECT_GE(values[2], std::abs(values[0] - values[1]) - RoundoffTolerance);
        }
      }

      void append(
        Real h, const ErrorNorms& represented, const LiftedErrorNorm::Result& lifted)
      {
        expectDecomposition(lifted);
        const auto errors = components(represented, lifted);
        for (size_t component = 0; component < errors.size(); ++component)
        {
          ASSERT_TRUE(errors[component].isFinite());
          ASSERT_GT(errors[component].getL2(), 0);
          ASSERT_GT(errors[component].getH1Seminorm(), 0);
          m_histories[component].append(h, errors[component]);
        }
      }

      void expectRates(size_t fieldDegree, size_t geometryDegree = 2) const
      {
        for (size_t component = 0; component < m_histories.size(); ++component)
        {
          const auto& history = m_histories[component];
          ASSERT_GE(history.getSize(), 3);
          const size_t degree = component < 2 ? fieldDegree
            : component == 2                  ? geometryDegree
                                              : std::min(fieldDegree, geometryDegree);
          for (size_t i = 1; i < history.getSize(); ++i)
          {
            const auto& coarse = history.getSample(i - 1).error;
            const auto& fine = history.getSample(i).error;
            const auto rate = history.getAlgebraicRates(i);
            SCOPED_TRACE(::testing::Message()
              << "component=" << component << " interval=" << i
              << " L2=" << coarse.getL2() << " -> " << fine.getL2()
              << " H1=" << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm()
              << " rates=" << rate.getL2() << "," << rate.getH1Seminorm());
            EXPECT_GT(coarse.getL2(), fine.getL2());
            EXPECT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
            EXPECT_GT(rate.getL2(), degree + 1 - L2Margin);
            EXPECT_LT(rate.getL2(), degree + 1 + L2Margin);
            EXPECT_GT(rate.getH1Seminorm(), degree - H1Margin);
            EXPECT_LT(rate.getH1Seminorm(), degree + H1Margin);
          }
        }
      }

      static void expectSensitivity(const Components& base, const Components& refined)
      {
        for (size_t component = 0; component < base.size(); ++component)
          for (const auto& pair :
            {std::pair{base[component].getL2(), refined[component].getL2()},
              std::pair{
                base[component].getH1Seminorm(), refined[component].getH1Seminorm()}})
          {
            SCOPED_TRACE(::testing::Message() << "component=" << component);
            ASSERT_GT(pair.first, 0);
            ASSERT_TRUE(std::isfinite(pair.second));
            EXPECT_LT(std::abs(pair.second / pair.first - 1), SensitivityTolerance);
          }
      }

    private:
      std::array<ErrorHistory, 4> m_histories;
  };
}
#endif
