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
  /** @brief Shared assertions for represented and lifted field studies.
   * @par Contract
   * Separate field and geometry defects must add pointwise. The triangle and
   * reverse-triangle inequalities are checked for their norms, not equality.
   * Every adjacent interval is checked. The caller supplies independent
   * physical data and a regular geometry map; finite-resolution rate windows
   * are acceptance policies, not asymptotic theorems.
   * @par Architecture
   * Four histories distinguish represented-domain, lifted-field, geometry,
   * and total errors. Ordinary approximation studies populate all histories.
   * Representable-field studies certify the first two against a caller-owned
   * absolute budget and populate only geometry and total histories. Both
   * paths require three levels and check every adjacent interval; relative
   * sensitivity is never applied to a vanishing field error.
   */
  class LiftedConvergence
  {
    public:
      using Components = std::array<ErrorNorms, 4>;
      static constexpr Real L2Margin = 0.55, H1Margin = 0.45;
      static constexpr Real RoundoffTolerance = 1e-11, SensitivityTolerance = 1e-6;

      LiftedConvergence()
        : m_histories()
      {}

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
        expectRatesFrom(0, fieldDegree, geometryDegree);
      }

      /** @brief Certify different field/geometry orders without assuming dominance.
       * Separate components retain their two-sided rate windows. For each norm,
       * the total error must decrease and satisfy the next-level sum envelope
       * @f$T_f\le F_c\rho^{s_F-\delta}+G_c\rho^{s_G-\delta}@f$,
       * where @f$\rho=h_f/h_c@f$. Cancellation and a large field constant do
       * not imply a two-sided total-error rate at the slower geometry order.
       */
      void expectMixedRates(size_t fieldDegree, size_t geometryDegree) const
      {
        expectRatesFrom(0, fieldDegree, geometryDegree, 3);
        const auto& total = m_histories[3];
        ASSERT_GE(total.getSize(), 3);
        for (size_t i = 1; i < total.getSize(); ++i)
        {
          SCOPED_TRACE(::testing::Message() << "mixed total interval=" << i);
          const auto& coarse = total.getSample(i - 1);
          const auto& fine = total.getSample(i);
          const Real ratio = fine.parameter / coarse.parameter;
          ASSERT_GT(ratio, 0);
          ASSERT_LT(ratio, 1);
          const auto& field = m_histories[1].getSample(i - 1).error;
          const auto& geometry = m_histories[2].getSample(i - 1).error;
          EXPECT_GT(coarse.error.getL2(), fine.error.getL2());
          EXPECT_GT(coarse.error.getH1Seminorm(), fine.error.getH1Seminorm());
          EXPECT_LE(fine.error.getL2(),
            field.getL2() * std::pow(ratio, fieldDegree + 1 - L2Margin) +
              geometry.getL2() * std::pow(ratio, geometryDegree + 1 - L2Margin) +
              RoundoffTolerance);
          EXPECT_LE(fine.error.getH1Seminorm(),
            field.getH1Seminorm() * std::pow(ratio, fieldDegree - H1Margin) +
              geometry.getH1Seminorm() * std::pow(ratio, geometryDegree - H1Margin) +
              RoundoffTolerance);
        }
      }

      /** @brief Check a representable field independently of geometry error.
       * The caller supplies an absolute field-error budget in the norm's units.
       * Field norms may vanish; geometry and total norms need not vanish.
       */
      static void expectRepresentable(const ErrorNorms& represented,
        const LiftedErrorNorm::Result& lifted, Real patchTolerance)
      {
        expectDecomposition(lifted);
        for (const auto& error : {represented, lifted.field})
        {
          EXPECT_TRUE(error.isFinite());
          EXPECT_LT(error.getL2(), patchTolerance);
          EXPECT_LT(error.getH1Seminorm(), patchTolerance);
        }
        EXPECT_NEAR(lifted.total.getL2(), lifted.geometry.getL2(), patchTolerance);
        EXPECT_NEAR(
          lifted.total.getH1Seminorm(), lifted.geometry.getH1Seminorm(), patchTolerance);
      }

      void appendRepresentable(Real h, const ErrorNorms& represented,
        const LiftedErrorNorm::Result& lifted, Real patchTolerance)
      {
        expectRepresentable(represented, lifted, patchTolerance);
        const auto errors = components(represented, lifted);
        for (size_t component : {2u, 3u})
        {
          ASSERT_TRUE(errors[component].isFinite());
          ASSERT_GT(errors[component].getL2(), 0);
          ASSERT_GT(errors[component].getH1Seminorm(), 0);
          m_histories[component].append(h, errors[component]);
        }
      }

      void expectGeometryRates(size_t geometryDegree) const
      {
        expectRatesFrom(2, geometryDegree, geometryDegree);
      }

      static void expectGeometrySensitivity(
        const LiftedErrorNorm::Result& base, const LiftedErrorNorm::Result& refined)
      {
        for (const auto& pair :
          {std::pair{base.geometry.getL2(), refined.geometry.getL2()},
            std::pair{base.geometry.getH1Seminorm(), refined.geometry.getH1Seminorm()},
            std::pair{base.total.getL2(), refined.total.getL2()},
            std::pair{base.total.getH1Seminorm(), refined.total.getH1Seminorm()}})
        {
          ASSERT_GT(pair.first, 0);
          ASSERT_TRUE(std::isfinite(pair.second));
          EXPECT_LT(std::abs(pair.second / pair.first - 1), SensitivityTolerance);
        }
      }
      static void expectSensitivity(const Components& base, const Components& refined)
      {
        for (size_t component = 0; component < base.size(); ++component)
        {
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
      }

    private:
      void expectRatesFrom(size_t firstComponent, size_t fieldDegree,
        size_t geometryDegree, size_t lastComponent = 4) const
      {
        for (size_t component = firstComponent; component < lastComponent; ++component)
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

      std::array<ErrorHistory, 4> m_histories;
  };
}
#endif
