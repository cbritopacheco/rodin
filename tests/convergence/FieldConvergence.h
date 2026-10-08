/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_FIELDCONVERGENCE_H
#define RODIN_TESTS_CONVERGENCE_FIELDCONVERGENCE_H

#include <array>
#include <gtest/gtest.h>

#include "Convergence.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Independent L2/H1 acceptance for a sequence of conforming fields.
   * @par Architecture
   * Each sample has one refinement parameter and one error pair per field.
   * ErrorHistory owns the samples and computes the rates. This study enforces
   * at least three levels, positive finite errors, strict reduction and the
   * prescribed floors for every field and adjacent interval. Exact patches
   * and discontinuous L2-only studies have different contracts and are not
   * admitted through this positive-error interface.
   */
  template <size_t Fields>
  class FieldConvergence
  {
      static_assert(Fields > 0);

    public:
      FieldConvergence& append(
        Real parameter, const std::array<ErrorNorms, Fields>& errors)
      {
        for (size_t field = 0; field < Fields; ++field)
          m_histories[field].append(parameter, errors[field]);
        return *this;
      }

      size_t getSize() const
      {
        return m_histories[0].getSize();
      }

      void expectExponentialFloor(const Rates& floor) const
      {
        expectFloor<true>(floor);
      }

      void expectAlgebraicFloor(const Rates& floor) const
      {
        expectFloor<false>(floor);
      }

    private:
      template <bool Exponential>
      void expectFloor(const Rates& floor) const
      {
        ASSERT_GE(m_histories[0].getSize(), 3u);
        for (size_t field = 0; field < Fields; ++field)
        {
          for (size_t i = 1; i < m_histories[field].getSize(); ++i)
          {
            SCOPED_TRACE(::testing::Message() << "field=" << field << " interval=" << i);
            const auto& coarse = m_histories[field].getSample(i - 1).error;
            const auto& fine = m_histories[field].getSample(i).error;
            ASSERT_TRUE(coarse.isFinite());
            ASSERT_TRUE(fine.isFinite());
            ASSERT_GT(coarse.getL2(), 0);
            ASSERT_GT(coarse.getH1Seminorm(), 0);
            ASSERT_GT(fine.getL2(), 0);
            ASSERT_GT(fine.getH1Seminorm(), 0);
            const Real coarseParameter = m_histories[field].getSample(i - 1).parameter;
            const Real fineParameter = m_histories[field].getSample(i).parameter;
            ASSERT_TRUE(std::isfinite(coarseParameter));
            ASSERT_TRUE(std::isfinite(fineParameter));
            if constexpr (Exponential)
            {
              ASSERT_GT(fineParameter, coarseParameter);
            }
            else
            {
              ASSERT_GT(fineParameter, 0);
              ASSERT_GT(coarseParameter, fineParameter);
            }
            const auto rate = Exponential ? m_histories[field].getExponentialRates(i)
                                          : m_histories[field].getAlgebraicRates(i);
            SCOPED_TRACE(::testing::Message()
              << "L2=" << coarse.getL2() << " -> " << fine.getL2()
              << " H1=" << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm()
              << " rates=" << rate.getL2() << "," << rate.getH1Seminorm());
            EXPECT_GT(coarse.getL2(), fine.getL2());
            EXPECT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
            EXPECT_GT(rate.getL2(), floor.getL2());
            EXPECT_GT(rate.getH1Seminorm(), floor.getH1Seminorm());
          }
        }
      }

      std::array<ErrorHistory, Fields> m_histories;
  };
}

#endif
