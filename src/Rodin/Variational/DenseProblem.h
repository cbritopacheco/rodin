/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DENSEPROBLEM_H
#define RODIN_VARIATIONAL_DENSEPROBLEM_H

#include "Problem.h"

namespace Rodin::Variational
{
  template <class ... Parameters>
  class DenseProblem;

  template <class TrialFES, class TestFES>
  DenseProblem(TrialFunction<TrialFES>&, TestFunction<TestFES>&)
    -> DenseProblem<TrialFES, TestFES,
          Math::Matrix<
            typename FormLanguage::Mult<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
          Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>;

  /**
   * @defgroup DenseProblemSpecializations DenseProblem Template Specializations
   * @brief Template specializations of the DenseProblem class.
   * @see DenseProblem
   */

  /**
   * @ingroup DenseProblemSpecializations
   * @brief General class to assemble linear systems with `Math::Matrix`
   * and `Math::Vector` types in a serial context.
   */
  template <class TrialFES, class TestFES>
  class DenseProblem<
    TrialFES, TestFES,
    Math::Matrix<
      typename FormLanguage::Mult<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
    Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>
    : public Problem<
        TrialFES, TestFES,
        Math::Matrix<
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TrialFES>::ScalarType>
          ::Type>,
        Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>
    {
      public:
        using Parent = Problem<
          TrialFES, TestFES,
          Math::Matrix<
            typename FormLanguage::Mult<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TrialFES>::ScalarType>
            ::Type>,
          Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>;

        using Parent::Parent;
        using Parent::operator=;
    };
}

#endif

