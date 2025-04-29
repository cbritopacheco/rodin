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
  class SparseProblem;

  template <class TrialFES, class TestFES>
  SparseProblem(TrialFunction<TrialFES>&, TestFunction<TestFES>&)
    -> SparseProblem<TrialFES, TestFES,
          Math::Matrix<
            typename FormLanguage::Mult<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
          Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>;

  /**
   * @defgroup SparseProblemSpecializations SparseProblem Template Specializations
   * @brief Template specializations of the SparseProblem class.
   * @see SparseProblem
   */

  /**
   * @ingroup SparseProblemSpecializations
   * @brief General class to assemble linear systems with `Math::SparseMatrix`
   * and `Math::Vector` types in a serial context.
   */
  template <class TrialFES, class TestFES>
  class SparseProblem<
    TrialFES, TestFES,
    Math::Matrix<
      typename FormLanguage::Mult<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
    Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>
    : public Problem<
        TrialFES, TestFES,
        Math::SparseMatrix<
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TrialFES>::ScalarType>
          ::Type>,
        Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>
    {
      public:
        using Parent = Problem<
          TrialFES, TestFES,
          Math::SparseMatrix<
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

