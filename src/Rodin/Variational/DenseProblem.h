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
  template <class LinearSystem, class U, class V>
  class DenseProblem<LinearSystem, U, V> : public Problem<LinearSystem, U, V>
  {
    public:
      using Parent = Problem<LinearSystem, U, V>;
      using Parent::Parent;
      using Parent::operator=;
  };

  template <class U, class V>
  DenseProblem(U& u, V& v)
    -> DenseProblem<
        Math::LinearSystem<
          Math::Matrix<
            typename FormLanguage::Traits<typename FormLanguage::Traits<U>::FESType>
            ::ScalarType>,
          Math::Vector<
            typename FormLanguage::Traits<typename FormLanguage::Traits<V>::FESType>
            ::ScalarType>>,
          U, V>;
}

#endif

