/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_NAMEDFORMSTORAGE_H
#define RODIN_VARIATIONAL_NAMEDFORMSTORAGE_H
/**
 * @file
 * @brief Operator selection and ownership for named bilinear forms.
 */
#include "Rodin/Math/SparseMatrix.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Selects named-form operators from the trial solution type.
   * @tparam Solution Trial solution type.
   * @tparam Scalar Operator entry type.
   */
  template <class Solution, class Scalar>
  struct NamedFormOperatorType
  {
      /// @brief Local sparse operator type.
      using Type = Math::SparseMatrix<Scalar>;
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Owns the operator of a named bilinear form.
   * @tparam Operator Matrix storage type.
   */
  template <class Operator>
  class NamedFormStorage
  {
    public:
      /**
       * @brief Gets the writable operator.
       * @returns Reference to the owned matrix.
       */
      Operator& get()
      {
        return m_operator;
      }
      /**
       * @brief Inspects the operator.
       * @returns Const reference to the owned matrix.
       */
      const Operator& get() const
      {
        return m_operator;
      }

    private:
      Operator m_operator{};
  };
}
#endif
