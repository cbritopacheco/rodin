/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GT_H
#define RODIN_VARIATIONAL_GT_H

#include "ForwardDecls.h"
#include "Exceptions.h"
#include "GridFunction.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup GTSpecializations GT Template Specializations
   * @brief Template specializations of the GT class.
   * @see GT
   */

  template <>
  class GT<FunctionBase, FunctionBase> : public BooleanFunctionBase
  {
    public:
      GT(const FunctionBase& lhs, const FunctionBase& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        if (lhs.getRangeType() != RangeType::Scalar)
          UnexpectedRangeTypeException(RangeType::Scalar, lhs.getRangeType()).raise();
        if (rhs.getRangeType() != RangeType::Scalar)
          UnexpectedRangeTypeException(RangeType::Scalar, rhs.getRangeType()).raise();
      }

      GT(const GT& other)
        : BooleanFunctionBase(other),
          m_lhs(other.m_lhs->copy()),
          m_rhs(other.m_rhs->copy())
      {}

      GT(GT&& other)
        : BooleanFunctionBase(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      FunctionValue getValue(const Geometry::Point& p) const override
      {
        return m_lhs->getValue(p).boolean() > m_rhs->getValue(p).boolean();
      }

      GT* copy() const noexcept override
      {
        return new GT(*this);
      }

    private:
      std::unique_ptr<FunctionBase> m_lhs;
      std::unique_ptr<FunctionBase> m_rhs;
  };
  GT(const FunctionBase&, const FunctionBase&) -> GT<FunctionBase, FunctionBase>;

  GT<FunctionBase, FunctionBase>
  operator>(const FunctionBase&, const FunctionBase&);

  template <class T>
  std::enable_if_t<std::is_arithmetic_v<T>, GT<FunctionBase, FunctionBase>>
  operator>(T lhs, const FunctionBase& rhs)
  {
    return GT(ScalarFunction(lhs), rhs);
  }

  template <class T>
  std::enable_if_t<std::is_arithmetic_v<T>, GT<FunctionBase, FunctionBase>>
  operator>(const FunctionBase& lhs, T rhs)
  {
    return GT(lhs, ScalarFunction(rhs));
  }
}

#endif

