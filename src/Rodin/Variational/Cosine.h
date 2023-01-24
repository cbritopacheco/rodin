/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_COSINE_H
#define RODIN_VARIATIONAL_COSINE_H

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"
#include "Exceptions.h"

namespace Rodin::Variational
{
  /**
   * @defgroup CosineSpecializations Cosine Template Specializations
   * @brief Template specializations of the Cosine class.
   * @see Cosine
   */

  /**
   * @ingroup CosSpecializations
   */
  template <>
  class Cosine<FunctionBase> : public ScalarFunctionBase
  {
    public:
      using Operand = FunctionBase;

      Cosine(const FunctionBase& v)
        : m_v(v.copy())
      {
        if (v.getRangeType() != RangeType::Scalar)
          throw UnexpectedRangeTypeException(RangeType::Scalar, v.getRangeType());
      }

      Cosine(const Cosine& other)
        :  ScalarFunctionBase(other),
          m_v(other.m_v->copy())
      {}

      Cosine(Cosine&& other)
        :  ScalarFunctionBase(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      Cosine& traceOf(Geometry::Attribute attrs) override
      {
        ScalarFunctionBase::traceOf(attrs);
        m_v->traceOf(attrs);
        return *this;
      }

      FunctionValue getValue(const Geometry::Point& p) const override
      {
        return Math::cos(m_v->getValue(p).scalar());
      }

      Cosine* copy() const noexcept override
      {
        return new Cosine(*this);
      }

    private:
      std::unique_ptr<FunctionBase> m_v;
  };
  Cosine(const FunctionBase&) -> Cosine<FunctionBase>;

  /**
   * @brief Convenience function to create objects of type
   * Cosine<FunctionBase>.
   * @param[in] op Scalar function type
   */
  Cosine<FunctionBase> cos(const FunctionBase& op);
}

#endif
