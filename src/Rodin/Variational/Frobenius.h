/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FROBENIUS_H
#define RODIN_VARIATIONAL_FROBENIUS_H

#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup FrobeniusSpecializations Frobenius Template Specializations
   * @brief Template specializations of the Frobenius class.
   * @see Frobenius
   */

  /**
   * @ingroup FrobeniusSpecializations
   */
  template <>
  class Frobenius<FunctionBase> : public ScalarFunctionBase
  {
    public:
      using Operand = FunctionBase;

      Frobenius(const FunctionBase& v)
        : m_v(v.copy())
      {}

      Frobenius(const Frobenius& other)
        :  ScalarFunctionBase(other),
          m_v(other.m_v->copy())
      {}

      Frobenius(Frobenius&& other)
        :  ScalarFunctionBase(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      Frobenius& traceOf(Geometry::Attribute attrs) override
      {
        ScalarFunctionBase::traceOf(attrs);
        m_v->traceOf(attrs);
        return *this;
      }

      FunctionValue getValue(const Geometry::Point& p) const override
      {
        switch (m_v->getRangeType())
        {
          case RangeType::Scalar:
          {
            return std::abs(m_v->getValue(p).scalar());
          }
          case RangeType::Vector:
          {
            return m_v->getValue(p).vector().norm();
          }
          case RangeType::Matrix:
          {
            return m_v->getValue(p).matrix().norm();
          }
          default:
          {
            assert(false);
            return 0.0;
          }
        }
      }

      Frobenius* copy() const noexcept override
      {
        return new Frobenius(*this);
      }

    private:
      std::unique_ptr<FunctionBase> m_v;
  };
  Frobenius(const FunctionBase&) -> Frobenius<FunctionBase>;
}

#endif

