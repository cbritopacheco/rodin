/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_POW_H
#define RODIN_VARIATIONAL_POW_H

#include <cmath>

#include "RangeShape.h"
#include "Exceptions.h"
#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup PowSpecializations Pow Template Specializations
   * @brief Template specializations of the Pow class.
   * @see Pow
   */

  /**
   * @ingroup PowSpecializations
   * @brief Represent the power function.
   *
   * This class represents the function @f$ f : \mathbb{R} \rightarrow
   * \mathbb{R} @f$ defined by the exponentiation of a base value @f$ x \in
   * \mathbb{R} @f$ to the power @f$ p @f$:
   * @f$
   *  f(x) = x^p \ .
   * @f$
   */
  template <class Number>
  class Pow<FunctionBase, Number> : public ScalarFunctionBase
  {
    static_assert(std::is_arithmetic_v<Number>, "T must be an arithmetic type");
    public:
      /**
       * @brief Constructs the power object
       * @param[in] s Base value
       * @param[in] p Power
       */
      Pow(const FunctionBase& s, Number p)
        : m_s(s.copy()),
          m_p(p)
      {
        if (s.getRangeType() != RangeType::Scalar)
          UnexpectedRangeTypeException(RangeType::Scalar, s.getRangeType());
      }

      Pow(const Pow& other)
        : ScalarFunctionBase(other),
          m_s(other.m_s->copy()),
          m_p(other.m_p)
      {}

      Pow(Pow&& other)
        : ScalarFunctionBase(std::move(other)),
          m_s(std::move(other.m_s)),
          m_p(other.m_p)
      {}

      Pow& traceOf(Geometry::Attribute attrs) override
      {
        ScalarFunctionBase::traceOf(attrs);
        m_s->traceOf(attrs);
        return *this;
      }

      FunctionValue getValue(const Geometry::Point& p) const override
      {
        return std::pow(m_s->getValue(p).scalar(), m_p);
      }

      Pow* copy() const noexcept override
      {
        return new Pow(*this);
      }
    private:
      std::unique_ptr<FunctionBase> m_s;
      Number m_p;
  };
  template <class Number>
  Pow(const FunctionBase&, Number) -> Pow<FunctionBase, Number>;
}

#endif
