/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "RangeShape.h"
#include "Exceptions.h"

#include "Max.h"

namespace Rodin::Variational
{
   Max::Max(const FunctionBase& a, double b)
      : Max(a, ScalarFunction(b))
   {}

   Max::Max(double a, const FunctionBase& b)
      : Max(ScalarFunction(a), b)
   {}

   /**
    * @bref Constructs the power object
    * @param[in] s Base value
    * @param[in] p Power
    */
   Max::Max(const FunctionBase& a, const FunctionBase& b)
      : m_a(a.copy()),
        m_b(b.copy())
   {
      if (a.getRangeType() != RangeType::Scalar)
         UnexpectedRangeTypeException(RangeType::Scalar, a.getRangeType());
      if (b.getRangeType() != RangeType::Scalar)
         UnexpectedRangeTypeException(RangeType::Scalar, b.getRangeType());
   }

   Max::Max(const Max& other)
      : ScalarFunctionBase(other),
        m_a(other.m_a->copy()),
        m_b(other.m_b->copy())
   {}

   Max::Max(Max&& other)
      : ScalarFunctionBase(std::move(other)),
        m_a(std::move(other.m_a)),
        m_b(std::move(other.m_b))
   {}

   Max& Max::traceOf(const std::set<int>& attrs)
   {
      ScalarFunctionBase::traceOf(attrs);
      m_a->traceOf(attrs);
      m_b->traceOf(attrs);
      return *this;
   }

   Max* Max::copy() const noexcept
   {
      return new Max(*this);
   }
}
