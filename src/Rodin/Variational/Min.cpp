/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Exceptions.h"
#include "RangeShape.h"

#include "Min.h"

namespace Rodin::Variational
{
   Min::Min(const FunctionBase& a, double b)
      : Min(a, ScalarFunction(b))
   {}

   Min::Min(double a, const FunctionBase& b)
      : Min(ScalarFunction(a), b)
   {}

   /**
    * @bref Constructs the power object
    * @param[in] s Base value
    * @param[in] p Power
    */
   Min::Min(const FunctionBase& a, const FunctionBase& b)
      : m_a(a.copy()),
        m_b(b.copy())
   {
      if (a.getRangeType() != RangeType::Scalar)
         UnexpectedRangeTypeException(RangeType::Scalar, a.getRangeType());
      if (b.getRangeType() != RangeType::Scalar)
         UnexpectedRangeTypeException(RangeType::Scalar, b.getRangeType());
   }

   Min::Min(const Min& other)
      : ScalarFunctionBase(other),
        m_a(other.m_a->copy()),
        m_b(other.m_b->copy())
   {}

   Min::Min(Min&& other)
      : ScalarFunctionBase(std::move(other)),
        m_a(std::move(other.m_a)),
        m_b(std::move(other.m_b))
   {}

   Min& Min::traceOf(const std::set<int>& attrs)
   {
      ScalarFunctionBase::traceOf(attrs);
      m_a->traceOf(attrs);
      m_b->traceOf(attrs);
      return *this;
   }

   Min* Min::copy() const noexcept
   {
      return new Min(*this);
   }
}
