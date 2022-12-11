/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "RangeShape.h"
#include "Exceptions.h"
#include "Division.h"

namespace Rodin::Variational
{
   Division<FunctionBase, FunctionBase>::Division(
         const FunctionBase& lhs, const FunctionBase& rhs)
      : m_lhs(lhs.copy()), m_rhs(rhs.copy())
   {
      if (rhs.getRangeType() != RangeType::Scalar)
      {
         // Right hand side must be a scalar!
         IncompatibleShapeException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
      }
   }

   Division<FunctionBase, FunctionBase>::Division(const Division& other)
      :  FunctionBase(other),
         m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
   {}

   Division<FunctionBase, FunctionBase>::Division(Division&& other)
      :  FunctionBase(std::move(other)),
         m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
   {}

   Division<FunctionBase, FunctionBase>&
   Division<FunctionBase, FunctionBase>::traceOf(const std::set<int>& attrs)
   {
      FunctionBase::traceOf(attrs);
      m_lhs->traceOf(attrs);
      m_rhs->traceOf(attrs);
      return *this;
   }

   RangeShape Division<FunctionBase, FunctionBase>::getRangeShape() const
   {
      return m_lhs->getRangeShape();
   }

   Division<FunctionBase, FunctionBase>
   operator/(const FunctionBase& lhs, const FunctionBase& rhs)
   {
      return Division(lhs, rhs);
   }
}

