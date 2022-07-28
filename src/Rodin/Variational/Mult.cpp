/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Exceptions.h"

#include "Mult.h"

namespace Rodin::Variational
{
   // ---- FunctionBase, FunctionBase ----------------------------------------
   Mult<FunctionBase, FunctionBase>::Mult(const FunctionBase& lhs, const FunctionBase& rhs)
      : m_lhs(lhs.copy()), m_rhs(rhs.copy())
   {
      if (lhs.getRangeType() != RangeType::Scalar && rhs.getRangeType() != RangeType::Scalar)
      {
         // One of the fields must be scalar valued!
         IncompatibleShapeException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
      }
   }

   Mult<FunctionBase, FunctionBase>::Mult(const Mult& other)
      :  FunctionBase(other),
         m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
   {}

   Mult<FunctionBase, FunctionBase>::Mult(Mult&& other)
      :  FunctionBase(std::move(other)),
         m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
   {}

   RangeShape Mult<FunctionBase, FunctionBase>::getRangeShape() const
   {
      assert(
         m_lhs->getRangeType() == RangeType::Scalar ||
         m_rhs->getRangeType() == RangeType::Scalar);

      if (m_lhs->getRangeType() == RangeType::Scalar)
         return m_rhs->getRangeShape();
      else
         return m_lhs->getRangeShape();
   }

   Mult<FunctionBase, FunctionBase>&
   Mult<FunctionBase, FunctionBase>::traceOf(const std::set<int>& attrs)
   {
      FunctionBase::traceOf(attrs);
      m_lhs->traceOf(attrs);
      m_rhs->traceOf(attrs);
      return *this;
   }

   void Mult<FunctionBase, FunctionBase>::getValue(
         mfem::DenseMatrix& value,
         mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip) const
   {
      double s;
      if (m_lhs->getRangeType() == RangeType::Scalar)
      {
         mfem::DenseMatrix tmp;
         m_lhs->getValue(tmp, trans, ip);
         s = tmp(0, 0);
         m_rhs->getValue(value, trans, ip);
      }
      else
      {
         mfem::DenseMatrix tmp;
         m_rhs->getValue(tmp, trans, ip);
         s = tmp(0, 0);
         m_lhs->getValue(value, trans, ip);
      }
      value *= s;
   }

   Mult<FunctionBase, FunctionBase>
   operator*(const FunctionBase& lhs, const FunctionBase& rhs)
   {
      return Mult(lhs, rhs);
   }
}
