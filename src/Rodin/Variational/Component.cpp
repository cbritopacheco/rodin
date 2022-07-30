/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Exceptions.h"

#include "Component.h"

namespace Rodin::Variational
{
   Component<FunctionBase>::Component(const FunctionBase& v, int component)
      :  m_v(v.copy()),
         m_idx(component)
   {
      if (v.getRangeType() != RangeType::Vector)
         UnexpectedRangeTypeException(RangeType::Vector, v.getRangeType()).raise();
   }

   Component<FunctionBase>::Component(const Component& other)
      :  ScalarFunctionBase(other),
         m_v(other.m_v->copy()),
         m_idx(other.m_idx)
   {}

   Component<FunctionBase>::Component(Component&& other)
      :  ScalarFunctionBase(std::move(other)),
         m_v(std::move(other.m_v)),
         m_idx(other.m_idx)
   {}

   int Component<FunctionBase>::getIndex() const
   {
      return m_idx;
   }

   Component<FunctionBase>& Component<FunctionBase>::traceOf(const std::set<int>& attrs)
   {
      ScalarFunctionBase::traceOf(attrs);
      m_v->traceOf(attrs);
      return *this;
   }

   double Component<FunctionBase>::getValue(
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const
   {
      mfem::DenseMatrix v;
      m_v->getValue(v, trans, ip);
      assert(m_idx < v.NumRows());
      return v(m_idx, 0);
   }

   Component<FunctionBase>* Component<FunctionBase>::copy() const noexcept
   {
      return new Component(*this);
   }
}

