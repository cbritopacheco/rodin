/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "ScalarFunction.h"

#include "Restriction.h"

namespace Rodin::Variational
{
   Restriction<ScalarFunctionBase>::Restriction(const ScalarFunctionBase& s)
      : m_s(s.copy())
   {}

   Restriction<ScalarFunctionBase>::Restriction(const Restriction& other)
      :  m_attr(other.m_attr),
         m_s(other.m_s->copy())
   {}

   Restriction<ScalarFunctionBase>&
   Restriction<ScalarFunctionBase>::to(const std::set<int>& attr)
   {
      m_attr = attr;
      return *this;
   }

   const std::set<int>& Restriction<ScalarFunctionBase>::getAttributes() const
   {
      assert(m_attr.size() > 0);
      return m_attr;
   }

   const ScalarFunctionBase&
   Restriction<ScalarFunctionBase>::getScalarFunction() const
   {
      assert(m_s);
      return *m_s;
   }
}
