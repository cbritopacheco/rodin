/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "ScalarCoefficient.h"

#include "Restriction.h"

namespace Rodin::Variational
{
   Restriction<ScalarCoefficientBase>::Restriction(const ScalarCoefficientBase& s)
      : m_s(s.copy())
   {}

   Restriction<ScalarCoefficientBase>::Restriction(const Restriction& other)
      :  m_attr(other.m_attr),
         m_s(other.m_s->copy())
   {}

   Restriction<ScalarCoefficientBase>&
   Restriction<ScalarCoefficientBase>::to(const std::set<int>& attr)
   {
      m_attr = attr;
      return *this;
   }

   const std::set<int>& Restriction<ScalarCoefficientBase>::getAttributes() const
   {
      assert(m_attr.size() > 0);
      return m_attr;
   }

   const ScalarCoefficientBase&
   Restriction<ScalarCoefficientBase>::getScalarCoefficient() const
   {
      assert(m_s);
      return *m_s;
   }
}
