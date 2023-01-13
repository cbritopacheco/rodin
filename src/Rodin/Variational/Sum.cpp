/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "RangeShape.h"
#include "Exceptions.h"

#include "BilinearFormIntegrator.h"

#include "Sum.h"

namespace Rodin::Variational
{
   Sum<FunctionBase, FunctionBase>::Sum(const FunctionBase& lhs, const FunctionBase& rhs)
      : m_lhs(lhs.copy()), m_rhs(rhs.copy())
   {
      if (lhs.getRangeShape() != rhs.getRangeShape())
         RangeShapeMismatchException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
   }

   Sum<FunctionBase, FunctionBase>::Sum(const Sum& other)
      :  FunctionBase(other),
         m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
   {}

   Sum<FunctionBase, FunctionBase>::Sum(Sum&& other)
      :  FunctionBase(std::move(other)),
         m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
   {}

   Sum<FunctionBase, FunctionBase>&
   Sum<FunctionBase, FunctionBase>::traceOf(Geometry::Attribute attrs)
   {
      FunctionBase::traceOf(attrs);
      m_lhs->traceOf(attrs);
      m_rhs->traceOf(attrs);
      return *this;
   }

   RangeShape Sum<FunctionBase, FunctionBase>::getRangeShape() const
   {
      assert(m_lhs->getRangeShape() == m_rhs->getRangeShape());
      return m_lhs->getRangeShape();
   }

   Sum<FunctionBase, FunctionBase>&
   Sum<FunctionBase, FunctionBase>::operator+=(const FunctionBase& lhs)
   {
      auto sum = new Sum(lhs, *m_rhs);
      m_rhs.reset(sum);
      return *this;
   }

   Sum<FunctionBase, FunctionBase>
   operator+(const FunctionBase& lhs, const FunctionBase& rhs)
   {
      return Sum(lhs, rhs);
   }

   FormLanguage::List<BilinearFormIntegratorBase>
   operator+(
         const BilinearFormIntegratorBase& lhs,
         const BilinearFormIntegratorBase& rhs)
   {
      FormLanguage::List<BilinearFormIntegratorBase> res;
      res.add(lhs);
      res.add(rhs);
      return res;
   }

   FormLanguage::List<BilinearFormIntegratorBase>
   operator+(
         const BilinearFormIntegratorBase& lhs,
         const FormLanguage::List<BilinearFormIntegratorBase>& rhs)
   {
      FormLanguage::List<BilinearFormIntegratorBase> res(rhs);
      res.add(lhs);
      return res;
   }

   FormLanguage::List<BilinearFormIntegratorBase>
   operator+(
         const FormLanguage::List<BilinearFormIntegratorBase>& lhs,
         const BilinearFormIntegratorBase& rhs)
   {
      FormLanguage::List<BilinearFormIntegratorBase> res(lhs);
      res.add(rhs);
      return res;
   }

   FormLanguage::List<BilinearFormIntegratorBase>
   operator+(
         const FormLanguage::List<BilinearFormIntegratorBase>& lhs,
         const FormLanguage::List<BilinearFormIntegratorBase>& rhs)
   {
      FormLanguage::List<BilinearFormIntegratorBase> res(lhs);
      for (const auto& p : rhs)
         res.add(p);
      return res;
   }
}
