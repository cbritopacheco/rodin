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

         Frobenius& traceOf(const std::set<int>& attrs) override
         {
            ScalarFunctionBase::traceOf(attrs);
            m_v->traceOf(attrs);
            return *this;
         }

         FunctionValue getValue(const Geometry::Point& p) const override
         {
            return m_v->getValue(p).matrix().FNorm();
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

