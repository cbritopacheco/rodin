/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TANGENT_H
#define RODIN_VARIATIONAL_TANGENT_H

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"
#include "Exceptions.h"

namespace Rodin::Variational
{
   /**
    * @defgroup TangentSpecializations Tangent Template Specializations
    * @brief Template specializations of the Tangent class.
    * @see Tangent
    */

   /**
    * @ingroup TangentSpecializations
    */
   template <>
   class Tangent<FunctionBase> : public ScalarFunctionBase
   {
      public:
         using Operand = FunctionBase;

         Tangent(const FunctionBase& v)
            : m_v(v.copy())
         {
            if (v.getRangeType() != RangeType::Scalar)
               throw UnexpectedRangeTypeException(RangeType::Scalar, v.getRangeType());
         }

         Tangent(const Tangent& other)
            :  ScalarFunctionBase(other),
               m_v(other.m_v->copy())
         {}

         Tangent(Tangent&& other)
            :  ScalarFunctionBase(std::move(other)),
               m_v(std::move(other.m_v))
         {}

         Tangent& traceOf(const std::set<int>& attrs) override
         {
            ScalarFunctionBase::traceOf(attrs);
            m_v->traceOf(attrs);
            return *this;
         }

         double getValue(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override
         {
            mfem::DenseMatrix s;
            m_v->getValue(s, trans, ip);
            return Math::tan(s(0, 0));
         }

         Tangent* copy() const noexcept override
         {
            return new Tangent(*this);
         }

      private:
         std::unique_ptr<FunctionBase> m_v;
   };
   Tangent(const FunctionBase&) -> Tangent<FunctionBase>;

   Tangent<FunctionBase> tan(const FunctionBase&);

}

#endif

