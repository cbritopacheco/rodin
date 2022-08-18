/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BOOLEANFUNCTION_H
#define RODIN_VARIATIONAL_BOOLEANFUNCTION_H

#include "ForwardDecls.h"
#include "Function.h"
#include "RangeShape.h"

namespace Rodin::Variational
{
   class BooleanFunctionBase : public FunctionBase
   {
      public:
         BooleanFunctionBase() = default;

         BooleanFunctionBase(const BooleanFunctionBase& other)
            : FunctionBase(other)
         {}

         BooleanFunctionBase(BooleanFunctionBase&& other)
            : FunctionBase(std::move(other))
         {}

         virtual ~BooleanFunctionBase() = default;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override
         {
            value.SetSize(1, 1);
            value(0, 0) = getValue(trans, ip);
         }

         RangeShape getRangeShape() const override
         {
            return {1, 1};
         }

         RangeType getRangeType() const override
         {
            return RangeType::Scalar;
         }

         virtual bool getValue(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const = 0;

         virtual BooleanFunctionBase* copy() const noexcept override = 0;
   };

   template <>
   class BooleanFunction<bool> : public BooleanFunctionBase
   {
      public:
         BooleanFunction(bool v)
            : m_v(v)
         {}

         BooleanFunction(const BooleanFunction& other)
            : BooleanFunctionBase(other),
              m_v(other.m_v)
         {}

         BooleanFunction(BooleanFunction&& other)
            : BooleanFunctionBase(std::move(other)),
              m_v(other.m_v)
         {}

         bool getValue(
               mfem::ElementTransformation&,
               const mfem::IntegrationPoint&) const override
         {
            return m_v;
         }

         BooleanFunction* copy() const noexcept override
         {
            return new BooleanFunction(*this);
         }

      private:
         bool m_v;
   };
   BooleanFunction(bool) -> BooleanFunction<bool>;
}

#endif

