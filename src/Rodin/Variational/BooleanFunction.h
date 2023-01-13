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
   /**
    * @defgroup BooleanFunctionSpecializations BooleanFunction Template Specializations
    * @brief Template specializations of the BooleanFunction class.
    * @see BooleanFunction
    */

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

         FunctionValue::Boolean operator()(const Geometry::Point& p) const
         {
            return getValue(p);
         }

         RangeShape getRangeShape() const override
         {
            return {1, 1};
         }

         RangeType getRangeType() const override
         {
            return RangeType::Scalar;
         }

         virtual BooleanFunctionBase* copy() const noexcept override = 0;
   };

   /**
    * @ingroup BooleanFunctionSpecializations
    */
   template <>
   class BooleanFunction<FunctionValue::Boolean> : public BooleanFunctionBase
   {
      public:
         BooleanFunction(FunctionValue::Boolean v)
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

         FunctionValue getValue(const Geometry::Point& p) const override
         {
            return m_v;
         }

         BooleanFunction* copy() const noexcept override
         {
            return new BooleanFunction(*this);
         }

      private:
         FunctionValue::Boolean m_v;
   };
   BooleanFunction(FunctionValue::Boolean) -> BooleanFunction<FunctionValue::Boolean>;
}

#endif

