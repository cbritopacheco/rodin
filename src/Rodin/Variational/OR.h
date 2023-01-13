/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_OR_H
#define RODIN_VARIATIONAL_OR_H

#include "ForwardDecls.h"
#include "Exceptions.h"
#include "GridFunction.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
   /**
    * @defgroup ORSpecializations OR Template Specializations
    * @brief Template specializations of the OR class.
    * @see OR
    */

   template <>
   class OR<BooleanFunctionBase, BooleanFunctionBase> : public BooleanFunctionBase
   {
      public:
         OR(const BooleanFunctionBase& lhs, const BooleanFunctionBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {
            if (lhs.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, lhs.getRangeType()).raise();
            if (rhs.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, rhs.getRangeType()).raise();
         }

         OR(const OR& other)
            : BooleanFunctionBase(other),
              m_lhs(other.m_lhs->copy()),
              m_rhs(other.m_rhs->copy())
         {}

         OR(OR&& other)
            : BooleanFunctionBase(std::move(other)),
              m_lhs(std::move(other.m_lhs)),
              m_rhs(std::move(other.m_rhs))
         {}

         FunctionValue getValue(const Geometry::Point& p) const override
         {
            return m_lhs->getValue(p).boolean() || m_rhs->getValue(p).boolean();
         }

         OR* copy() const noexcept override
         {
            return new OR(*this);
         }

      private:
         std::unique_ptr<BooleanFunctionBase> m_lhs;
         std::unique_ptr<BooleanFunctionBase> m_rhs;
   };
   OR(const BooleanFunctionBase&, const BooleanFunctionBase&)
      -> OR<BooleanFunctionBase, BooleanFunctionBase>;

   OR<BooleanFunctionBase, BooleanFunctionBase>
   operator||(const BooleanFunctionBase&, const BooleanFunctionBase&);

   OR<BooleanFunctionBase, BooleanFunctionBase>
   operator||(bool lhs, const BooleanFunctionBase& rhs);

   OR<BooleanFunctionBase, BooleanFunctionBase>
   operator||(const BooleanFunctionBase& lhs, bool rhs);
}

#endif



