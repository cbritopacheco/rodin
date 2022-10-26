/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_AND_H
#define RODIN_VARIATIONAL_AND_H

#include "ForwardDecls.h"
#include "Exceptions.h"
#include "GridFunction.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
   /**
    * @defgroup ANDSpecializations AND Template Specializations
    * @brief Template specializations of the AND class.
    * @see AND
    */

   /**
    * @ingroup ANDSpecializations
    * @brief Logical AND operator between two instances of BooleanFunctionBase
    */
   template <>
   class AND<BooleanFunctionBase, BooleanFunctionBase> : public BooleanFunctionBase
   {
      public:
         AND(const BooleanFunctionBase& lhs, const BooleanFunctionBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {
            if (lhs.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, lhs.getRangeType()).raise();
            if (rhs.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, rhs.getRangeType()).raise();
         }

         AND(const AND& other)
            : BooleanFunctionBase(other),
              m_lhs(other.m_lhs->copy()),
              m_rhs(other.m_rhs->copy())
         {}

         AND(AND&& other)
            : BooleanFunctionBase(std::move(other)),
              m_lhs(std::move(other.m_lhs)),
              m_rhs(std::move(other.m_rhs))
         {}

         bool getValue(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override
         {
            return m_lhs->getValue(trans, ip) && m_rhs->getValue(trans, ip);
         }

         AND* copy() const noexcept override
         {
            return new AND(*this);
         }

      private:
         std::unique_ptr<BooleanFunctionBase> m_lhs;
         std::unique_ptr<BooleanFunctionBase> m_rhs;
   };
   AND(const BooleanFunctionBase&, const BooleanFunctionBase&)
      -> AND<BooleanFunctionBase, BooleanFunctionBase>;

   AND<BooleanFunctionBase, BooleanFunctionBase>
   operator&&(const BooleanFunctionBase&, const BooleanFunctionBase&);

   AND<BooleanFunctionBase, BooleanFunctionBase>
   operator&&(bool lhs, const BooleanFunctionBase& rhs);

   AND<BooleanFunctionBase, BooleanFunctionBase>
   operator&&(const BooleanFunctionBase& lhs, bool rhs);
}

#endif


