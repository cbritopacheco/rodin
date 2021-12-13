/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_FORMEXPRBCEXPRSUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_FORMEXPRBCEXPRSUM_H

#include <memory>

#include "Rodin/Utility/IsSpecialization.h"

#include "TypeTraits.h"
#include "ForwardDecls.h"

#include "../Problem.h"

namespace Rodin::Variational::FormLanguage
{
   template <class Lhs, class Rhs>
   class FormExprBCExprSum
      : public ProblemBody<FormExprBCExprSum<Lhs, Rhs>>
   {
      public:
         FormExprBCExprSum(const Lhs& lhs, const Rhs& rhs)
         {
            m_lhs = std::unique_ptr<Lhs>(lhs.copy());
            m_rhs = std::unique_ptr<Rhs>(rhs.copy());
         }

         FormExprBCExprSum(FormExprBCExprSum&&) = default;

         FormExprBCExprSum(const FormExprBCExprSum& other)
         {
            assert(other.m_lhs);
            m_lhs = std::unique_ptr<Lhs>(other.m_lhs->copy());

            assert(other.m_rhs);
            m_rhs = std::unique_ptr<Rhs>(other.m_rhs->copy());
         }

         virtual FormExprBCExprSum& setProblem(ProblemBase& problem) override
         {
            static_assert(Utility::IsSpecialization<Lhs, BFExprLFExprSum>::value);
            m_lhs->setLinearForm(problem.getLinearForm());
            m_lhs->setBilinearForm(problem.getBilinearForm());

            m_rhs->setProblem(problem);
            return *this;
         }

         virtual void eval() override
         {
            m_lhs->eval();
            m_rhs->eval();
         }

         const Lhs& lhs() const
         {
            return *m_lhs;
         }

         const Rhs& rhs() const
         {
            return *m_rhs;
         }

         virtual FormExprBCExprSum* copy() const noexcept override
         {
            return new FormExprBCExprSum(*this);
         }

      private:
         std::unique_ptr<Lhs> m_lhs;
         std::unique_ptr<Rhs> m_rhs;
   };

   template <class LhsNested, class RhsNested, class Rhs>
   auto
   operator+(const BFExprLFExprSum<LhsNested, RhsNested>& lhs,
         const BCExprList<Rhs>& rhs)
   {
      return FormExprBCExprSum(lhs, rhs);
   }

   template <class LhsNested, class RhsNested, class Rhs>
   auto
   operator+(const FormBCExprSum<LhsNested, RhsNested>& lhs,
         const BCExprList<Rhs>& rhs)
   {
      return FormExprBCExprSum(lhs, rhs);
   }
}

#endif

