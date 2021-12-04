/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BFEXPRBCEXPRSUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BFEXPRBCEXPRSUM_H

#include <memory>

#include "Rodin/Utility/IsSpecialization.h"

#include "TypeTraits.h"
#include "ForwardDecls.h"


#include "../Problem.h"

namespace Rodin::Variational::FormLanguage
{
   // template <class Lhs, class Rhs>
   // struct TypeTraits<BFExprBCExprSum<Lhs, Rhs>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Constructor;
   //    using Rule = ProblemBody<BFExprBCExprSum<Lhs, Rhs>>;
   // };

   template <class Lhs, class Rhs>
   class BFExprBCExprSum
      : public ProblemBody<BFExprBCExprSum<Lhs, Rhs>>
   {
      public:
         BFExprBCExprSum(const Lhs& lhs, const Rhs& rhs)
         {
            m_lhs = std::unique_ptr<Lhs>(lhs.copy());
            m_rhs = std::unique_ptr<Rhs>(rhs.copy());
         }

         BFExprBCExprSum(BFExprBCExprSum&&) = default;

         BFExprBCExprSum(const BFExprBCExprSum& other)
         {
            assert(other.m_lhs);
            m_lhs = std::unique_ptr<Lhs>(other.m_lhs->copy());

            assert(other.m_rhs);
            m_rhs = std::unique_ptr<Rhs>(other.m_rhs->copy());
         }

         virtual BFExprBCExprSum& setProblem(ProblemBase& problem) override
         {
            if constexpr (Utility::IsSpecialization<Lhs, BilinearFormExpr>::value)
            {
               m_lhs->setBilinearForm(problem.getBilinearForm());
            }
            else
            {
               // We should be hitting the rest of the BFExprBCExprSum instead
               // of the trailing BCExprList
               m_lhs->setProblem(problem);
            }

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

         template <class ... Args>
         static BFExprBCExprSum* create(Args&&... args) noexcept
         {
            return new BFExprBCExprSum(std::forward<Args>(args)...);
         }

         virtual BFExprBCExprSum* copy() const noexcept override
         {
            return new BFExprBCExprSum(*this);
         }

      private:
         std::unique_ptr<Lhs> m_lhs;
         std::unique_ptr<Rhs> m_rhs;
   };

   template <class Lhs, class Rhs>
   auto
   operator+(const BilinearFormExpr<Lhs>& lhs, const BCExprList<Rhs>& rhs)
   {
      return BFExprBCExprSum(lhs, rhs);
   }

   template <class LhsNested, class RhsNested, class Rhs>
   auto
   operator+(const BFExprBCExprSum<LhsNested, RhsNested>& lhs,
         const BCExprList<Rhs>& rhs)
   {
      return BFExprBCExprSum(lhs, rhs);
   }
}

#endif

