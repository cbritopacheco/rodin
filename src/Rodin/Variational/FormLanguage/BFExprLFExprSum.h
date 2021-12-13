/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BFEXPRLFEXPRSUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BFEXPRLFEXPRSUM_H

#include <memory>

#include "TypeTraits.h"
#include "ForwardDecls.h"

#include "../Problem.h"

namespace Rodin::Variational::FormLanguage
{
   template <class Lhs, class Rhs>
   class BFExprLFExprSum
      : public ProblemBody<BFExprLFExprSum<Lhs, Rhs>>
   {
      public:
         BFExprLFExprSum(const Lhs& lhs, const Rhs& rhs)
         {
            m_lhs = std::unique_ptr<Lhs>(lhs.copy());
            m_rhs = std::unique_ptr<Rhs>(rhs.copy());
         }

         BFExprLFExprSum(BFExprLFExprSum&&) = default;

         BFExprLFExprSum(const BFExprLFExprSum& other)
         {
            assert(other.m_lhs);
            m_lhs = std::unique_ptr<Lhs>(other.m_lhs->copy());

            assert(other.m_rhs);
            m_rhs = std::unique_ptr<Rhs>(other.m_rhs->copy());
         }

         virtual BFExprLFExprSum& setBilinearForm(BilinearFormBase& bf)
         {
            m_lhs->setBilinearForm(bf);
            return *this;
         }

         virtual BFExprLFExprSum& setLinearForm(LinearFormBase& lf)
         {
            m_rhs->setLinearForm(lf);
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

         virtual BFExprLFExprSum* copy() const noexcept override
         {
            return new BFExprLFExprSum(*this);
         }

      private:
         std::unique_ptr<Lhs> m_lhs;
         std::unique_ptr<Rhs> m_rhs;
   };

   template <class Lhs, class Rhs>
   auto
   operator+(const BilinearFormExpr<Lhs>& lhs, const LinearFormExpr<Rhs>& rhs)
   {
      return BFExprLFExprSum(lhs, rhs);
   }

   template <class Lhs, class Rhs>
   auto
   operator-(const BilinearFormExpr<Lhs>& lhs, const LinearFormExpr<Rhs>& rhs)
   {
      return BFExprLFExprSum(lhs, LinearFormExprUnaryMinus(rhs));
   }
}

#endif
