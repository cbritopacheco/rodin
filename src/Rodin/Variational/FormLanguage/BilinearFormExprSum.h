/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BILINEARFORMEXPRSUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BILINEARFORMEXPRSUM_H

#include <memory>
#include <variant>

#include "ForwardDecls.h"

#include "BilinearFormExpr.h"

namespace Rodin::Variational::FormLanguage
{
   /**
    * @brief Represents an expression between two objects of type
    * BilinearFormExpr.
    *
    * @code{.unparsed}
    *    BilinearFormExpr + BilinearFormExpr
    * @endcode
    */
   template <class Lhs, class Rhs>
   class BilinearFormExprSum
      : public BilinearFormExpr<BilinearFormExprSum<Lhs, Rhs>>
   {
      public:
         /**
          * Constructs a BilinearFormExprSum from two BilinearFormExpr.
          *
          * @todo Add static_assert to make sure Lhs and Rhs are of the
          * appropiate type.
          */
         BilinearFormExprSum(const Lhs& lhs, const Rhs& rhs)
         {
            m_lhs = std::unique_ptr<Lhs>(lhs.copy());
            m_rhs = std::unique_ptr<Rhs>(rhs.copy());
         }

         /**
          * @brief Default move constructor.
          */
         BilinearFormExprSum(BilinearFormExprSum&& other) = default;

         /**
          * @brief Copies the left and right children.
          */
         BilinearFormExprSum(const BilinearFormExprSum& other)
         {
            assert(other.m_lhs);
            m_lhs = std::unique_ptr<Lhs>(other.m_lhs->copy());

            assert(other.m_rhs);
            m_rhs = std::unique_ptr<Rhs>(other.m_rhs->copy());
         }

         /**
          * @brief Sets the bilinear form to build from the expression.
          * @returns Reference to self (for method chaining).
          */
         virtual BilinearFormExprSum& setBilinearForm(BilinearFormBase& bf) override
         {
            m_lhs->setBilinear(bf);
            m_rhs->setBilinear(bf);
            return *this;
         }

         /**
          * @brief Evaluates the BilinearFormExpr
          */
         virtual void eval() override
         {
            m_lhs->eval();
            m_rhs->eval();
         }

         /**
          * @brief Toggles the sign that is in front of the sum.
          *
          * This function has the semantic operation to toggle the sign of the
          * whole expression. That is, it will apply the distributive operation
          * of the minus operator so that the evaluation satisfies:
          *
          * @code{.unparsed}
          *    Eval[-(lhs + rhs)] = Eval[(-lhs) + (-rhs)]
          * @endcode
          */
         virtual BilinearFormExprSum& toggleSign()
         {
            m_lhs->toggleSign();
            m_rhs->toggleSign();
            return *this;
         }

         const Lhs& lhs() const
         {
            return *m_lhs;
         }

         const Rhs& rhs() const
         {
            return *m_rhs;
         }

         virtual BilinearFormExprSum* copy() const noexcept override
         {
            return new BilinearFormExprSum(*this);
         }

      private:
         std::unique_ptr<Lhs> m_lhs;
         std::unique_ptr<Rhs> m_rhs;
   };

   template <class Lhs, class Rhs>
   auto
   operator+(const BilinearFormExpr<Lhs>& lhs, const BilinearFormExpr<Rhs>& rhs)
   {
      return BilinearFormExprSum(lhs, rhs);
   }
}

#endif

