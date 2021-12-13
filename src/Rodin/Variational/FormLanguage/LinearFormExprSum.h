/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_LINEARFORMEXPRSUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_LINEARFORMEXPRSUM_H

#include <memory>
#include <variant>

#include "ForwardDecls.h"

#include "LinearFormExpr.h"

namespace Rodin::Variational::FormLanguage
{
   /**
    * @brief Represents an expression between two objects of type
    * LinearFormExpr.
    *
    * @code{.unparsed}
    *    LinearFormExpr + LinearFormExpr
    * @endcode
    */
   template <class Lhs, class Rhs>
   class LinearFormExprSum
      : public LinearFormExpr<LinearFormExprSum<Lhs, Rhs>>
   {
      public:
         /**
          * Constructs a LinearFormExprSum from two LinearFormExpr.
          *
          * @todo Add static_assert to make sure Lhs and Rhs are of the
          * appropiate type.
          */
         LinearFormExprSum(const Lhs& lhs, const Rhs& rhs)
         {
            m_lhs = std::unique_ptr<Lhs>(lhs.copy());
            m_rhs = std::unique_ptr<Rhs>(rhs.copy());
         }

         /**
          * @brief Default move constructor.
          */
         LinearFormExprSum(LinearFormExprSum&& other) = default;

         /**
          * @brief Copies the left and right children.
          */
         LinearFormExprSum(const LinearFormExprSum& other)
         {
            assert(other.m_lhs);
            m_lhs = std::unique_ptr<Lhs>(other.m_lhs->copy());

            assert(other.m_rhs);
            m_rhs = std::unique_ptr<Rhs>(other.m_rhs->copy());
         }

         /**
          * @brief Sets the Linear form to build from the expression.
          * @returns Reference to self (for method chaining).
          */
         virtual LinearFormExprSum& setLinearForm(LinearFormBase& bf) override
         {
            m_lhs->setLinear(bf);
            m_rhs->setLinear(bf);
            return *this;
         }

         /**
          * @brief Evaluates the LinearFormExpr
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
         virtual LinearFormExprSum& toggleSign()
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

         virtual LinearFormExprSum* copy() const noexcept override
         {
            return new LinearFormExprSum(*this);
         }

      private:
         std::unique_ptr<Lhs> m_lhs;
         std::unique_ptr<Rhs> m_rhs;
   };

   template <class Lhs, class Rhs>
   auto
   operator+(const LinearFormExpr<Lhs>& lhs, const LinearFormExpr<Rhs>& rhs)
   {
      return LinearFormExprSum(lhs, rhs);
   }
}

#endif


