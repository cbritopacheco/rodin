/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BCEXPRLISTCONS_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BCEXPRLISTCONS_H

#include <memory>
#include <variant>

#include "TypeTraits.h"
#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   // template <class Head>
   // struct TypeTraits<BCExprListCons<Head>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Constructor;
   //    using Rule = BCExprList<BCExprListCons<Head>>;
   // };

   // template <class Head, class Tail>
   // struct TypeTraits<BCExprListCons<Head, Tail>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Constructor;
   //    using Rule = BCExprList<BCExprListCons<Head, Tail>>;
   // };

   template <class Head>
   class BCExprListCons<Head, void>
      : public BCExprList<BCExprListCons<Head, void>>
   {
      public:
         BCExprListCons(const Head& head)
         {
            m_head = std::unique_ptr<Head>(head.copy());
         }

         BCExprListCons(BCExprListCons&& other) = default;

         BCExprListCons(const BCExprListCons& other)
         {
            assert(other.m_head);
            m_head = std::unique_ptr<Head>(other.m_head.copy());
         }

         virtual void eval() override
         {
            m_head->eval();
         }

         virtual BCExprListCons& setProblem(ProblemBase& problem) override
         {
            m_head->setProblem(problem);
            return *this;
         }

         const Head& head() const
         {
            assert(m_head);
            return *m_head;
         }

         template <class ... Args>
         static BCExprListCons* create(Args&&... args) noexcept
         {
            return new BCExprListCons(std::forward<Args>(args)...);
         }

         virtual BCExprListCons* copy() const noexcept override
         {
            return new BCExprListCons(*this);
         }

      private:
         std::unique_ptr<Head> m_head;
   };

   template <class Head, class Tail>
   class BCExprListCons
      : public BCExprList<BCExprListCons<Head, Tail>>
   {
      public:
         BCExprListCons(const Head& head, const Tail& tail)
         {
            m_head = std::unique_ptr<Head>(head.copy());
            m_tail = std::unique_ptr<Tail>(tail.copy());
         }

         BCExprListCons(BCExprListCons&& other) = default;

         BCExprListCons(const BCExprListCons& other)
         {
            assert(other.m_head);
            m_head = std::unique_ptr<Head>(other.m_head.copy());

            assert(other.m_tail);
            m_tail = std::unique_ptr<Tail>(other.m_tail.copy());
         }

         const Head& head() const
         {
            assert(m_head);
            return *m_head;
         }

         const Tail& tail() const
         {
            assert(m_head);
            return *m_tail;
         }

         virtual void eval() override
         {
            m_head->eval();
         }

         virtual BCExprListCons& setProblem(ProblemBase& problem) override
         {
            m_head->setProblem(problem);
            m_tail->setProblem(problem);
            return *this;
         }

         template <class ... Args>
         static BCExprListCons* create(Args&&... args) noexcept
         {
            return new BCExprListCons(std::forward<Args>(args)...);
         }

         virtual BCExprListCons* copy() const noexcept override
         {
            return new BCExprListCons(*this);
         }

      private:
         std::unique_ptr<Head> m_head;
         std::unique_ptr<Tail> m_tail;
   };

   template <class Lhs, class Rhs, class OtherLhs>
   auto
   operator+(const BCExpr<OtherLhs>& lhs, const BCExprListCons<Lhs, Rhs>& rhs)
   {
      return BCExprListCons(lhs, rhs);
   }

   template <class Lhs, class Rhs, class OtherRhs>
   auto
   operator+(const BCExprListCons<Lhs, Rhs>& lhs, const BCExpr<OtherRhs>& rhs)
   {
      return BCExprListCons(rhs, lhs);
   }

   template <class LhsNested, class RhsNested>
   auto
   operator+(const BCExpr<LhsNested>& lhs, const BCExpr<RhsNested>& rhs)
   {
      return BCExprListCons(lhs, rhs);
   }
}

#endif

