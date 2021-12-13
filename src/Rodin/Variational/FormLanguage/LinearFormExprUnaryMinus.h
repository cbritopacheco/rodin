/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_LINEARFORMEXPRUNARYMINUS_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_LINEARFORMEXPRUNARYMINUS_H

#include <memory>
#include <variant>

#include "ForwardDecls.h"

#include "LinearFormExpr.h"

namespace Rodin::Variational::FormLanguage
{
   template <class T>
   class LinearFormExprUnaryMinus
      : public LinearFormExpr<LinearFormExprUnaryMinus<T>>
   {
      public:
         LinearFormExprUnaryMinus(const T& v)
         {
            m_v = std::unique_ptr<T>(v.copy());
         }

         LinearFormExprUnaryMinus(LinearFormExprUnaryMinus&& other) = default;

         LinearFormExprUnaryMinus(const LinearFormExprUnaryMinus& other)
         {
            assert(other.m_v);
            m_v = std::unique_ptr<T>(other.m_v->copy());
         }

         virtual LinearFormExprUnaryMinus& setLinearForm(LinearFormBase& bf) override
         {
            assert(m_v);
            m_v->setLinearForm(bf);
            return *this;
         }

         virtual void eval() override
         {
            assert(m_v);
            m_v->toggleSign().eval();
         }

         const T& v() const
         {
            assert(m_v);
            return *m_v;
         }

         virtual LinearFormExprUnaryMinus* copy() const noexcept override
         {
            return new LinearFormExprUnaryMinus(*this);
         }

      private:
         std::unique_ptr<T> m_v;
   };

   template <class T>
   auto
   operator-(const LinearFormExpr<T>& v)
   {
      return LinearFormExprUnaryMinus(v);
   }

   template <class Lhs, class Rhs>
   auto
   operator-(const LinearFormExpr<Lhs>& lhs, const LinearFormExpr<Rhs>& rhs)
   {
      return LinearFormExprSum(lhs, LinearFormExprUnaryMinus(rhs));
   }
}

#endif

