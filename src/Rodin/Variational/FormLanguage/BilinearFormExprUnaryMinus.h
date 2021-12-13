/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BILINEARFORMEXPRUNARYMINUS_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BILINEARFORMEXPRUNARYMINUS_H

#include <memory>
#include <variant>

#include "ForwardDecls.h"

#include "BilinearFormExpr.h"

namespace Rodin::Variational::FormLanguage
{
   template <class T>
   class BilinearFormExprUnaryMinus
      : public BilinearFormExpr<BilinearFormExprUnaryMinus<T>>
   {
      public:
         BilinearFormExprUnaryMinus(const T& v)
         {
            m_v = std::unique_ptr<T>(v.copy());
         }

         BilinearFormExprUnaryMinus(BilinearFormExprUnaryMinus&& other) = default;

         BilinearFormExprUnaryMinus(const BilinearFormExprUnaryMinus& other)
         {
            assert(other.m_v);
            m_v = std::unique_ptr<T>(other.m_v->copy());
         }

         virtual BilinearFormExprUnaryMinus& setBilinearForm(BilinearFormBase& bf)
         override
         {
            m_v->setBilinearForm(bf);
            return *this;
         }

         virtual void eval() override
         {
            m_v->toggleSign().eval();
         }

         const T& v() const
         {
            return *m_v;
         }

         template <class ... Args>
         static BilinearFormExprUnaryMinus* create(Args&&... args) noexcept
         {
            return new BilinearFormExprUnaryMinus(std::forward<Args>(args)...);
         }

         virtual BilinearFormExprUnaryMinus* copy() const noexcept override
         {
            return new BilinearFormExprUnaryMinus(*this);
         }

      private:
         std::unique_ptr<T> m_v;
   };

   template <class T>
   auto
   operator-(const BilinearFormExpr<T>& v)
   {
      return BilinearFormExprUnaryMinus(v);
   }

   template <class Lhs, class Rhs>
   auto
   operator-(const BilinearFormExpr<Lhs>& lhs, const BilinearFormExpr<Rhs>& rhs)
   {
      return BilinearFormExprSum(lhs, BilinearFormExprUnaryMinus(rhs));
   }
}

#endif
