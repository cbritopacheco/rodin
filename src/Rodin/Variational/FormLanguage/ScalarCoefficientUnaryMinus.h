/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_COEFFUNARYMINUS_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_COEFFUNARYMINUS_H

#include "ForwardDecls.h"

#include "../ScalarCoefficient.h"

namespace Rodin::Variational::FormLanguage
{
   template <class T>
   class ScalarCoefficientUnaryMinus
      : public ScalarCoefficient<ScalarCoefficientUnaryMinus<T>>
   {
      public:
         ScalarCoefficientUnaryMinus(const T& v)
         {
            m_v = std::unique_ptr<T>(v.copy());
         }

         ScalarCoefficientUnaryMinus(const ScalarCoefficientUnaryMinus& other)
         {
            assert(other.m_v);
            m_v = std::unique_ptr<T>(other.m_v->copy());
         }

         ScalarCoefficientUnaryMinus(ScalarCoefficientUnaryMinus&&) = default;

         const T& v() const
         {
            return *m_v;
         }

         template <class ... Args>
         static ScalarCoefficientUnaryMinus* create(Args&&... args) noexcept
         {
            return new ScalarCoefficientUnaryMinus(std::forward<Args>(args)...);
         }

         virtual ScalarCoefficientUnaryMinus* copy() const noexcept override
         {
            return new ScalarCoefficientUnaryMinus(*this);
         }

      private:
         std::unique_ptr<T> m_v;
   };

   template <class T>
   auto operator-(const CoeffExpr<T>& v)
   {
      return ScalarCoefficientUnaryMinus(v);
   }

   template <class Lhs, class Rhs>
   auto operator-(const ScalarCoefficient<Lhs>& lhs, const ScalarCoefficient<Rhs>& rhs)
   {
      return ScalarCoefficientSum(lhs, ScalarCoefficientUnaryMinus(rhs));
   }
}

#endif
