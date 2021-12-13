/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_SCALARCOEFFSUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_SCALARCOEFFSUM_H

#include <memory>
#include <utility>

#include "ForwardDecls.h"

#include "../ScalarCoefficient.h"

namespace Rodin::Variational::FormLanguage
{
   template <class Lhs, class Rhs>
   class ScalarCoefficientSum : public ScalarCoefficient<ScalarCoefficientSum<Lhs, Rhs>>
   {
      public:
         ScalarCoefficientSum(
               const ScalarCoefficient<Lhs>& lhs, const ScalarCoefficient<Rhs>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         ScalarCoefficientSum(const ScalarCoefficientSum& other)
            : m_lhs(other.m_lhs.copy()), m_rhs(other.m_rhs.copy())
         {}

         const ScalarCoefficient<Lhs>& getLHS() const
         {
            assert(m_lhs);
            return *m_lhs;
         }

         const ScalarCoefficient<Rhs>& getRHS() const
         {
            assert(m_rhs);
            return *m_rhs;
         }

         template <class ... Args>
         static ScalarCoefficientSum* create(Args&&... args) noexcept
         {
            return new ScalarCoefficientSum(std::forward<Args>(args)...);
         }

         virtual ScalarCoefficientSum* copy() const noexcept override
         {
            return new ScalarCoefficientSum(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficient<Lhs>> m_lhs;
         std::unique_ptr<ScalarCoefficient<Rhs>> m_rhs;
   };

   template <class Lhs, class Rhs>
   auto operator+(const ScalarCoefficient<Lhs>& lhs, const ScalarCoefficient<Rhs>& rhs)
   {
      return ScalarCoefficientSum(lhs, rhs);
   }
}

#endif
