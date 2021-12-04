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

#include "TypeTraits.h"
#include "ForwardDecls.h"

#include "../ScalarCoefficient.h"

namespace Rodin::Variational::FormLanguage
{
   // template <class Lhs, class Rhs>
   // struct FormLanguage::TypeTraits<ScalarCoefficientSum<Lhs, Rhs>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Constructor;
   //    using Rule = ScalarCoefficient<ScalarCoefficientSum<Lhs, Rhs>>;
   // };

   template <class Lhs, class Rhs>
   class ScalarCoefficientSum : public ScalarCoefficient<ScalarCoefficientSum<Lhs, Rhs>>
   {
      public:
         ScalarCoefficientSum(const Lhs& lhs, const Rhs& rhs)
         {
            m_lhs = std::unique_ptr<Lhs>(lhs.copy());
            m_rhs = std::unique_ptr<Rhs>(rhs.copy());
         }

         ScalarCoefficientSum(const ScalarCoefficientSum& other)
         {
            assert(other.m_lhs);
            m_lhs = std::unique_ptr<Lhs>(other.m_lhs->copy());

            assert(other.m_rhs);
            m_rhs = std::unique_ptr<Rhs>(other.m_rhs->copy());
         }

         const Lhs& lhs() const
         {
            assert(m_lhs);
            return *m_lhs;
         }

         const Rhs& rhs() const
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
         std::unique_ptr<Lhs> m_lhs;
         std::unique_ptr<Rhs> m_rhs;
   };

   template <class Lhs, class Rhs>
   auto operator+(const ScalarCoefficient<Lhs>& lhs, const ScalarCoefficient<Rhs>& rhs)
   {
      return ScalarCoefficientSum(lhs, rhs);
   }
}

/* Hijack the mfem namespace to add the syntactic sugar for mfem::Coefficient
 * If you are getting errors related to the operator+ maybe you want to check
 * here.
 */
namespace mfem
{
   template <class Lhs, class Rhs>
   std::enable_if_t<
      std::is_base_of_v<mfem::Coefficient, Lhs>
      && std::is_base_of_v<mfem::Coefficient, Rhs>,
      Rodin::Variational::FormLanguage::ScalarCoefficientSum<
         Rodin::Variational::ScalarCoefficient<Lhs>,
         Rodin::Variational::ScalarCoefficient<Rhs>
            >>
   operator+(const Lhs& lhs, const Rhs& rhs)
   {
      return
         Rodin::Variational::FormLanguage::ScalarCoefficientSum<
            Rodin::Variational::ScalarCoefficient<Lhs>,
            Rodin::Variational::ScalarCoefficient<Rhs>
               >(lhs, rhs);
   }
}

#endif
