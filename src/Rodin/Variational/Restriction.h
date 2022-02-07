/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_RESTRICTION_H
#define RODIN_VARIATIONAL_RESTRICTION_H

#include <set>
#include "ForwardDecls.h"
#include "ScalarCoefficient.h"
#include "FormLanguage/Base.h"

namespace Rodin::Variational
{
   class RestrictionBase : public FormLanguage::Base
   {
      public:
         virtual const std::set<int>& getAttributes() const = 0;
   };

   template <>
   class Restriction<RestrictionBase>
   {};

   template <>
   class Restriction<ScalarCoefficientBase> : public RestrictionBase
   {
      public:
         Restriction(const ScalarCoefficientBase& s);

         Restriction(const Restriction& other);

         Restriction& to(int attr)
         {
            return to(std::set<int>{attr});
         }

         Restriction& to(const std::set<int>& attrs);

         const ScalarCoefficientBase& getScalarCoefficient() const;

         const std::set<int>& getAttributes() const override;

         Restriction* copy() const noexcept override
         {
            return new Restriction(*this);
         }

      private:
         std::set<int> m_attr;
         std::unique_ptr<ScalarCoefficientBase> m_s;
   };

   using RestrictedScalarCoefficient = Restriction<ScalarCoefficientBase>;
}

#endif
