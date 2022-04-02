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
#include "ScalarFunction.h"
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
   class Restriction<ScalarFunctionBase> : public RestrictionBase
   {
      public:
         Restriction(const ScalarFunctionBase& s);

         Restriction(const Restriction& other);

         Restriction& to(int attr)
         {
            return to(std::set<int>{attr});
         }

         Restriction& to(const std::set<int>& attrs);

         const ScalarFunctionBase& getScalarFunction() const;

         const std::set<int>& getAttributes() const override;

         Restriction* copy() const noexcept override
         {
            return new Restriction(*this);
         }

      private:
         std::set<int> m_attr;
         std::unique_ptr<ScalarFunctionBase> m_s;
   };

   using RestrictedScalarFunction = Restriction<ScalarFunctionBase>;
}

#endif
