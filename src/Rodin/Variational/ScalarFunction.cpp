/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Restriction.h"

#include "ScalarFunction.h"

namespace Rodin::Variational
{
   std::unique_ptr<mfem::Coefficient> ScalarFunctionBase::build() const
   {
      return std::unique_ptr<mfem::Coefficient>(new Internal::ProxyScalarFunction(*this));
   }

   Restriction<ScalarFunctionBase> ScalarFunctionBase::restrictTo(int attr)
   {
      return restrictTo(std::set<int>{attr});
   }

   Restriction<ScalarFunctionBase> ScalarFunctionBase::restrictTo(
         const std::set<int>& attrs)
   {
      return Restriction<ScalarFunctionBase>(*this).to(attrs);
   }
}
