/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Restriction.h"

#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   std::unique_ptr<mfem::Coefficient> ScalarCoefficientBase::build() const
   {
      return std::unique_ptr<mfem::Coefficient>(new Internal::ProxyScalarCoefficient(*this));
   }

   Restriction<ScalarCoefficientBase> ScalarCoefficientBase::restrictTo(int attr)
   {
      return restrictTo(std::set<int>{attr});
   }

   Restriction<ScalarCoefficientBase> ScalarCoefficientBase::restrictTo(
         const std::set<int>& attrs)
   {
      return Restriction<ScalarCoefficientBase>(*this).to(attrs);
   }
}
