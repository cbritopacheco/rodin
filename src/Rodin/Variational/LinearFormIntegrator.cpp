/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "LinearFormIntegrator.h"

namespace Rodin::Variational
{
   std::unique_ptr<mfem::LinearFormIntegrator> LinearFormIntegratorBase::build() const
   {
      return std::make_unique<Internal::ProxyLinearFormIntegrator>(*this);
   }
}
