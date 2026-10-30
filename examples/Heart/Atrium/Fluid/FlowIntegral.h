/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or copy at https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_FLUID_FLOWINTEGRAL_H
#define RODIN_ATRIUM_FLUID_FLOWINTEGRAL_H

#include <utility>
#include <Rodin/Variational/Integral.h>

namespace Rodin::Examples::Heart::Atrium
{
  /// @brief Matches Oseen and material-derivative volume quadrature at order four.
  template <class... Terms>
  auto FlowIntegral(Terms&&... terms)
  {
    auto form = Variational::Integral(std::forward<Terms>(terms)...);
    form.setOrder(4);
    return form;
  }
}

#endif
