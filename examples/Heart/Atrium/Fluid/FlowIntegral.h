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
  /**
   * @brief Volume quadrature shared by the Oseen residual and its material derivative.
   *
   * Viscosity, cell size, and stabilization parameters are constant on a
   * cell, so every polynomial integrand of the affine-P1 Oseen form has
   * degree at most two and the four-point rule integrates it exactly. The
   * fibrin drag is the only non-polynomial coefficient; the residual and its
   * derivative sample it with the same rule.
   */
  template <class... Terms>
  auto FlowIntegral(Terms&&... terms)
  {
    auto form = Variational::Integral(std::forward<Terms>(terms)...);
    form.setOrder(2);
    return form;
  }

  /**
   * @brief Volume quadrature shared by every term of the species transport.
   *
   * Order two is exact for the affine-P1 mass, advection, and diffusion terms.
   * Using the same rule for the history and the unknown makes the pointwise
   * stabilization terms cancel for a constant state, so the transport
   * preserves constants and the invariant of the kinetics.
   */
  template <class... Terms>
  auto SpeciesIntegral(Terms&&... terms)
  {
    auto form = Variational::Integral(std::forward<Terms>(terms)...);
    form.setOrder(2);
    return form;
  }
}

#endif
