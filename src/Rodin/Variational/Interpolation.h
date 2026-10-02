/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTERPOLATION_H
#define RODIN_VARIATIONAL_INTERPOLATION_H

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @brief Selects and evaluates finite-element interpolation functionals.
   *
   * For a space @f$ V_h @f$, interpolation determines coefficients
   * @f$ c_i = \ell_i(f) @f$ from its degree-of-freedom functionals.
   * The space supplies the scalar field, DOF mappings, pullbacks, and mesh
   * context. Context-specific implementations specialize this template;
   * coefficient storage is handled separately by the consuming backend.
   *
   * The MPI specialization is provided by
   * <a href="_m_p_i_2_variational_2_interpolation_8h.html">MPI/Variational/Interpolation.h</a>.
   * The primary template has no implementation for other contexts.
   *
   * @tparam FES Finite element space defining the interpolation functionals.
   */
  template <class FES>
  class Interpolation;
}

#endif
