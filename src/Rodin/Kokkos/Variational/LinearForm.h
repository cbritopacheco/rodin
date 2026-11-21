/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_LINEARFORM_H
#define RODIN_KOKKOS_VARIATIONAL_LINEARFORM_H

/**
 * @file
 * @brief Backend algebra deduction for linear forms.
 */

#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/TestFunction.h"
#include "DOFs.h"
#include "GridFunctionData.h"

namespace Rodin::Variational
{
  /// @param v Test function on a Kokkos space.
  template <::Rodin::Kokkos::Variational::Space FES>
  LinearForm(const TestFunction<FES>& v)
    -> LinearForm<FES, ::Rodin::Kokkos::Variational::Data<FES>>;

}

#endif
