/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_KOKKOS_VARIATIONAL_GRIDFUNCTION_H

/**
 * @file
 * @brief Backend Data deduction for grid functions.
 */

#include "Rodin/Variational/GridFunction.h"
#include "DOFs.h"
#include "GridFunctionData.h"
// Retain the historical convenience includes for existing clients.
#include "Integral.h"
#include "BilinearForm.h"
#include "LinearForm.h"
#include "TrialFunction.h"
#include "Problem.h"

namespace Rodin::Variational
{
  /// @param fes Kokkos space. @return Grid function with backend coefficients.
  template <::Rodin::Kokkos::Variational::Space FES>
  GridFunction(
    const FES& fes) -> GridFunction<FES, ::Rodin::Kokkos::Variational::Data<FES>>;

}

#endif
