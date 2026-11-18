/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_KOKKOS_VARIATIONAL_TRIALFUNCTION_H

/**
 * @file
 * @brief Backend Data deduction for trial functions.
 */

#include "Rodin/Variational/TrialFunction.h"
#include "Rodin/Variational/GridFunction.h"
#include "DOFs.h"
#include "GridFunctionData.h"

namespace Rodin::Variational
{
  /// @param fes Kokkos space. @return Trial function with backend coefficients.
  template <::Rodin::Kokkos::Variational::Space FES>
  TrialFunction(const FES& fes)
    -> TrialFunction<GridFunction<FES, ::Rodin::Kokkos::Variational::Data<FES>>, FES>;

}

#endif
