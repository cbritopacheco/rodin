/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_CONVERGENCE_PETSC_POISSONBOUNDARYPROBLEM_H
#define RODIN_TESTS_CONVERGENCE_PETSC_POISSONBOUNDARYPROBLEM_H

#include "PETScDiffusionBoundaryProblem.h"

namespace Rodin::Tests::Convergence
{
  /** Poisson is the constant-coefficient scalar diffusion specialization. */
  template <size_t K, class MeshType>
  using PETScPoissonBoundaryProblem = PETScDiffusionBoundaryProblem<K, MeshType, true>;
}

#endif
