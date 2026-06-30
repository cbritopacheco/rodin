/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ADAPTATION_WNGIRBACKEND_H
#define RODIN_PETSC_ADAPTATION_WNGIRBACKEND_H

#include <petscvec.h>
#include <type_traits>

#include "Rodin/Adaptation/WNGIRBackend.h"
#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/PETSc/Variational.h"

namespace Rodin::Adaptation
{
  /**
   * @brief PETSc specialization of @ref WNGIRBackend for a displacement stored
   *        as a @c ::Vec -backed GridFunction.
   *
   * Selects the PETSc trial/test functions and a @c PETSc::Math::LinearSystem
   * step problem, so a WNGIR solver constructed from a PETSc displacement runs
   * its whole pipeline on PETSc objects with no per-backend code in the solver.
   */
  namespace Internal
  {
    template <class FES>
    struct PETScWNGIRBackend
    {
      using FESType = FES;

      using TrialFunctionType = std::decay_t<
        decltype(PETSc::Variational::TrialFunction(std::declval<const FES&>()))>;

      using TestFunctionType = std::decay_t<
        decltype(PETSc::Variational::TestFunction(std::declval<const FES&>()))>;

      using ProblemType = Variational::Problem<
        PETSc::Math::LinearSystem, TrialFunctionType, TestFunctionType>;
    };
  }

  /// PETSc backend for a displacement stored as the base @c ::Vec -backed
  /// GridFunction specialization.
  template <class FES>
  struct WNGIRBackend<Variational::GridFunction<FES, ::Vec>>
    : Internal::PETScWNGIRBackend<FES>
  {};

  /// PETSc backend for a displacement declared with the convenience
  /// @ref PETSc::Variational::GridFunction class (derived from the base above).
  template <class FES>
  struct WNGIRBackend<PETSc::Variational::GridFunction<FES>>
    : Internal::PETScWNGIRBackend<FES>
  {};
}

#endif
