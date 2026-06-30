/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRBACKEND_H
#define RODIN_ADAPTATION_WNGIRBACKEND_H

#include <type_traits>

#include "Rodin/Variational/TrialFunction.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Variational/Problem.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Maps a displacement GridFunction type to the backend-matched
   *        trial/test function and step-problem types used by @ref WNGIR.
   *
   * The WNGIR solver is constructed from a displacement GridFunction and infers
   * its entire backend (trial/test functions, step Problem, solver) from that
   * single type. The primary template covers the default Eigen-backed
   * GridFunctions; backend headers (e.g. PETSc) provide specializations.
   *
   * The contract a specialization must satisfy:
   *  - @c TrialFunctionType / @c TestFunction Type : backend-matched
   *    trial/test functions over the displacement's finite element space.
   *  - @c ProblemType : the @ref Rodin::Variational::Problem the WNGIR step is
   *    assembled into and solved with @ref Rodin::Solver::CG.
   */
  template <class Displacement>
  struct WNGIRBackend
  {
    using FESType = std::decay_t<
      decltype(std::declval<Displacement&>().getFiniteElementSpace())>;

    using TrialFunctionType = std::decay_t<
      decltype(Variational::TrialFunction(std::declval<const FESType&>()))>;

    using TestFunctionType = std::decay_t<
      decltype(Variational::TestFunction(std::declval<const FESType&>()))>;

    using ProblemType = std::decay_t<decltype(Variational::Problem(
        std::declval<TrialFunctionType&>(),
        std::declval<TestFunctionType&>()))>;
  };
}

#endif
