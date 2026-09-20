/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Info.h
 * @brief Status reported by the factorization-based solvers.
 *
 * The solvers that factorize their operator have no common ancestor: each one
 * wraps a library of its own. They share this vocabulary instead, and expose
 * it uniformly through `getInfo()` and `success()`.
 *
 * A factorization can fail for reasons that are not programming errors, such
 * as a singular matrix or a factor that does not fit in memory. Solving with a
 * failed factorization yields no solution and, depending on the backend,
 * either leaves the solution vector untouched or reads state that was never
 * computed. A solver therefore stops before the solve and reports what
 * happened through its @ref Rodin::Solver::Info "Info", leaving the caller to
 * act on it. A violated contract, such as a matrix that is not square, remains
 * an exception: it is a defect in the calling code rather than an outcome.
 */
#ifndef RODIN_SOLVER_INFO_H
#define RODIN_SOLVER_INFO_H

#include "Rodin/Types.h"

namespace Rodin::Solver
{
  /**
   * @brief Stage of a factorization.
   *
   * Symbolic analysis depends only on the sparsity pattern, whereas the
   * numeric factorization also depends on the values. A solver that retains
   * Factorization::Numeric therefore retains Factorization::Symbolic too.
   */
  enum class Factorization
  {
    Symbolic,
    Numeric
  };

  /**
   * @brief Outcome of the most recent operation of a factorization-based
   * solver.
   *
   * Obtained from the solver's `getInfo()`, whose `success` is also returned
   * on its own by `success()`.
   */
  struct Info
  {
      /// @brief Whether the most recent operation succeeded.
      Boolean success = true;

      /**
       * @brief The most advanced factorization stage currently retained.
       *
       * Empty when no factorization is available, either because none has
       * been computed yet or because the last one failed or was released.
       */
      Optional<Factorization> factorization;

      /**
       * @brief Status reported by the underlying solver.
       *
       * Zero denotes success. The meaning of a nonzero value is the backend's
       * own: each solver documents which status it reports here.
       */
      Integer status = 0;
  };

}

#endif
