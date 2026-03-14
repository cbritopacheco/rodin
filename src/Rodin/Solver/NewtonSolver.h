/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_NEWTONSOLVER_H
#define RODIN_SOLVER_NEWTONSOLVER_H

#include "Rodin/Copyable.h"

namespace Rodin::Solver
{
  /**
   * @brief Base interface for Newton-type nonlinear solvers.
   * @tparam NonlinearSystem Type of nonlinear system handled by the solver.
   */
  template <class NonlinearSystem>
  class NewtonSolverBase : public Copyable
  {
    public:
      using SystemType = NonlinearSystem;
      using Parent = Copyable;

      virtual ~NewtonSolverBase() = default;

      /**
       * @brief Solve a nonlinear system.
       * @param[in,out] system Nonlinear system to solve.
       */
      virtual void solve(SystemType& system) = 0;
  };
}

#endif
