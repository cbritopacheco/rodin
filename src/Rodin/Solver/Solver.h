/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SOLVER_H
#define RODIN_SOLVER_SOLVER_H

#include "Rodin/Configure.h"

#include "ForwardDecls.h"

namespace Rodin::Solver
{
  /**
   * @brief Base class for solving linear algebra systems.
   */
  template <class OperatorType, class VectorType>
  class SolverBase
  {
    public:
      /**
       * @brief Default virtual destructor.
       */
      virtual ~SolverBase() = default;

      /**
       * @brief Solves the specified Variational::Problem.
       * @param[in,out] problem Variational problem to solve.
       */
      virtual void solve(OperatorType& A, VectorType& X, VectorType& B) = 0;
  };
}

#endif
