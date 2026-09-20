/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Info.h
 * @brief Status reporting shared by the factorization-based solvers.
 */
#ifndef RODIN_SOLVER_INFO_H
#define RODIN_SOLVER_INFO_H

#include "Rodin/Types.h"

#include "LinearSolver.h"

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
   * @see FactorizationSolverBase
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

  /**
   * @brief Base class of the solvers that factorize their operator.
   *
   * A factorization can fail for reasons that are not programming errors, such
   * as a singular matrix or a factor that does not fit in memory. Solving with
   * a failed factorization yields no solution and, depending on the backend,
   * either leaves the solution vector untouched or reads state that was never
   * computed. These solvers therefore stop before the solve and report what
   * happened through @ref getInfo, leaving the caller to act on it.
   *
   * A violated contract, such as a matrix that is not square, remains an
   * exception: it is a defect in the calling code rather than an outcome.
   *
   * @tparam LinearSystem Type of linear system to solve.
   */
  template <class LinearSystem>
  class FactorizationSolverBase : public LinearSolverBase<LinearSystem>
  {
    public:
      /// @brief Parent class type.
      using Parent = LinearSolverBase<LinearSystem>;

      /// @brief Problem type solved by this solver.
      using ProblemBaseType = Variational::ProblemBase<LinearSystem>;

      using Parent::solve;

      /// @brief Constructs the solver for the given problem.
      FactorizationSolverBase(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /// @brief Copy constructor. The copy retains no factorization.
      FactorizationSolverBase(const FactorizationSolverBase& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      FactorizationSolverBase(FactorizationSolverBase&& other) noexcept
        : Parent(std::move(other)),
          m_info(other.m_info)
      {}

      /// @brief Default virtual destructor.
      virtual ~FactorizationSolverBase() = default;

      /**
       * @brief Returns the outcome of the most recent operation.
       *
       * The returned @ref Info is updated by every factorization and every
       * solve, so a caller reads it after the call it wants to check.
       */
      const Info& getInfo() const noexcept
      {
        return m_info;
      }

      /**
       * @brief Checks whether the most recent operation succeeded.
       * @returns true if the solver succeeded, false otherwise.
       */
      Boolean success() const noexcept
      {
        return m_info.success;
      }

    protected:
      /// @brief Outcome of the most recent operation.
      Info m_info;
  };
}

#endif
