/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LDLT.h
 * @brief Robust LDLT Cholesky factorization for dense matrices.
 *
 * This header provides the LDLT solver class for dense matrices, implementing
 * Cholesky decomposition with pivoting for numerical stability.
 *
 * ## Algorithm
 * The solver computes:
 * @f[
 *   PAP^T = LDL^T
 * @f]
 * where @f$ P @f$ is a permutation matrix, @f$ L @f$ is unit lower triangular,
 * and @f$ D @f$ is diagonal. Pivoting improves numerical stability.
 *
 * ## Applicability
 * - Symmetric indefinite dense matrices
 * - When numerical stability is critical
 * - Small to medium-sized systems
 * - Systems requiring robust factorization
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * 
 * Solver::LDLT solver(problem);
 * solver.solve();
 * ```
 *
 * @see LDLT for the solver implementation
 */
#ifndef RODIN_SOLVER_LDLT_H
#define RODIN_SOLVER_LDLT_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"

#include "ForwardDecls.h"
#include "LinearSolver.h"

namespace Rodin::FormLanguage
{
  template <class LinearSystem>
  struct Traits<Solver::LDLT<LinearSystem>>
  {
    using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup LDLTSpecializations LDLT Template Specializations
   * @brief Template specializations of the LDLT class.
   * @see LDLT
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for LDLT
   */
  template <class LinearSystem>
  LDLT(Variational::ProblemBase<LinearSystem>&) -> LDLT<LinearSystem>;

  /**
   * @ingroup LDLTSpecializations
   * @brief Robust LDLT Cholesky factorization with pivoting for dense matrices.
   *
   * This solver performs LDLT decomposition with robust pivoting for numerical
   * stability. It handles symmetric indefinite matrices and provides reliable
   * solutions for dense systems.
   *
   * **Characteristics:**
   * - **Type**: Direct solver
   * - **Matrix requirements**: Symmetric dense (positive definite or indefinite)
   * - **Memory**: @f$ O(n^2) @f$ for dense storage
   * - **Stability**: Excellent, uses pivoting
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class LDLT<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public LinearSolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// Type of scalar values in the system
      using ScalarType = Scalar;

      /// Type of solution and right-hand side vectors
      using VectorType = Math::Vector<ScalarType>;

      /// Type of system matrix (operator)
      using OperatorType = Math::Matrix<ScalarType>;

      /// Type of linear system being solved
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      /// Type of problem being solved
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      /// Parent class type
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /**
       * @brief Constructs an LDLT solver for the given problem.
       * @param pb Reference to the problem to solve
       */
      LDLT(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      LDLT(const LDLT& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      LDLT(LDLT&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Solves the linear system using LDLT factorization with pivoting.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new LDLT instance
       */
      LDLT* copy() const noexcept override
      {
        return new LDLT(*this);
      }

    private:
      /// Underlying Eigen LDLT solver
      Eigen::LDLT<OperatorType> m_solver;
  };
}

#endif



