/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SparseQR.h
 * @brief Sparse QR factorization with column pivoting.
 *
 * This header provides the SparseQR solver class, which implements QR decomposition
 * with numerical column pivoting for sparse matrices.
 *
 * ## Algorithm
 * The solver computes:
 * @f[
 *   AP = QR
 * @f]
 * where @f$ P @f$ is a permutation matrix, @f$ Q @f$ is orthogonal, and @f$ R @f$ is
 * upper triangular. The left-looking algorithm proceeds column by column.
 *
 * ## Applicability
 * - General sparse matrices (no symmetry required)
 * - Rank-deficient systems
 * - Least-squares problems
 * - Overdetermined systems
 * - Systems where LU factorization is numerically unstable
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * 
 * Solver::SparseQR solver(problem);
 * solver.solve();
 * ```
 *
 * @see SparseQR for the solver implementation
 */
#ifndef RODIN_SOLVER_SPARSEQR_H
#define RODIN_SOLVER_SPARSEQR_H

#include <Eigen/SparseQR>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "LinearSolver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup SparseQRSpecializations SparseQR Template Specializations
   * @brief Template specializations of the SparseQR class.
   * @see SparseQR
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for SparseQR
   */
  template <class LinearSystem>
  SparseQR(Variational::ProblemBase<LinearSystem>&) -> SparseQR<LinearSystem>;

  /**
   * @ingroup SparseQRSpecializations
   * @brief Sparse QR factorization with numerical column pivoting.
   *
   * This solver performs sparse QR decomposition using a left-looking algorithm
   * with COLAMD ordering for column pivoting. It is suitable for general sparse
   * matrices and handles rank-deficient systems.
   *
   * **Characteristics:**
   * - **Type**: Direct solver
   * - **Matrix requirements**: General sparse (no symmetry required)
   * - **Memory**: @f$ O(nnz) @f$ for the factorization
   * - **Stability**: Excellent numerical stability
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class SparseQR<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
    : public LinearSolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// Type of scalar values in the system
      using ScalarType = Scalar;

      /// Type of solution and right-hand side vectors
      using VectorType = Math::Vector<ScalarType>;

      /// Type of system matrix (operator)
      using OperatorType = Math::SparseMatrix<ScalarType>;

      /// Type of linear system being solved
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      /// Type of problem being solved
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      /// Parent class type
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /**
       * @brief Constructs a SparseQR solver for the given problem.
       * @param pb Reference to the problem to solve
       */
      SparseQR(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      SparseQR(const SparseQR& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      SparseQR(SparseQR&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Solves the linear system using sparse QR factorization.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       *
       * This method computes the QR factorization with column pivoting and
       * uses it to solve the system via triangular solves.
       */
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new SparseQR instance
       */
      SparseQR* copy() const noexcept override
      {
        return new SparseQR(*this);
      }

    private:
      /// Underlying Eigen SparseQR solver with COLAMD ordering
      Eigen::SparseQR<OperatorType, Eigen::COLAMDOrdering<int>> m_solver;
  };
}

#endif



