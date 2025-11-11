/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SPQR.h
 * @brief SuiteSparseQR multifrontal sparse QR solver wrapper.
 *
 * This header provides a wrapper for SPQR (SuiteSparseQR) from SuiteSparse,
 * a high-performance direct solver using multifrontal sparse QR factorization.
 *
 * ## Algorithm
 * SPQR computes:
 * @f[
 *   AP = QR
 * @f]
 * where @f$ P @f$ is a fill-reducing permutation, @f$ Q @f$ is orthogonal,
 * and @f$ R @f$ is upper triangular. The multifrontal method processes the
 * matrix efficiently using an elimination tree.
 *
 * ## Applicability
 * - General sparse matrices
 * - Rank-deficient systems
 * - Least-squares problems
 * - Overdetermined systems
 * - Problems requiring excellent numerical stability
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * 
 * Solver::SPQR solver(problem);
 * solver.setTolerance(1e-10).solve();
 * ```
 *
 * @note This solver requires SPQR from SuiteSparse to be installed
 * and RODIN_USE_SPQR to be defined at compile time.
 *
 * @see SPQR for the solver implementation
 */
#ifndef RODIN_SOLVER_SPQR_H
#define RODIN_SOLVER_SPQR_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_SPQR

#include <optional>
#include <functional>

#include <Eigen/SPQRSupport>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup SPQRSpecializations SPQR Template Specializations
   * @brief Template specializations of the SPQR class.
   * @see SPQR
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for SPQR
   */
  template <class LinearSystem>
  SPQR(Variational::ProblemBase<LinearSystem>&) -> SPQR<LinearSystem>;

  /**
   * @ingroup SPQRSpecializations
   * @brief SuiteSparseQR multifrontal sparse QR solver.
   *
   * This solver provides access to SPQR (SuiteSparseQR), a high-performance
   * QR factorization solver from SuiteSparse. It uses multifrontal methods
   * and is particularly effective for rank-deficient and least-squares problems.
   *
   * **Characteristics:**
   * - **Type**: Direct solver
   * - **Matrix requirements**: General sparse
   * - **Memory**: Efficient for sparse systems
   * - **Stability**: Excellent, uses QR factorization
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class SPQR<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
    : public SolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
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
      using Parent = SolverBase<LinearSystemType>;

      using Parent::solve;

      /**
       * @brief Constructs a SPQR solver for the given problem.
       * @param pb Reference to the problem to solve
       */
      SPQR(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      SPQR(const SPQR& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      SPQR(SPQR&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Default destructor.
       */
      ~SPQR() = default;

      /**
       * @brief Solves the linear system using SPQR.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      /**
       * @brief Sets the tolerance for rank detection.
       * @param tol Tolerance for determining rank
       * @returns Reference to this solver for method chaining
       */
      SPQR& setTolerance(Real tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /**
       * @brief Sets the maximum number of iterations.
       * @param maxIt Maximum number of iterations
       * @returns Reference to this solver for method chaining
       */
      SPQR& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new SPQR instance
       */
      SPQR* copy() const noexcept override
      {
        return new SPQR(*this);
      }

    private:
      /// Underlying Eigen SPQR solver
      Eigen::SPQR<OperatorType> m_solver;
  };
}

#endif
#endif




