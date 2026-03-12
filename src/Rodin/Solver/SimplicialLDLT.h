/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SimplicialLDLT.h
 * @brief Simplicial LDLT Cholesky factorization for sparse SPD matrices.
 *
 * This header provides the SimplicialLDLT solver class, which implements Cholesky
 * decomposition without square root for symmetric positive definite sparse matrices.
 *
 * ## Algorithm
 * The solver computes:
 * @f[
 *   A = LDL^T
 * @f]
 * where @f$ L @f$ is unit lower triangular and @f$ D @f$ is diagonal. This variant
 * avoids square root operations, improving numerical stability.
 *
 * ## Applicability
 * - Symmetric positive definite sparse matrices
 * - When numerical stability is a concern
 * - Problems requiring higher precision
 * - Systems where LLT might encounter numerical issues
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * 
 * Solver::SimplicialLDLT solver(problem);
 * solver.solve();
 * ```
 *
 * @see SimplicialLDLT for the solver implementation
 */
#ifndef RODIN_SOLVER_SIMPLICIALLDLT_H
#define RODIN_SOLVER_SIMPLICIALLDLT_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "LinearSolver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup SimplicialLDLTSpecializations SimplicialLDLT Template Specializations
   * @brief Template specializations of the SimplicialLDLT class.
   * @see SimplicialLDLT
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for SimplicialLDLT
   */
  template <class LinearSystemType>
  SimplicialLDLT(Variational::ProblemBase<LinearSystemType>&) -> SimplicialLDLT<LinearSystemType>;

  /**
   * @ingroup SimplicialLDLTSpecializations
   * @brief Simplicial LDLT Cholesky factorization without square root.
   *
   * This solver performs Cholesky decomposition @f$ A = LDL^T @f$ using a simplicial
   * algorithm. It avoids square root operations, providing better numerical stability
   * than LLT for ill-conditioned matrices.
   *
   * **Characteristics:**
   * - **Type**: Direct solver
   * - **Matrix requirements**: Symmetric positive definite sparse
   * - **Memory**: @f$ O(nnz) @f$ for the factorization
   * - **Stability**: Better than LLT, avoids square roots
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class SimplicialLDLT<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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
       * @brief Constructs a SimplicialLDLT solver for the given problem.
       * @param pb Reference to the problem to solve
       */
      SimplicialLDLT(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      SimplicialLDLT(const SimplicialLDLT& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      SimplicialLDLT(SimplicialLDLT&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Solves the linear system using simplicial LDLT factorization.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new SimplicialLDLT instance
       */
      inline
      SimplicialLDLT* copy() const noexcept override
      {
        return new SimplicialLDLT(*this);
      }

    private:
      /// Underlying Eigen SimplicialLDLT solver
      Eigen::SimplicialLDLT<OperatorType> m_solver;
  };
}

#endif


