/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SimplicialLLT.h
 * @brief Simplicial LLT Cholesky factorization for sparse SPD matrices.
 *
 * This header provides the SimplicialLLT solver class, which implements Cholesky
 * decomposition using simplicial factorization for symmetric positive definite
 * sparse matrices.
 *
 * ## Algorithm
 * The solver computes:
 * @f[
 *   A = LL^T
 * @f]
 * where @f$ L @f$ is a lower triangular matrix. The simplicial algorithm processes
 * the matrix column by column without forming supernodes.
 *
 * ## Applicability
 * - Symmetric positive definite sparse matrices
 * - Problems from elliptic PDEs
 * - Structural mechanics
 * - Heat transfer problems
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * 
 * Solver::SimplicialLLT solver(problem);
 * solver.solve();
 * ```
 *
 * @see SimplicialLLT for the solver implementation
 */
#ifndef RODIN_SOLVER_SIMPLICIALLLT_H
#define RODIN_SOLVER_SIMPLICIALLLT_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "LinearSolver.h"

namespace Rodin::FormLanguage
{
  template <class LinearSystem>
  struct Traits<Solver::SimplicialLLT<LinearSystem>>
  {
    using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup SimplicialLLTSpecializations SimplicialLLT Template Specializations
   * @brief Template specializations of the SimplicialLLT class.
   * @see SimplicialLLT
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for SimplicialLLT
   */
  template <class LinearSystem>
  SimplicialLLT(Variational::ProblemBase<LinearSystem>&) -> SimplicialLLT<LinearSystem>;

  /**
   * @ingroup SimplicialLLTSpecializations
   * @brief Simplicial LLT Cholesky factorization for sparse SPD matrices.
   *
   * This solver performs Cholesky decomposition @f$ A = LL^T @f$ using a simplicial
   * algorithm. It is efficient for sparse symmetric positive definite matrices
   * arising from finite element discretizations.
   *
   * **Characteristics:**
   * - **Type**: Direct solver
   * - **Matrix requirements**: Symmetric positive definite sparse
   * - **Memory**: @f$ O(nnz) @f$ for the factorization
   * - **Stability**: Requires positive definiteness
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class SimplicialLLT<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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
       * @brief Constructs a SimplicialLLT solver for the given problem.
       * @param pb Reference to the problem to solve
       */
      SimplicialLLT(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      SimplicialLLT(const SimplicialLLT& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      SimplicialLLT(SimplicialLLT&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Solves the linear system using simplicial LLT factorization.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new SimplicialLLT instance
       */
      SimplicialLLT* copy() const noexcept override
      {
        return new SimplicialLLT(*this);
      }

    private:
      /// Underlying Eigen SimplicialLLT solver
      Eigen::SimplicialLLT<OperatorType> m_solver;
  };
}

#endif



