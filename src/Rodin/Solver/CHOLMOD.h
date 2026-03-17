/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CHOLMOD.h
 * @brief CHOLMOD supernodal Cholesky factorization solver wrapper.
 *
 * This header provides a wrapper for CHOLMOD from SuiteSparse, a high-performance
 * solver for sparse symmetric positive definite systems using supernodal LLT
 * Cholesky factorization.
 *
 * ## Algorithm
 * CHOLMOD computes:
 * @f[
 *   A = LL^T
 * @f]
 * where @f$ L @f$ is lower triangular. The supernodal method groups columns
 * with similar sparsity patterns to exploit dense matrix operations (BLAS)
 * for improved performance.
 *
 * ## Applicability
 * - Symmetric positive definite sparse matrices
 * - Large-scale finite element problems
 * - Structural mechanics
 * - Heat transfer problems
 * - Graph Laplacian systems
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * 
 * Solver::CHOLMOD::SupernodalLLT solver(problem);
 * solver.solve();
 * ```
 *
 * @note This solver requires CHOLMOD from SuiteSparse to be installed
 * and RODIN_USE_CHOLMOD to be defined at compile time.
 *
 * @see CHOLMOD::SupernodalLLT for the solver implementation
 */
#ifndef RODIN_SOLVER_CHOLMOD_H
#define RODIN_SOLVER_CHOLMOD_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_CHOLMOD

#include <optional>
#include <functional>

#include <Eigen/CholmodSupport>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "LinearSolver.h"

namespace Rodin::FormLanguage
{
  template <class LinearSystem>
  struct Traits<Solver::CHOLMOD::SupernodalLLT<LinearSystem>>
  {
    using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver::CHOLMOD
{
  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for SupernodalLLT
   */
  template <class LinearSystemType>
  SupernodalLLT(Variational::ProblemBase<LinearSystemType>&) -> SupernodalLLT<LinearSystemType>;

  /**
   * @brief CHOLMOD supernodal LLT Cholesky factorization solver.
   *
   * This solver provides access to CHOLMOD's supernodal LLT factorization,
   * a high-performance direct solver from SuiteSparse for symmetric positive
   * definite sparse matrices. The supernodal approach exploits dense substructures
   * for improved performance using BLAS operations.
   *
   * **Characteristics:**
   * - **Type**: Direct solver
   * - **Matrix requirements**: Symmetric positive definite sparse
   * - **Memory**: Efficient supernodal representation
   * - **Performance**: Excellent for large SPD systems
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class SupernodalLLT<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
    : public LinearSolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>>
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
       * @brief Constructs a CHOLMOD SupernodalLLT solver for the given problem.
       * @param pb Reference to the problem to solve
       */
      SupernodalLLT(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      SupernodalLLT(const SupernodalLLT& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      SupernodalLLT(SupernodalLLT&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Default destructor.
       */
      ~SupernodalLLT() = default;

      /**
       * @brief Solves the linear system using CHOLMOD supernodal LLT.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      /**
       * @brief Checks if the factorization succeeded.
       * @returns true if the solver succeeded, false otherwise
       */
      Real success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new SupernodalLLT instance
       */
      SupernodalLLT* copy() const noexcept override
      {
        return new SupernodalLLT(*this);
      }

    private:
      /// Underlying Eigen CHOLMOD supernodal LLT solver
      Eigen::CholmodSupernodalLLT<OperatorType> m_solver;
  };
}

#endif // #ifdef RODIN_USE_CHOLMOD
#endif



