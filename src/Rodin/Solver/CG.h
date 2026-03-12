/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CG.h
 * @brief Conjugate gradient solver for symmetric positive definite systems.
 *
 * This header provides the CG (Conjugate Gradient) solver class, an iterative
 * method for solving linear systems with symmetric positive definite matrices.
 *
 * ## Algorithm
 * The conjugate gradient method solves:
 * @f[
 *   Ax = b
 * @f]
 * where @f$ A @f$ is symmetric positive definite. It minimizes the quadratic
 * form @f$ \frac{1}{2}x^TAx - b^Tx @f$ using conjugate search directions.
 *
 * ## Convergence
 * CG converges in at most @f$ n @f$ iterations in exact arithmetic. The error
 * decreases as:
 * @f[
 *   \|x_k - x^*\|_A \leq 2 \left( \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right)^k \|x_0 - x^*\|_A
 * @f]
 * where @f$ \kappa @f$ is the condition number of @f$ A @f$.
 *
 * ## Applicability
 * - Symmetric positive definite matrices
 * - Large sparse systems
 * - Elliptic PDEs (Poisson, heat equation, elasticity)
 * - Systems amenable to preconditioning
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * 
 * Solver::CG solver(problem);
 * solver.setTolerance(1e-10).setMaxIterations(1000).solve();
 * 
 * if (solver.success())
 *   std::cout << "Converged!\n";
 * ```
 *
 * @see CG for the solver implementation
 */
#ifndef RODIN_SOLVER_CG_H
#define RODIN_SOLVER_CG_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Rodin/Types.h"
#include "LinearSolver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup CGSpecializations CG Template Specializations
   * @brief Template specializations of the CG class.
   * @see CG
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for CG
   */
  template <class LinearSystem>
  CG(Variational::ProblemBase<LinearSystem>&) -> CG<LinearSystem>;

  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for symmetric positive definite sparse systems.
   *
   * The conjugate gradient method is an iterative solver for linear systems
   * @f$ Ax = b @f$ where @f$ A @f$ is symmetric positive definite. It constructs
   * a sequence of approximations by minimizing the error in the @f$ A @f$-norm.
   *
   * **Characteristics:**
   * - **Type**: Iterative solver
   * - **Matrix requirements**: Symmetric positive definite
   * - **Memory**: @f$ O(n) @f$ for vectors
   * - **Convergence**: Depends on condition number @f$ \kappa(A) @f$
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class CG<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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
       * @brief Constructs a CG solver with default parameters.
       * @param pb Reference to the problem to solve
       */
      CG(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      CG(const CG& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      CG(CG&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Default destructor.
       */
      ~CG() = default;

      /**
       * @brief Sets the convergence tolerance.
       * @param tol Convergence tolerance
       * @returns Reference to this solver for method chaining
       *
       * The solver stops when @f$ \|r_k\| < \text{tol} \times \|b\| @f$.
       */
      CG& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /**
       * @brief Sets the maximum number of iterations.
       * @param maxIt Maximum number of iterations
       * @returns Reference to this solver for method chaining
       */
      CG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Solves the linear system using the conjugate gradient method.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      /**
       * @brief Checks if the solver converged successfully.
       * @returns true if the solver converged, false otherwise
       */
      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new CG instance
       */
      CG* copy() const noexcept override
      {
        return new CG(*this);
      }

    private:
      /// Underlying Eigen conjugate gradient solver
      Eigen::ConjugateGradient<OperatorType, Eigen::Lower | Eigen::Upper> m_solver;
  };

  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for symmetric positive definite dense systems.
   *
   * This specialization provides the conjugate gradient method for dense matrices.
   * See the sparse specialization for detailed algorithm description.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class CG<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
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
      using ProblemType = Variational::ProblemBase<LinearSystemType>;

      /// Parent class type
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /**
       * @brief Constructs a CG solver with default parameters.
       * @param pb Reference to the problem to solve
       */
      CG(ProblemType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      CG(const CG& other)
        : Parent(other)
      {}

      /**
       * @brief Default destructor.
       */
      ~CG() = default;

      /**
       * @brief Sets the convergence tolerance.
       * @param tol Convergence tolerance
       * @returns Reference to this solver for method chaining
       */
      CG& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /**
       * @brief Sets the maximum number of iterations.
       * @param maxIt Maximum number of iterations
       * @returns Reference to this solver for method chaining
       */
      CG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Solves the linear system using the conjugate gradient method.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solveWithGuess(axb.getVector(), axb.getSolution());
      }

      /**
       * @brief Checks if the solver converged successfully.
       * @returns true if the solver converged, false otherwise
       */
      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new CG instance
       */
      CG* copy() const noexcept override
      {
        return new CG(*this);
      }

    private:
      /// Underlying Eigen conjugate gradient solver
      Eigen::ConjugateGradient<OperatorType, Eigen::Lower | Eigen::Upper> m_solver;
  };
}

#endif

