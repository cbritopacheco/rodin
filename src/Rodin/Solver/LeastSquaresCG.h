/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LeastSquaresCG.h
 * @brief Least-squares conjugate gradient solver.
 *
 * This header provides the LeastSquaresCG solver class for solving least-squares
 * problems using the conjugate gradient method applied to the normal equations.
 *
 * ## Algorithm
 * Solves the least-squares problem:
 * @f[
 *   \min_x \|Ax - b\|^2
 * @f]
 * by applying CG to the normal equations @f$ A^TAx = A^Tb @f$.
 *
 * ## Applicability
 * - Overdetermined systems (more equations than unknowns)
 * - Rectangular matrices
 * - Least-squares fitting problems
 * - Problems where @f$ A^TA @f$ is well-conditioned
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * 
 * Solver::LeastSquaresCG solver(problem);
 * solver.setTolerance(1e-10).setMaxIterations(1000).solve();
 * ```
 *
 * @see LeastSquaresCG for the solver implementation
 */
#ifndef RODIN_SOLVER_LEASTSQUARESCG_H
#define RODIN_SOLVER_LEASTSQUARESCG_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup LeastSquaresCGSpecializations LeastSquaresCG Template Specializations
   * @brief Template specializations of the LeastSquaresCG class.
   * @see LeastSquaresCG
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for LeastSquaresCG
   */
  template <class LinearSystem>
  LeastSquaresCG(Variational::ProblemBase<LinearSystem>&) -> LeastSquaresCG<LinearSystem>;

  /**
   * @ingroup LeastSquaresCGSpecializations
   * @brief Least-squares conjugate gradient solver for sparse systems.
   *
   * The Least-Squares Conjugate Gradient method solves the least-squares problem
   * @f$ \min_x \|Ax - b\|^2 @f$ by applying the conjugate gradient method to the
   * normal equations @f$ A^TAx = A^Tb @f$. This is particularly useful for:
   * - Overdetermined systems (more equations than unknowns)
   * - Rectangular matrices
   * - Least-squares fitting problems
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   *
   * **Example usage:**
   * @code{.cpp}
   * Problem problem(u, v);
   * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
   * 
   * Solver::LeastSquaresCG solver(problem);
   * solver.setTolerance(1e-10)
   *       .setMaxIterations(1000)
   *       .solve();
   * @endcode
   *
   * @note For well-conditioned square systems, prefer CG which is more efficient.
   */
  template <class Scalar>
  class LeastSquaresCG<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
    : public SolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// Type of scalar values
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
       * @brief Constructs a LeastSquaresCG solver with default parameters.
       * @param pb Reference to the problem to solve
       */
      LeastSquaresCG(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      LeastSquaresCG(const LeastSquaresCG& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      LeastSquaresCG(LeastSquaresCG&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Default destructor.
       */
      ~LeastSquaresCG() = default;

      /**
       * @brief Sets the convergence tolerance for the iterative solver.
       * @param tol Convergence tolerance
       * @returns Reference to this solver for method chaining
       *
       * The solver stops when @f$ \|r_k\| < \text{tol} \times \|b\| @f$ where
       * @f$ r_k @f$ is the residual at iteration @f$ k @f$.
       */
      LeastSquaresCG& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /**
       * @brief Sets the maximum number of iterations.
       * @param maxIt Maximum number of iterations
       * @returns Reference to this solver for method chaining
       *
       * The solver will stop after this many iterations even if convergence
       * has not been achieved.
       */
      LeastSquaresCG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Solves the least-squares problem using the conjugate gradient method.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new LeastSquaresCG instance
       */
      LeastSquaresCG* copy() const noexcept override
      {
        return new LeastSquaresCG(*this);
      }

    private:
      /// Underlying Eigen solver implementation
      Eigen::LeastSquaresConjugateGradient<OperatorType> m_solver;
  };

  /**
   * @ingroup LeastSquaresCGSpecializations
   * @brief Least-squares conjugate gradient solver for dense systems.
   *
   * This specialization provides the least-squares conjugate gradient method
   * for dense matrices. See the sparse specialization for detailed algorithm
   * description.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   *
   * @see LeastSquaresCG<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
   */
  template <class Scalar>
  class LeastSquaresCG<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public SolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// Type of scalar values
      using ScalarType = Scalar;

      /// Type of solution and right-hand side vectors
      using VectorType = Math::Vector<Scalar>;

      /// Type of system matrix (operator)
      using OperatorType = Math::Matrix<Scalar>;

      /// Type of linear system being solved
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      /// Type of problem being solved
      using ProblemType = Variational::ProblemBase<LinearSystemType>;

      /// Parent class type
      using Parent = SolverBase<LinearSystemType>;

      using Parent::solve;

      /**
       * @brief Constructs a LeastSquaresCG solver with default parameters.
       * @param pb Reference to the problem to solve
       */
      LeastSquaresCG(ProblemType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      LeastSquaresCG(const LeastSquaresCG& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      LeastSquaresCG(LeastSquaresCG&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Default destructor.
       */
      ~LeastSquaresCG() = default;

      /**
       * @brief Sets the convergence tolerance for the iterative solver.
       * @param tol Convergence tolerance
       * @returns Reference to this solver for method chaining
       */
      LeastSquaresCG& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /**
       * @brief Sets the maximum number of iterations.
       * @param maxIt Maximum number of iterations
       * @returns Reference to this solver for method chaining
       */
      LeastSquaresCG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Solves the least-squares problem using the conjugate gradient method.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getOperator());
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new LeastSquaresCG instance
       */
      LeastSquaresCG* copy() const noexcept override
      {
        return new LeastSquaresCG(*this);
      }

    private:
      /// Underlying Eigen solver implementation
      Eigen::LeastSquaresConjugateGradient<OperatorType> m_solver;
  };
}

#endif


