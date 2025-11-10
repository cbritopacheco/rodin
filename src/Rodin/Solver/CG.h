/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_CG_H
#define RODIN_SOLVER_CG_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Rodin/Types.h"
#include "Solver.h"

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
   * The conjugate gradient (CG) method is an iterative algorithm for solving
   * linear systems @f$ Ax = b @f$ where @f$ A @f$ is symmetric positive definite.
   * It is particularly efficient for large sparse systems and provides
   * convergence in at most @f$ n @f$ iterations for an @f$ n \times n @f$ system,
   * though in practice it often converges much faster.
   *
   * The method minimizes the quadratic form @f$ \frac{1}{2}x^TAx - b^Tx @f$
   * by constructing conjugate search directions.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   *
   * **Example usage:**
   * @code{.cpp}
   * Problem problem(u, v);
   * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
   * 
   * Solver::CG solver(problem);
   * solver.setTolerance(1e-10)
   *       .setMaxIterations(1000)
   *       .solve();
   * @endcode
   *
   * @note This solver requires the system matrix to be symmetric positive definite.
   */
  template <class Scalar>
  class CG<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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
       * @brief Constructs a CG solver with default parameters.
       * @param pb Reference to the problem to solve
       *
       * Default parameters:
       * - Maximum iterations: Determined by Eigen (typically system size)
       * - Tolerance: Machine epsilon dependent
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
       * @brief Sets the convergence tolerance for the iterative solver.
       * @param tol Convergence tolerance
       * @returns Reference to this solver for method chaining
       *
       * The solver stops when @f$ \|r_k\| < \text{tol} \times \|b\| @f$ where
       * @f$ r_k @f$ is the residual at iteration @f$ k @f$.
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
       *
       * The solver will stop after this many iterations even if convergence
       * has not been achieved.
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
      /// Underlying Eigen solver implementation
      Eigen::ConjugateGradient<OperatorType, Eigen::Lower | Eigen::Upper> m_solver;
  };

  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for symmetric positive definite dense systems.
   *
   * This specialization provides the conjugate gradient method for dense
   * matrices. See the sparse specialization for detailed algorithm description.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   *
   * @see CG<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
   */
  template <class Scalar>
  class CG<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public SolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// Type of scalar values
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
      using Parent = SolverBase<LinearSystemType>;

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
       * @brief Sets the convergence tolerance for the iterative solver.
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
      /// Underlying Eigen solver implementation
      Eigen::ConjugateGradient<OperatorType, Eigen::Lower | Eigen::Upper> m_solver;
  };
}

#endif

