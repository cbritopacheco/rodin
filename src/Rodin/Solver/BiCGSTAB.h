/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file BiCGSTAB.h
 * @brief Bi-conjugate gradient stabilized solver for non-symmetric systems.
 *
 * This header provides the BiCGSTAB solver class, an iterative method for
 * solving non-symmetric linear systems with better convergence properties
 * than the standard Bi-Conjugate Gradient method.
 *
 * ## Algorithm
 * BiCGSTAB solves:
 * @f[
 *   Ax = b
 * @f]
 * for non-symmetric matrices @f$ A @f$. It combines BiCG with stabilization
 * to avoid the irregular convergence patterns often seen in BiCG.
 *
 * ## Applicability
 * - Non-symmetric sparse matrices
 * - Large systems where direct solvers are expensive
 * - Convection-dominated problems
 * - Advection-diffusion equations
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) + Integral(beta * Grad(u), v) - Integral(f, v);
 * 
 * Solver::BiCGSTAB solver(problem);
 * solver.setTolerance(1e-10).setMaxIterations(1000).solve();
 * ```
 *
 * @see BiCGSTAB for the solver implementation
 */
#ifndef RODIN_SOLVER_BiCGSTAB_H
#define RODIN_SOLVER_BiCGSTAB_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "LinearSolver.h"

namespace Rodin::FormLanguage
{
  template <class LinearSystem>
  struct Traits<Solver::BiCGSTAB<LinearSystem>>
  {
    using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup BiCGSTABSpecializations BiCGSTAB Template Specializations
   * @brief Template specializations of the BiCGSTAB class.
   * @see BiCGSTAB
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for BiCGSTAB
   */
  template <class LinearSystem>
  BiCGSTAB(Variational::ProblemBase<LinearSystem>&) -> BiCGSTAB<LinearSystem>;

  /**
   * @ingroup BiCGSTABSpecializations
   * @brief Bi-conjugate gradient stabilized solver for non-symmetric sparse systems.
   *
   * The BiCGSTAB (Bi-Conjugate Gradient Stabilized) method is an iterative solver
   * for non-symmetric linear systems @f$ Ax = b @f$. It combines ideas from the
   * BiCG algorithm with stabilization to avoid the irregular convergence patterns
   * often seen in the Bi-Conjugate Gradient method.
   *
   * BiCGSTAB is particularly suitable for:
   * - Non-symmetric matrices
   * - Large sparse systems
   * - Problems where direct solvers are too expensive
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   *
   * **Example usage:**
   * @code{.cpp}
   * Problem problem(u, v);
   * problem = Integral(Grad(u), Grad(v)) + Integral(beta * Grad(u), v) - Integral(f, v);
   * 
   * Solver::BiCGSTAB solver(problem);
   * solver.setTolerance(1e-10)
   *       .setMaxIterations(1000)
   *       .solve();
   * @endcode
   */
  template <class Scalar>
  class BiCGSTAB<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
    : public LinearSolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
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
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /**
       * @brief Constructs a BiCGSTAB solver with default parameters.
       * @param pb Reference to the problem to solve
       *
       * Default parameters:
       * - Maximum iterations: Determined by Eigen (typically system size)
       * - Tolerance: Machine epsilon dependent
       */
      BiCGSTAB(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      BiCGSTAB(const BiCGSTAB& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      BiCGSTAB(BiCGSTAB&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Default destructor.
       */
      ~BiCGSTAB() = default;

      /**
       * @brief Sets the convergence tolerance for the iterative solver.
       * @param tol Convergence tolerance
       * @returns Reference to this solver for method chaining
       *
       * The solver stops when @f$ \|r_k\| < \text{tol} \times \|b\| @f$ where
       * @f$ r_k @f$ is the residual at iteration @f$ k @f$.
       */
      BiCGSTAB& setTolerance(const Real& tol)
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
      BiCGSTAB& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Solves the linear system using the BiCGSTAB method.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
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
       * @returns Pointer to a new BiCGSTAB instance
       */
      BiCGSTAB* copy() const noexcept override
      {
        return new BiCGSTAB(*this);
      }

    private:
      /// Underlying Eigen solver implementation
      Eigen::BiCGSTAB<OperatorType> m_solver;
  };
}

#endif


