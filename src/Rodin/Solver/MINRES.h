/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MINRES.h
 * @brief MINRES solver for symmetric (possibly indefinite) linear systems.
 *
 * This header provides the MINRES (Minimum Residual) solver class, an iterative
 * method for solving linear systems:
 * @f[
 *   Ax = b
 * @f]
 * where @f$ A @f$ is symmetric (and may be indefinite).
 *
 * ## Algorithm
 * MINRES is a Lanczos-based Krylov method that minimizes the residual norm
 * over the generated Krylov subspace, while preserving short recurrences for
 * symmetric operators.
 *
 * ## Applicability
 * - Symmetric positive definite (works, but CG is usually preferable)
 * - Symmetric indefinite systems (e.g., some mixed formulations, saddle-point after
 *   suitable transformations, Helmholtz-type operators, etc.)
 *
 * ## Notes
 * - This implementation uses Eigen's unsupported MINRES.
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 *
 * Solver::MINRES solver(problem);
 * solver.setTolerance(1e-10).setMaxIterations(1000).solve();
 *
 * if (solver.success())
 *   std::cout << "Converged!\n";
 * ```
 *
 * @see MINRES for the solver implementation
 */
#ifndef RODIN_SOLVER_MINRES_H
#define RODIN_SOLVER_MINRES_H

#include <Eigen/Core>
#include <unsupported/Eigen/IterativeSolvers>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Matrix.h"

#include "ForwardDecls.h"
#include "Rodin/Types.h"
#include "LinearSolver.h"

namespace Rodin::FormLanguage
{
  /// @brief Form-language traits for MINRES solvers.
  template <class LinearSystem>
  struct Traits<Solver::MINRES<LinearSystem>>
  {
    /// @brief Linear system type.
      using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup MINRESSpecializations MINRES Template Specializations
   * @brief Template specializations of the MINRES class.
   * @see MINRES
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for MINRES
   */
  template <class LinearSystem>
  MINRES(Variational::ProblemBase<LinearSystem>&) -> MINRES<LinearSystem>;

  /**
   * @ingroup MINRESSpecializations
   * @brief MINRES solver for symmetric sparse systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class MINRES<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
    : public LinearSolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Scalar;
      /// @brief Vector type of the linear system.
      using VectorType = Math::Vector<ScalarType>;
      /// @brief Assembled operator type.
      using OperatorType = Math::SparseMatrix<ScalarType>;
      /// @brief Linear system type.
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      /// @brief Associated problem base type.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      /// @brief Parent class type.
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /// @brief Constructs the solver from the problem to be solved.
      MINRES(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /// @brief Copy constructor.
      MINRES(const MINRES& other)
        : Parent(other),
          m_solver(other.m_solver)
      {}

      /// @brief Move constructor.
      MINRES(MINRES&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      /// @brief Destructor.
      ~MINRES() = default;

      /// @brief Sets the convergence tolerance; returns a reference to this solver.
      MINRES& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /// @brief Sets the maximum number of iterations; returns a reference to this solver.
      MINRES& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /// @brief Solves the assembled linear system.
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        if (axb.getSolution().size() == axb.getVector().size())
          axb.getSolution() = m_solver.solveWithGuess(axb.getVector(), axb.getSolution());
        else
          axb.getSolution() = m_solver.solve(axb.getVector());
      }

      /// @brief Returns whether the most recent solve converged successfully.
      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      /// @brief Returns a polymorphic copy of this solver.
      MINRES* copy() const noexcept override
      {
        return new MINRES(*this);
      }

    private:
      Eigen::MINRES<OperatorType> m_solver;
  };

  /**
   * @ingroup MINRESSpecializations
   * @brief MINRES solver for symmetric dense systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class MINRES<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public LinearSolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Scalar;
      /// @brief Vector type of the linear system.
      using VectorType = Math::Vector<ScalarType>;
      /// @brief Assembled operator type.
      using OperatorType = Math::Matrix<ScalarType>;
      /// @brief Linear system type.
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      /// @brief Problem type solved by this solver.
      using ProblemType = Variational::ProblemBase<LinearSystemType>;
      /// @brief Parent class type.
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /// @brief Constructs the solver from the problem to be solved.
      MINRES(ProblemType& pb)
        : Parent(pb)
      {}

      /// @brief Copy constructor.
      MINRES(const MINRES& other)
        : Parent(other),
          m_solver(other.m_solver)
      {}

      /// @brief Move constructor.
      MINRES(MINRES&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      /// @brief Destructor.
      ~MINRES() = default;

      /// @brief Sets the convergence tolerance; returns a reference to this solver.
      MINRES& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /// @brief Sets the maximum number of iterations; returns a reference to this solver.
      MINRES& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /// @brief Solves the assembled linear system.
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        if (axb.getSolution().size() == axb.getVector().size())
          axb.getSolution() = m_solver.solveWithGuess(axb.getVector(), axb.getSolution());
        else
          axb.getSolution() = m_solver.solve(axb.getVector());
      }

      /// @brief Returns whether the most recent solve converged successfully.
      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      /// @brief Returns a polymorphic copy of this solver.
      MINRES* copy() const noexcept override
      {
        return new MINRES(*this);
      }

    private:
      Eigen::MINRES<OperatorType> m_solver;
  };
}

#endif
