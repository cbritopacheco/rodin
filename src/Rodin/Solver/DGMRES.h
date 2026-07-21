/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file DGMRES.h
 * @brief DGMRES solver for (generally) non-symmetric linear systems.
 *
 * This header provides the DGMRES (Deflated GMRES) solver class, an iterative
 * method for solving linear systems:
 * @f[
 *   Ax = b
 * @f]
 * for general (possibly non-symmetric) matrices.
 *
 * DGMRES accelerates restarted GMRES by deflating (approximations of) invariant
 * subspaces associated with slow-to-converge eigencomponents.
 *
 * ## Notes
 * - This implementation uses Eigen's unsupported DGMRES.
 *
 * ## Usage Example
 * ```cpp
 * Solver::DGMRES solver(problem);
 * solver.setTolerance(1e-10).setMaxIterations(1000).setRestart(50).solve();
 * ```
 */
#ifndef RODIN_SOLVER_DGMRES_H
#define RODIN_SOLVER_DGMRES_H

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
  /// @brief Form-language traits for DGMRES solvers.
  template <class LinearSystem>
  struct Traits<Solver::DGMRES<LinearSystem>>
  {
    /// @brief Linear system type.
      using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup DGMRESSpecializations DGMRES Template Specializations
   * @brief Template specializations of the DGMRES class.
   * @see DGMRES
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for DGMRES
   */
  template <class LinearSystem>
  DGMRES(Variational::ProblemBase<LinearSystem>&) -> DGMRES<LinearSystem>;

  /**
   * @ingroup DGMRESSpecializations
   * @brief DGMRES solver for sparse systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class DGMRES<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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
      DGMRES(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /// @brief Copy constructor.
      DGMRES(const DGMRES& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      DGMRES(DGMRES&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      /// @brief Destructor.
      ~DGMRES() = default;

      /// @brief Sets the convergence tolerance; returns a reference to this solver.
      DGMRES& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /// @brief Sets the maximum number of iterations; returns a reference to this solver.
      DGMRES& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Sets the restart parameter (dimension of Krylov subspace before restart).
       */
      DGMRES& setRestart(size_t restart)
      {
        m_solver.set_restart(restart);
        return *this;
      }

      /**
       * @brief Sets the number of vectors used for deflation (Eigen calls this "d").
       *
       * Larger values may improve convergence but increase overhead.
       */
      DGMRES& setEigenv(size_t d)
      {
        m_solver.setEigenv(d);
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
      DGMRES* copy() const noexcept override
      {
        return new DGMRES(*this);
      }

    private:
      Eigen::DGMRES<OperatorType> m_solver;
  };

  /**
   * @ingroup DGMRESSpecializations
   * @brief DGMRES solver for dense systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class DGMRES<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
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
      DGMRES(ProblemType& pb)
        : Parent(pb)
      {}

      /// @brief Copy constructor.
      DGMRES(const DGMRES& other)
        : Parent(other),
          m_solver(other.m_solver)
      {}

      /// @brief Move constructor.
      DGMRES(DGMRES&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      /// @brief Destructor.
      ~DGMRES() = default;

      /// @brief Sets the convergence tolerance; returns a reference to this solver.
      DGMRES& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /// @brief Sets the maximum number of iterations; returns a reference to this solver.
      DGMRES& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /// @brief Sets the Krylov restart dimension; returns a reference to this solver.
      DGMRES& setRestart(size_t restart)
      {
        m_solver.set_restart(restart);
        return *this;
      }

      /// @brief Sets the number of deflation vectors; returns a reference to this solver.
      DGMRES& setDeflationSize(size_t d)
      {
        m_solver.set_d(d);
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
      DGMRES* copy() const noexcept override
      {
        return new DGMRES(*this);
      }

    private:
      Eigen::DGMRES<OperatorType> m_solver;
  };
}

#endif
