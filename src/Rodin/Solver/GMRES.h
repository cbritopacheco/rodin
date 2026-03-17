/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file GMRES.h
 * @brief GMRES solver for (generally) non-symmetric linear systems.
 *
 * This header provides the GMRES (Generalized Minimal Residual) solver class, an
 * iterative method for solving linear systems:
 * @f[
 *   Ax = b
 * @f]
 * for general (possibly non-symmetric) matrices.
 *
 * ## Algorithm
 * GMRES builds a Krylov subspace and computes the iterate that minimizes the
 * residual norm over that subspace.
 *
 * ## Applicability
 * - Non-symmetric or indefinite matrices
 * - Large sparse systems
 * - Convection-diffusion, saddle-point systems (when used appropriately)
 *
 * ## Notes
 * - Restarted GMRES is typically used to bound memory; Eigen's GMRES supports restart.
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 *
 * Solver::GMRES solver(problem);
 * solver.setTolerance(1e-10).setMaxIterations(1000).setRestart(50).solve();
 *
 * if (solver.success())
 *   std::cout << "Converged!\n";
 * ```
 *
 * @see GMRES for the solver implementation
 */
#ifndef RODIN_SOLVER_GMRES_H
#define RODIN_SOLVER_GMRES_H

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
  template <class LinearSystem>
  struct Traits<Solver::GMRES<LinearSystem>>
  {
    using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup GMRESSpecializations GMRES Template Specializations
   * @brief Template specializations of the GMRES class.
   * @see GMRES
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for GMRES
   */
  template <class LinearSystem>
  GMRES(Variational::ProblemBase<LinearSystem>&) -> GMRES<LinearSystem>;

  /**
   * @ingroup GMRESSpecializations
   * @brief GMRES solver for sparse systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class GMRES<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
    : public LinearSolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;
      using VectorType = Math::Vector<ScalarType>;
      using OperatorType = Math::SparseMatrix<ScalarType>;
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      GMRES(ProblemBaseType& pb)
        : Parent(pb)
      {}

      GMRES(const GMRES& other)
        : Parent(other)
      {}

      GMRES(GMRES&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      ~GMRES() = default;

      GMRES& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      GMRES& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Sets the GMRES restart parameter (dimension of Krylov subspace before restart).
       *
       * Larger values can improve convergence but increase memory and cost per iteration.
       */
      GMRES& setRestart(size_t restart)
      {
        m_solver.set_restart(restart);
        return *this;
      }

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }

    private:
      Eigen::GMRES<OperatorType> m_solver;
  };

  /**
   * @ingroup GMRESSpecializations
   * @brief GMRES solver for dense systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class GMRES<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public LinearSolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;
      using VectorType = Math::Vector<ScalarType>;
      using OperatorType = Math::Matrix<ScalarType>;
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      using ProblemType = Variational::ProblemBase<LinearSystemType>;
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      GMRES(ProblemType& pb)
        : Parent(pb)
      {}

      GMRES(const GMRES& other)
        : Parent(other)
      {}

      GMRES(GMRES&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      ~GMRES() = default;

      GMRES& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      GMRES& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      GMRES& setRestart(size_t restart)
      {
        m_solver.set_restart(restart);
        return *this;
      }

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }

    private:
      Eigen::GMRES<OperatorType> m_solver;
  };
}

#endif
