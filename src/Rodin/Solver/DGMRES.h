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
  template <class LinearSystem>
  struct Traits<Solver::DGMRES<LinearSystem>>
  {
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
      using ScalarType = Scalar;
      using VectorType = Math::Vector<ScalarType>;
      using OperatorType = Math::SparseMatrix<ScalarType>;
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      DGMRES(ProblemBaseType& pb)
        : Parent(pb)
      {}

      DGMRES(const DGMRES& other)
        : Parent(other)
      {}

      DGMRES(DGMRES&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      ~DGMRES() = default;

      DGMRES& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

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

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

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
      using ScalarType = Scalar;
      using VectorType = Math::Vector<ScalarType>;
      using OperatorType = Math::Matrix<ScalarType>;
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      using ProblemType = Variational::ProblemBase<LinearSystemType>;
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      DGMRES(ProblemType& pb)
        : Parent(pb)
      {}

      DGMRES(const DGMRES& other)
        : Parent(other),
          m_solver(other.m_solver)
      {}

      DGMRES(DGMRES&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      ~DGMRES() = default;

      DGMRES& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      DGMRES& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      DGMRES& setRestart(size_t restart)
      {
        m_solver.set_restart(restart);
        return *this;
      }

      DGMRES& setDeflationSize(size_t d)
      {
        m_solver.set_d(d);
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

      DGMRES* copy() const noexcept override
      {
        return new DGMRES(*this);
      }

    private:
      Eigen::DGMRES<OperatorType> m_solver;
  };
}

#endif
