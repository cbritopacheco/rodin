/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file IDRSTABL.h
 * @brief IDRstab(l) solver for (generally) non-symmetric linear systems.
 *
 * This header provides the IDRstab(l) solver class, an iterative method for solving:
 * @f[
 *   Ax = b
 * @f]
 * for general (possibly non-symmetric) matrices.
 *
 * IDRstab(l) is a stabilized variant in the IDR family (similar spirit to BiCGSTAB(l)).
 *
 * ## Notes
 * - This implementation uses Eigen's unsupported IDRSTABL.
 *
 * ## Usage Example
 * ```cpp
 * Solver::IDRSTABL solver(problem);
 * solver.setTolerance(1e-10).setMaxIterations(1000).setL(4).solve();
 * ```
 */
#ifndef RODIN_SOLVER_IDRSTABL_H
#define RODIN_SOLVER_IDRSTABL_H

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
  struct Traits<Solver::IDRSTABL<LinearSystem>>
  {
    /// @brief Linear system type.
      using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup IDRSTABLSpecializations IDRSTABL Template Specializations
   * @brief Template specializations of the IDRSTABL class.
   * @see IDRSTABL
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for IDRSTABL
   */
  template <class LinearSystem>
  IDRSTABL(Variational::ProblemBase<LinearSystem>&) -> IDRSTABL<LinearSystem>;

  /**
   * @ingroup IDRSTABLSpecializations
   * @brief IDRstab(l) solver for sparse systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class IDRSTABL<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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
      IDRSTABL(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /// @brief Copy constructor.
      IDRSTABL(const IDRSTABL& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      IDRSTABL(IDRSTABL&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      /// @brief Destructor.
      ~IDRSTABL() = default;

      /// @brief Sets the convergence tolerance; returns a reference to this solver.
      IDRSTABL& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /// @brief Sets the maximum number of iterations; returns a reference to this solver.
      IDRSTABL& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Sets the IDRstab(l) parameter l.
       *
       * Typical values: 2, 4, 8. Larger l may improve robustness but increases
       * work and storage.
       */
      IDRSTABL& setL(size_t l)
      {
        m_solver.setL(l);
        return *this;
      }

      /// @brief Solves the assembled linear system.
      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      /// @brief Returns whether the most recent solve converged successfully.
      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      /// @brief Returns a polymorphic copy of this solver.
      IDRSTABL* copy() const noexcept override
      {
        return new IDRSTABL(*this);
      }

    private:
      Eigen::IDRSTABL<OperatorType> m_solver;
  };

  /**
   * @ingroup IDRSTABLSpecializations
   * @brief IDRstab(l) solver for dense systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class IDRSTABL<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
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
      IDRSTABL(ProblemType& pb)
        : Parent(pb)
      {}

      /// @brief Copy constructor.
      IDRSTABL(const IDRSTABL& other)
        : Parent(other),
          m_solver(other.m_solver)
      {}

      /// @brief Move constructor.
      IDRSTABL(IDRSTABL&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      /// @brief Destructor.
      ~IDRSTABL() = default;

      /// @brief Sets the convergence tolerance; returns a reference to this solver.
      IDRSTABL& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /// @brief Sets the maximum number of iterations; returns a reference to this solver.
      IDRSTABL& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      IDRSTABL& setL(size_t l)
      {
        m_solver.setL(l);
        return *this;
      }

      /// @brief Solves the assembled linear system.
      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      /// @brief Returns whether the most recent solve converged successfully.
      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      /// @brief Returns a polymorphic copy of this solver.
      IDRSTABL* copy() const noexcept override
      {
        return new IDRSTABL(*this);
      }

    private:
      Eigen::IDRSTABL<OperatorType> m_solver;
  };
}

#endif
