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
 * # Algorithm
 * GMRES builds a Krylov subspace and computes the iterate that minimizes the
 * residual norm over that subspace.
 *
 * # Applicability
 * - Non-symmetric or indefinite matrices
 * - Large sparse systems
 * - Convection-diffusion, saddle-point systems (when used appropriately)
 *
 * # Notes
 * - Restarted GMRES is typically used to bound memory; Eigen's GMRES supports restart.
 *
 * # Usage Example
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
 * @par Custom native preconditioners
 * `GMRES<LinearSystem, Preconditioner>` accepts an Eigen-compatible
 * preconditioner and exposes its owned configuration through getPreconditioner().
 * The default void parameter preserves the existing native diagonal
 * preconditioner and PETSc KSP specialization. Custom types must support copy
 * assignment of valid configuration. Types containing noncopyable numerical
 * factors should copy settings and immutable owned inputs into a fresh state;
 * compute rebuilds factors independently for every solver instance and solve.
 * A failed preconditioner setup leaves the supplied solution unchanged.
 *
 * @see <a href="class_rodin_1_1_solver_1_1_g_m_r_e_s.html">GMRES</a> for the solver implementation
 */
#ifndef RODIN_SOLVER_GMRES_H
#define RODIN_SOLVER_GMRES_H

#include <type_traits>
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
  /// @brief Form-language traits for GMRES solvers.
  template <class LinearSystem, class Preconditioner>
  struct Traits<Solver::GMRES<LinearSystem, Preconditioner>>
  {
    /// @brief Linear system type.
    using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup GMRESSpecializations GMRES Template Specializations
   * @brief Template specializations of the GMRES class.
   * @see @ref GMRES
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref GMRES "GMRES<LinearSystem<SparseMatrix<Scalar>, Vector<Scalar>>>" | Iterative GMRES solver for sparse systems. |
   * | @ref GMRES "GMRES<LinearSystem<Matrix<Scalar>, Vector<Scalar>>>" | Iterative GMRES solver for dense systems. |
   * | @ref GMRES "GMRES<PETSc::LinearSystem>" | PETSc-backed GMRES solver for distributed systems. |
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for GMRES
   * @param pb Variational problem to operate on.
   */
  template <class LinearSystem>
  GMRES(Variational::ProblemBase<LinearSystem>& pb) -> GMRES<LinearSystem>;

  /**
   * @ingroup GMRESSpecializations
   * @brief GMRES solver for sparse systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar, class Preconditioner>
  class GMRES<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>,
    Preconditioner>
    final : public LinearSolverBase<
        Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Scalar;

      /// @brief Eigen-compatible preconditioner; void retains the backend default.
      using PreconditionerType = std::conditional_t<std::is_void_v<Preconditioner>,
        Eigen::DiagonalPreconditioner<ScalarType>, Preconditioner>;

      static_assert(std::is_copy_assignable_v<PreconditionerType>,
        "GMRES preconditioners must support copying their configuration.");

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

      /**
       * @brief Constructs the solver from the problem to be solved.
       * @param pb Variational problem to operate on.
       */
      GMRES(ProblemBaseType& pb)
        : Parent(pb),
          m_maxIterations(),
          m_restart()
      {}

      /**
       * @brief Copies the problem and configured settings into a fresh backend.
       * @param other Object to copy from.
       */
      GMRES(const GMRES& other)
        : Parent(other),
          m_maxIterations(other.m_maxIterations),
          m_restart(other.m_restart)
      {
        if constexpr (!std::is_void_v<Preconditioner>)
          m_solver.preconditioner() = other.m_solver.preconditioner();
        setTolerance(other.m_solver.tolerance());
        if (m_maxIterations)
          setMaxIterations(*m_maxIterations);
        if (m_restart)
          setRestart(*m_restart);
      }

      /**
       * @brief Transfers the problem and configured settings into a fresh backend.
       * @param other Object to move from.
       */
      GMRES(GMRES&& other)
        : Parent(std::move(other)),
          m_maxIterations(other.m_maxIterations),
          m_restart(other.m_restart)
      {
        if constexpr (!std::is_void_v<Preconditioner>)
          m_solver.preconditioner() = other.m_solver.preconditioner();
        setTolerance(other.m_solver.tolerance());
        if (m_maxIterations)
          setMaxIterations(*m_maxIterations);
        if (m_restart)
          setRestart(*m_restart);
      }

      /// @brief Destructor.
      ~GMRES() = default;

      /**
       * @brief Accesses the preconditioner configuration before the next solve.
       * @returns Mutable preconditioner owned by this solver.
       *
       * Custom preconditioners follow Eigen's compute/solve/info protocol.
       * Copy assignment must preserve configuration with valid ownership;
       * expensive backend factors should be rebuilt by compute, not shared
       * between solver copies. Every solve recomputes the preconditioner.
       */
      PreconditionerType& getPreconditioner()
      {
        return m_solver.preconditioner();
      }

      /// @brief Returns the owned preconditioner for inspection.
      const PreconditionerType& getPreconditioner() const
      {
        return m_solver.preconditioner();
      }

      /**
       * @brief Sets the convergence tolerance; returns a reference to this solver.
       * @param tol Convergence tolerance.
       * @returns Reference to this object after the operation.
       */
      GMRES& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /**
       * @brief Sets the maximum number of iterations; returns a reference to this solver.
       * @param maxIt Maximum number of solver iterations.
       * @returns Reference to this object after the operation.
       */
      GMRES& setMaxIterations(size_t maxIt)
      {
        m_maxIterations = maxIt;
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Sets the GMRES restart parameter (dimension of Krylov subspace before restart).
       *
       * Larger values can improve convergence but increase memory and cost per iteration.
       * @param restart Number of iterations between solver restarts.
       * @returns Reference to this object after the operation.
       */
      GMRES& setRestart(size_t restart)
      {
        m_restart = restart;
        m_solver.set_restart(restart);
        return *this;
      }

      /**
       * @brief Solves the assembled linear system.
       * @param axb Variational expression defining the problem.
       */
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        if (m_solver.info() != Eigen::Success)
          return;
        if (axb.getSolution().size() == axb.getVector().size())
          axb.getSolution() = m_solver.solveWithGuess(axb.getVector(), axb.getSolution());
        else
          axb.getSolution() = m_solver.solve(axb.getVector());
      }

      /**
       * @brief Returns whether the most recent solve converged successfully.
       * @returns Whether the operation completed successfully.
       */
      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      /**
       * @brief Returns Eigen's iteration count for the most recent solve.
       * @returns Solver-reported iteration count, without reinterpretation.
       */
      size_t getIterationNumber() const
      {
        return m_solver.iterations();
      }

      /**
       * @brief Returns Eigen's relative error estimate for the most recent solve.
       *
       * This backend estimate is not an independently evaluated true residual.
       * Both diagnostics require a completed solve.
       *
       * @returns Solver-reported relative error estimate.
       */
      Real getError() const
      {
        return m_solver.error();
      }

      /**
       * @brief Returns a polymorphic copy of this solver.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }

    private:
      Optional<size_t>
        m_maxIterations; ///< Explicit iteration limit; absence retains Eigen's size-dependent default.
      Optional<size_t> m_restart; ///< Explicit Krylov restart dimension.
      Eigen::GMRES<OperatorType, PreconditionerType> m_solver;
  };

  /**
   * @ingroup GMRESSpecializations
   * @brief GMRES solver for dense systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar, class Preconditioner>
  class GMRES<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>,
    Preconditioner>
    final : public LinearSolverBase<
        Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Scalar;

      /// @brief Eigen-compatible preconditioner; void retains the backend default.
      using PreconditionerType = std::conditional_t<std::is_void_v<Preconditioner>,
        Eigen::DiagonalPreconditioner<ScalarType>, Preconditioner>;

      static_assert(std::is_copy_assignable_v<PreconditionerType>,
        "GMRES preconditioners must support copying their configuration.");

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

      /**
       * @brief Constructs the solver from the problem to be solved.
       * @param pb Variational problem to operate on.
       */
      GMRES(ProblemType& pb)
        : Parent(pb),
          m_maxIterations(),
          m_restart()
      {}

      /**
       * @brief Copies the problem and configured settings into a fresh backend.
       * @param other Object to copy from.
       */
      GMRES(const GMRES& other)
        : Parent(other),
          m_maxIterations(other.m_maxIterations),
          m_restart(other.m_restart)
      {
        if constexpr (!std::is_void_v<Preconditioner>)
          m_solver.preconditioner() = other.m_solver.preconditioner();
        setTolerance(other.m_solver.tolerance());
        if (m_maxIterations)
          setMaxIterations(*m_maxIterations);
        if (m_restart)
          setRestart(*m_restart);
      }

      /**
       * @brief Transfers the problem and configured settings into a fresh backend.
       * @param other Object to move from.
       */
      GMRES(GMRES&& other)
        : Parent(std::move(other)),
          m_maxIterations(other.m_maxIterations),
          m_restart(other.m_restart)
      {
        if constexpr (!std::is_void_v<Preconditioner>)
          m_solver.preconditioner() = other.m_solver.preconditioner();
        setTolerance(other.m_solver.tolerance());
        if (m_maxIterations)
          setMaxIterations(*m_maxIterations);
        if (m_restart)
          setRestart(*m_restart);
      }

      /// @brief Destructor.
      ~GMRES() = default;

      /**
       * @brief Accesses the preconditioner configuration before the next solve.
       * @returns Mutable preconditioner owned by this solver.
       *
       * Custom preconditioners follow Eigen's compute/solve/info protocol.
       * Copy assignment must preserve configuration with valid ownership;
       * expensive backend factors should be rebuilt by compute, not shared
       * between solver copies. Every solve recomputes the preconditioner.
       */
      PreconditionerType& getPreconditioner()
      {
        return m_solver.preconditioner();
      }

      /// @brief Returns the owned preconditioner for inspection.
      const PreconditionerType& getPreconditioner() const
      {
        return m_solver.preconditioner();
      }

      /**
       * @brief Sets the convergence tolerance; returns a reference to this solver.
       * @param tol Convergence tolerance.
       * @returns Reference to this object after the operation.
       */
      GMRES& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      /**
       * @brief Sets the maximum number of iterations; returns a reference to this solver.
       * @param maxIt Maximum number of solver iterations.
       * @returns Reference to this object after the operation.
       */
      GMRES& setMaxIterations(size_t maxIt)
      {
        m_maxIterations = maxIt;
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Sets the Krylov restart dimension; returns a reference to this solver.
       * @param restart Number of iterations between solver restarts.
       * @returns Reference to this object after the operation.
       */
      GMRES& setRestart(size_t restart)
      {
        m_restart = restart;
        m_solver.set_restart(restart);
        return *this;
      }

      /**
       * @brief Solves the assembled linear system.
       * @param axb Variational expression defining the problem.
       */
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        if (m_solver.info() != Eigen::Success)
          return;
        if (axb.getSolution().size() == axb.getVector().size())
          axb.getSolution() = m_solver.solveWithGuess(axb.getVector(), axb.getSolution());
        else
          axb.getSolution() = m_solver.solve(axb.getVector());
      }

      /**
       * @brief Returns whether the most recent solve converged successfully.
       * @returns Whether the operation completed successfully.
       */
      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      /**
       * @brief Returns Eigen's iteration count for the most recent solve.
       * @returns Solver-reported iteration count, without reinterpretation.
       */
      size_t getIterationNumber() const
      {
        return m_solver.iterations();
      }

      /**
       * @brief Returns Eigen's relative error estimate for the most recent solve.
       *
       * This backend estimate is not an independently evaluated true residual.
       * Both diagnostics require a completed solve.
       *
       * @returns Solver-reported relative error estimate.
       */
      Real getError() const
      {
        return m_solver.error();
      }

      /**
       * @brief Returns a polymorphic copy of this solver.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }

    private:
      Optional<size_t>
        m_maxIterations; ///< Explicit iteration limit; absence retains Eigen's size-dependent default.
      Optional<size_t> m_restart; ///< Explicit Krylov restart dimension.
      Eigen::GMRES<OperatorType, PreconditionerType> m_solver;
  };
}

#endif
