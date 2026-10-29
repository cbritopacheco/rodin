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
 * Both GMRES(problem) and GMRES(problem, pc) produce GMRES<LinearSystem>.
 * The native specialization owns a private backend whose Eigen solver retains
 * the concrete preconditioner type throughout its numerical loop. Preconditioner
 * configuration is copied at construction; it must own or safely share its
 * inputs. Solver copies clone configuration into fresh, independent backends.
 * No dispatch is added within Krylov iterations or preconditioner applications.
 * The default diagonal preconditioner and PETSc KSP specialization are unchanged.
 * A failed setup leaves the supplied solution unchanged.
 *
 * @see <a href="class_rodin_1_1_solver_1_1_g_m_r_e_s.html">GMRES</a> for the solver implementation
 */
#ifndef RODIN_SOLVER_GMRES_H
#define RODIN_SOLVER_GMRES_H

#include <limits>
#include <memory>
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
  template <class LinearSystem>
  struct Traits<Solver::GMRES<LinearSystem>>
  {
    using LinearSystemType = LinearSystem;
  };
}

/// @cond INTERNAL
namespace Rodin::Solver::Internal
{
  /**
   * @brief Owns an Eigen GMRES backend independently of its preconditioner type.
   * @tparam LinearSystem Native sparse or dense linear system.
   *
   * Settings and result inspection require no backend dispatch. One virtual call
   * performs the complete compute/solve operation; Eigen's iterations remain
   * statically typed. Configuration copies never copy the computed Eigen solver.
   */
  template <class LinearSystem>
  class EigenGMRES
  {
    private:
      using Operator = typename FormLanguage::Traits<LinearSystem>::OperatorType;
      using Scalar = typename Operator::Scalar;

      struct Settings
      {
        Optional<Real> tolerance;
        Optional<size_t> maximumIterations;
        Optional<size_t> restart;
      };

      struct Result
      {
        Eigen::ComputationInfo status;
        size_t iterations;
        Real error;
      };

      struct Backend
      {
        virtual ~Backend() = default;
        virtual std::unique_ptr<Backend> cloneConfiguration() const = 0;
        virtual Result solve(LinearSystem& system, const Settings& settings) = 0;
      };

      template <class Preconditioner>
      class EigenBackend final : public Backend
      {
        public:
          explicit EigenBackend(const Preconditioner& configuration)
            : m_configuration(configuration)
          {
            m_solver.preconditioner() = m_configuration;
          }

          std::unique_ptr<Backend> cloneConfiguration() const override
          {
            return std::make_unique<EigenBackend>(m_configuration);
          }

          Result solve(LinearSystem& system, const Settings& settings) override
          {
            if (settings.tolerance)
              m_solver.setTolerance(*settings.tolerance);
            if (settings.maximumIterations)
              m_solver.setMaxIterations(*settings.maximumIterations);
            if (settings.restart)
              m_solver.set_restart(*settings.restart);
            m_solver.compute(system.getOperator());
            if (m_solver.info() != Eigen::Success)
              return {m_solver.info(), 0, std::numeric_limits<Real>::infinity()};
            if (system.getSolution().size() == system.getVector().size())
            {
              system.getSolution() =
                m_solver.solveWithGuess(system.getVector(), system.getSolution());
            }
            else
              system.getSolution() = m_solver.solve(system.getVector());
            return {m_solver.info(), size_t(m_solver.iterations()), m_solver.error()};
          }

        private:
          const Preconditioner m_configuration;
          Eigen::GMRES<Operator, Preconditioner> m_solver;
      };

    public:
      /// @brief Constructs the unchanged Eigen default diagonal backend.
      EigenGMRES()
        : EigenGMRES(Eigen::DiagonalPreconditioner<Scalar>{})
      {}

      /// @brief Copies preconditioner configuration into a typed owned backend.
      template <class Preconditioner>
      explicit EigenGMRES(const Preconditioner& configuration)
        : m_backend(std::make_unique<EigenBackend<Preconditioner>>(configuration)),
          m_result{Eigen::InvalidInput, 0, std::numeric_limits<Real>::infinity()}
      {
        static_assert(std::is_copy_constructible_v<Preconditioner> &&
            std::is_copy_assignable_v<Preconditioner>,
          "GMRES preconditioners must support copying their configuration.");
      }

      /// @brief Copies settings and configuration, without copying factors or results.
      EigenGMRES(const EigenGMRES& other)
        : m_backend(other.m_backend->cloneConfiguration()),
          m_settings(other.m_settings),
          m_result{Eigen::InvalidInput, 0, std::numeric_limits<Real>::infinity()}
      {}

      /// @brief Moves settings into a fresh backend while leaving factors independent.
      EigenGMRES(EigenGMRES&& other)
        : EigenGMRES(static_cast<const EigenGMRES&>(other))
      {}

      /// @brief Records the backend tolerance for the next solve.
      void setTolerance(Real tolerance)
      {
        m_settings.tolerance = tolerance;
      }

      /// @brief Records the explicit iteration limit.
      void setMaxIterations(size_t count)
      {
        m_settings.maximumIterations = count;
      }

      /// @brief Records the explicit restart dimension.
      void setRestart(size_t count)
      {
        m_settings.restart = count;
      }

      /// @brief Dispatches one complete backend operation and records its diagnostics.
      void solve(LinearSystem& system)
      {
        m_result = m_backend->solve(system, m_settings);
      }

      /// @brief Returns whether the most recent operation succeeded.
      Boolean success() const
      {
        return m_result.status == Eigen::Success;
      }

      /// @brief Returns the last completed backend iteration count.
      size_t getIterationNumber() const
      {
        return m_result.iterations;
      }

      /// @brief Returns the backend estimate, not an independently measured residual.
      Real getError() const
      {
        return m_result.error;
      }

    private:
      std::unique_ptr<Backend> m_backend;
      Settings m_settings;
      Result m_result;
  };
}

/// @endcond

namespace Rodin::Solver
{
  /**
   * @defgroup GMRESSpecializations GMRES Template Specializations
   * @ingroup RodinSolvers
   * @brief Native sparse and dense GMRES implementations.
   */
  /**
   * @ingroup RodinCTAD
   * @brief Deduces the linear system when selecting the default backend.
   */
  template <class LinearSystem>
  GMRES(Variational::ProblemBase<LinearSystem>& problem) -> GMRES<LinearSystem>;

  /**
   * @ingroup RodinCTAD
   * @brief Deduces the linear system; the preconditioner type remains private.
   */
  template <class LinearSystem, class Preconditioner>
  GMRES(Variational::ProblemBase<LinearSystem>& problem, const Preconditioner& pc)
    -> GMRES<LinearSystem>;

  /**
   * @ingroup GMRESSpecializations
   * @brief Native SparseMatrix GMRES with an owned, type-erased Eigen backend.
   * @tparam Scalar Real or complex scalar field.
   */
  template <class Scalar>
  class GMRES<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
    : public LinearSolverBase<
        Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// @brief Assembled native system type.
      using LinearSystemType =
        Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;

      /// @brief Native operator type.
      using OperatorType = Math::SparseMatrix<Scalar>;

      /// @brief Native vector type.
      using VectorType = Math::Vector<Scalar>;

      /// @brief Coefficient scalar type.
      using ScalarType = Scalar;

      /// @brief Solver base type.
      using Parent = LinearSolverBase<LinearSystemType>;

      /// @brief Associated variational problem type.
      using ProblemBaseType = typename Parent::ProblemBaseType;

      /// @brief Associated problem type, retained for dense callers.
      using ProblemType = ProblemBaseType;
      using Parent::solve;

      /// @brief Constructs the existing default-preconditioned native solver.
      explicit GMRES(ProblemBaseType& problem)
        : Parent(problem)
      {}

      /**
       * @brief Copies a configured preconditioner into the private typed backend.
       * @param problem Problem whose system is solved.
       * @param pc Eigen-compatible, copy-configurable preconditioner.
       */
      template <class Preconditioner>
      GMRES(ProblemBaseType& problem, const Preconditioner& pc)
        : Parent(problem),
          m_native(pc)
      {}

      /// @brief Copies the problem and settings into fresh independent factors.
      GMRES(const GMRES& other)
        : Parent(other),
          m_native(other.m_native)
      {}

      /// @brief Transfers the problem while rebuilding an independent backend.
      GMRES(GMRES&& other)
        : Parent(std::move(other)),
          m_native(std::move(other.m_native))
      {}

      /// @brief Sets the relative backend tolerance.
      GMRES& setTolerance(const Real& tolerance)
      {
        m_native.setTolerance(tolerance);
        return *this;
      }

      /// @brief Sets the explicit iteration limit.
      GMRES& setMaxIterations(size_t count)
      {
        m_native.setMaxIterations(count);
        return *this;
      }

      /// @brief Sets the restart dimension controlling Krylov storage.
      GMRES& setRestart(size_t count)
      {
        m_native.setRestart(count);
        return *this;
      }

      /// @brief Solves the system, preserving the supplied initial-guess behavior.
      void solve(LinearSystemType& system) override
      {
        m_native.solve(system);
      }

      /// @brief Returns whether the last setup and solve succeeded.
      Boolean success() const
      {
        return m_native.success();
      }

      /// @brief Returns the last backend iteration count.
      size_t getIterationNumber() const
      {
        return m_native.getIterationNumber();
      }

      /// @brief Returns the backend error estimate, not the true residual.
      Real getError() const
      {
        return m_native.getError();
      }

      /// @brief Copies configuration into a newly allocated solver owned by the caller.
      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }

    private:
      Internal::EigenGMRES<LinearSystemType> m_native;
  };

  /**
   * @ingroup GMRESSpecializations
   * @brief Native Matrix GMRES with an owned, type-erased Eigen backend.
   * @tparam Scalar Real or complex scalar field.
   */
  template <class Scalar>
  class GMRES<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public LinearSolverBase<
        Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// @brief Assembled native system type.
      using LinearSystemType =
        Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>;

      /// @brief Native operator type.
      using OperatorType = Math::Matrix<Scalar>;

      /// @brief Native vector type.
      using VectorType = Math::Vector<Scalar>;

      /// @brief Coefficient scalar type.
      using ScalarType = Scalar;

      /// @brief Solver base type.
      using Parent = LinearSolverBase<LinearSystemType>;

      /// @brief Associated variational problem type.
      using ProblemBaseType = typename Parent::ProblemBaseType;

      /// @brief Associated problem type, retained for dense callers.
      using ProblemType = ProblemBaseType;
      using Parent::solve;

      /// @brief Constructs the existing default-preconditioned native solver.
      explicit GMRES(ProblemBaseType& problem)
        : Parent(problem)
      {}

      /**
       * @brief Copies a configured preconditioner into the private typed backend.
       * @param problem Problem whose system is solved.
       * @param pc Eigen-compatible, copy-configurable preconditioner.
       */
      template <class Preconditioner>
      GMRES(ProblemBaseType& problem, const Preconditioner& pc)
        : Parent(problem),
          m_native(pc)
      {}

      /// @brief Copies the problem and settings into fresh independent factors.
      GMRES(const GMRES& other)
        : Parent(other),
          m_native(other.m_native)
      {}

      /// @brief Transfers the problem while rebuilding an independent backend.
      GMRES(GMRES&& other)
        : Parent(std::move(other)),
          m_native(std::move(other.m_native))
      {}

      /// @brief Sets the relative backend tolerance.
      GMRES& setTolerance(const Real& tolerance)
      {
        m_native.setTolerance(tolerance);
        return *this;
      }

      /// @brief Sets the explicit iteration limit.
      GMRES& setMaxIterations(size_t count)
      {
        m_native.setMaxIterations(count);
        return *this;
      }

      /// @brief Sets the restart dimension controlling Krylov storage.
      GMRES& setRestart(size_t count)
      {
        m_native.setRestart(count);
        return *this;
      }

      /// @brief Solves the system, preserving the supplied initial-guess behavior.
      void solve(LinearSystemType& system) override
      {
        m_native.solve(system);
      }

      /// @brief Returns whether the last setup and solve succeeded.
      Boolean success() const
      {
        return m_native.success();
      }

      /// @brief Returns the last backend iteration count.
      size_t getIterationNumber() const
      {
        return m_native.getIterationNumber();
      }

      /// @brief Returns the backend error estimate, not the true residual.
      Real getError() const
      {
        return m_native.getError();
      }

      /// @brief Copies configuration into a newly allocated solver owned by the caller.
      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }

    private:
      Internal::EigenGMRES<LinearSystemType> m_native;
  };
}
#endif
