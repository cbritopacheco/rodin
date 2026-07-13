/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_NEWTONSOLVER_H
#define RODIN_SOLVER_NEWTONSOLVER_H

#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Copyable.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Solver
{
  /**
   * @brief Base interface for Newton-type nonlinear solvers.
   *
   * Newton solver class for solving nonlinear systems through a
   * Variational::ProblemBase that assembles the tangent system at each iterate:
   * @f[
   * J(u^k)\,\delta u^k = -F(u^k), \qquad
   * u^{k + 1} = u^k + \delta u^k.
   * @f]
   *
   * The linear solver is passed by reference at construction time.
   * The associated ProblemBase is obtained from the solver via
   * LinearSolverBase::getProblem().
   *
   * @tparam LinearSolver Type of the linear solver used at each Newton step.
   *   Must have a FormLanguage::Traits specialization providing LinearSystemType.
   */
  template <class LinearSolver>
  class NewtonSolverBase : public Copyable
  {
    public:
      using LinearSolverType = LinearSolver;
      /// @brief Linear system type.
      using LinearSystemType = typename FormLanguage::Traits<LinearSolver>::LinearSystemType;
      /// @brief Associated problem base type.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      /// @brief Solution vector type.
      using SolutionType = typename FormLanguage::Traits<LinearSystemType>::VectorType;

      /// @brief Parent class type.
      using Parent = Copyable;

      virtual ~NewtonSolverBase() = default;

      /**
       * @brief Constructs the base from a linear solver reference.
       * @param solver The linear solver whose associated ProblemBase will
       *   be used for tangent assembly at each Newton iteration.
       */
      explicit NewtonSolverBase(LinearSolver& solver)
        : m_solver(solver)
      {}

      /**
       * @brief Returns the associated linear solver.
       */
      const LinearSolver& getLinearSolver() const noexcept
      {
        return m_solver.get();
      }

      /**
       * @brief Returns the associated linear solver.
       */
      LinearSolver& getLinearSolver() noexcept
      {
        return m_solver.get();
      }

      /**
       * @brief Solve a nonlinear system starting from the initial guess stored in @p x.
       * @param[in,out] x On entry: initial guess. On exit: final Newton iterate.
       */
      virtual void solve(SolutionType& x) = 0;

      /**
       * @brief Solve a nonlinear system using a GridFunction as the solution container.
       *
       * The GridFunction's internal data vector is used both as the initial
       * guess and as the storage for the converged solution.  The backend
       * data type is inferred from @c FormLanguage::Traits<GridFunctionType>::DataType
       * and must match the solver's @c SolutionType.
       *
       * This is equivalent to calling @c solve(gf.getData()).
       *
       * @tparam GridFunctionType GridFunction type whose Traits must expose DataType.
       * @param[in,out] gf  GridFunction holding the initial guess / final solution.
       */
      template <class GridFunctionType,
        class = std::enable_if_t<
          std::is_same_v<
            typename FormLanguage::Traits<std::decay_t<GridFunctionType>>::DataType,
            SolutionType>>>
      void solve(GridFunctionType& gf)
      {
        static_assert(
          std::is_same_v<
            typename FormLanguage::Traits<std::decay_t<GridFunctionType>>::DataType,
            SolutionType>,
          "GridFunction DataType must match the solver's SolutionType.");
        solve(gf.getData());
      }

    protected:
      ProblemBaseType& getProblem() noexcept
      {
        return m_solver.get().getProblem();
      }

      const ProblemBaseType& getProblem() const noexcept
      {
        return m_solver.get().getProblem();
      }

    private:
      std::reference_wrapper<LinearSolver> m_solver;
  };

  /**
   * @brief Newton solver for state-dependent tangential problems.
   *
   * This class implements a damped Newton method for nonlinear problems whose
   * current linearization is represented by a tangential
   * @c Variational::ProblemBase owned or referenced by the associated linear
   * solver.
   *
   * At Newton iteration @f$ k @f$, the solver reassembles the current
   * tangential problem and solves the linear system
   * @f[
   *   J(x^k)\,\delta x^k = -F(x^k),
   * @f]
   * then updates the nonlinear state according to
   * @f[
   *   x^{k+1} = x^k + \alpha\,\delta x^k,
   * @f]
   * where @f$ \alpha > 0 @f$ is a constant damping factor.
   *
   * @par Globalization
   * The damping factor is @e constant: this solver performs no line search
   * and no step-acceptance test, so convergence is only guaranteed close to
   * a solution (the classical local quadratic convergence of Newton's
   * method). Practical consequences:
   * - For problems solved far from equilibrium (large loads, strong
   *   nonlinearity), use @e incremental @e loading — ramp the load over
   *   several solves so each Newton run starts inside its basin of
   *   attraction (see @c examples/Solid/BlockGravity.cpp).
   * - Merit functions with barriers (e.g. terms that blow up as an element
   *   Jacobian determinant approaches zero) cannot be minimized reliably
   *   with a constant step: they require a backtracking line search with a
   *   validity test, which callers must currently implement around this
   *   solver.
   * - A damping factor @f$ \alpha < 1 @f$ trades convergence speed for
   *   robustness and can stabilize mildly nonconvex problems, but it is not
   *   a substitute for a line search.
   *
   * The nonlinear state vector passed to @c solve() is updated in place and is
   * used both as the initial guess and as the storage for the final iterate.
   *
   * @par Associated tangential problem
   * The solver obtains the tangential @c Variational::ProblemBase from the
   * supplied linear solver through @c LinearSolverBase::getProblem().
   * At each Newton iteration:
   * - the tangential problem is reassembled by calling @c assemble(),
   * - the assembled linear system is accessed through
   *   @c ProblemBase::getLinearSystem(),
   * - the linear correction is computed by the associated linear solver.
   *
   * @par Contract
   * This solver assumes that the associated problem is already assembled in
   * Newton form:
   * @f[
   *   J(x)\,\delta x = -F(x).
   * @f]
   * Consequently:
   * - the norm of the assembled right-hand side is interpreted as the
   *   nonlinear residual norm,
   * - the solution of the assembled linear system is interpreted as the Newton
   *   correction,
   * - the nonlinear state object passed to @c solve(x) must be the same state
   *   observed by the state-dependent integrators used to define the tangential
   *   problem.
   *
   * In typical Rodin usage, state-dependent integrators store a reference to
   * the nonlinear state object, so that reassembly automatically reflects the
   * updated iterate after each Newton step.
   *
   * @par Convergence criteria
   * Residual convergence is declared when either
   * @f[
   *   \|F(x^k)\| \le \mathrm{atol}
   * @f]
   * or
   * @f[
   *   \|F(x^k)\| \le \mathrm{rtol}\,\|F(x^0)\|.
   * @f]
   *
   * An optional step criterion may also be enabled. When the step tolerance is
   * positive, convergence is additionally declared if
   * @f[
   *   \alpha\,\|\delta x^k\| \le \mathrm{stol}.
   * @f]
   *
   * @par Diagnostics
   * The solver stores the result of the most recent solve in a @c Report object,
   * accessible through @c getReport(). An optional iteration monitor may also
   * be provided. When set, the monitor is called once per Newton iteration with
   * the current report.
   *
   * @par Default values
   * A newly constructed solver uses the following defaults:
   * - maximum iterations: 100
   * - absolute tolerance: 1e-12
   * - relative tolerance: 1e-8
   * - step tolerance: 0.0 (disabled)
   * - damping factor: 1.0
   * - monitor: none
   *
   * @par Usage
   * A typical nonlinear solve pairs a Newton problem in the increment
   * @f$ \delta u @f$ with a linear solver and iterates into the state
   * @f$ u @f$:
   * @code{.cpp}
   * auto ivw = Solid::InternalVirtualWork(law, u);  // residual + tangent
   *
   * Problem newton(du, v);
   * newton = ivw(du, v)
   *        - Integral(load, v)
   *        + DirichletBC(du, zero).on(fixedBoundary);
   *
   * SparseLU linearSolver(newton);
   * NewtonSolver solver(linearSolver);
   * solver.setMaxIterations(50)
   *       .setAbsoluteTolerance(1e-10)
   *       .setRelativeTolerance(1e-8);
   * solver.solve(u);   // u is the initial guess and receives the solution
   * @endcode
   *
   * @tparam LinearSolver
   *   Type of the linear solver used at each Newton step. It must have a
   *   @c FormLanguage::Traits specialization exposing @c LinearSystemType.
   */
  template <class LinearSolver>
  class NewtonSolver final
    : public NewtonSolverBase<LinearSolver>
  {
    public:
      /**
       * @brief Reason why the most recent solve terminated.
       */
      enum class ConvergedReason
      {
        /**
         * @brief The absolute residual tolerance was satisfied.
         */
        AbsoluteTolerance,

        /**
         * @brief The relative residual tolerance was satisfied.
         */
        RelativeTolerance,

        /**
         * @brief The step tolerance was satisfied.
         */
        StepTolerance,

        /**
         * @brief The maximum number of Newton iterations was reached.
         */
        MaxIterations,

        ResidualNormIsNotFinite,

        StepNormIsNotFinite,

        StepRejected,

        UserCriterion
      };

      /**
       * @brief Summary of the most recent call to @c solve().
       *
       * This structure is reset at the beginning of each solve and filled
       * progressively during the Newton iterations.
       */
      struct Report
      {
        /**
         * @brief Number of completed Newton iterations before termination.
         *
         * This is the zero-based index of the final processed Newton iteration
         * for converged solves. On termination by maximum iterations, it is set
         * to @c getMaxIterations().
         */
        size_t iterations = 0;

        /**
         * @brief Residual norm at the initial iterate.
         */
        Real initialResidual = 0.0;

        /**
         * @brief Residual norm of the last assembled tangential system.
         */
        Real finalResidual = 0.0;

        /**
         * @brief Norm of the last damped Newton correction.
         *
         * This quantity is
         * @f[
         *   \alpha\,\|\delta x\|.
         * @f]
         * It is zero until a linear correction has been computed.
         */
        Real finalStepNorm = 0.0;

        /**
         * @brief Damping factor used during the solve.
         */
        Real dampingFactor = 1.0;

        /**
         * @brief Reason for termination.
         *
         * The default value corresponds to non-convergence by exhaustion of the
         * iteration budget.
         */
        ConvergedReason reason = ConvergedReason::MaxIterations;

        /**
         * @brief Whether the solve terminated by a convergence criterion.
         */
        bool converged = false;
      };

      /**
       * @brief Type of the optional iteration monitor.
       *
       * The monitor is called with the current @c Report after the residual has
       * been evaluated and, when applicable, after the step norm has been
       * computed.
       */
      using Monitor = std::function<void(const Report&)>;

      /// @brief Parent class type.
      using Parent = NewtonSolverBase<LinearSolver>;
      /// @brief Linear system type.
      using LinearSystemType = typename Parent::LinearSystemType;
      /// @brief Associated problem base type.
      using ProblemBaseType = typename Parent::ProblemBaseType;
      /// @brief Solution vector type.
      using SolutionType = typename Parent::SolutionType;
      using LinearSolverType = LinearSolver;

      /**
       * @brief Outcome of a single Newton step.
       */
      struct StepResult
      {
        bool accepted = true;
        bool converged = false;
        Real stepNorm = 0.0;
      };

      using StepPolicy =
        std::function<StepResult(SolutionType&, LinearSystemType&, Report&)>;

      using Parent::solve;

      /**
       * @brief Constructs a Newton solver from a linear solver reference.
       *
       * The associated tangential problem is obtained from the linear solver and
       * reused throughout the Newton iterations.
       *
       * Default values:
       * - maximum iterations: 100
       * - absolute tolerance: 1e-12
       * - relative tolerance: 1e-8
       * - step tolerance: 0.0
       * - damping factor: 1.0
       * - monitor: none
       *
       * @param solver Linear solver used for each tangential solve.
       */
      explicit NewtonSolver(LinearSolver& solver)
        : Parent(solver),
          m_maxIt(100),
          m_atol(1e-12),
          m_rtol(1e-8),
          m_stol(0.0),
          m_alpha(1.0),
          m_monitor(std::nullopt),
          m_stepPolicy(std::nullopt)
      {}

      ~NewtonSolver() override = default;

      NewtonSolver* copy() const noexcept override
      {
        return new NewtonSolver(*this);
      }

      /**
       * @brief Sets the absolute residual tolerance.
       *
       * The absolute tolerance must be finite and non-negative.
       *
       * @param atol Absolute residual tolerance.
       * @returns Reference to this solver.
       */
      NewtonSolver& setAbsoluteTolerance(Real atol)
      {
        if (!std::isfinite(atol) || atol < 0.0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Absolute tolerance must be finite and non-negative."
            << Alert::Raise;
        }
        m_atol = atol;
        return *this;
      }

      /**
       * @brief Gets the absolute residual tolerance.
       *
       * Default value: 1e-12.
       *
       * @returns Current absolute tolerance.
       */
      Real getAbsoluteTolerance() const noexcept
      {
        return m_atol;
      }

      /**
       * @brief Sets the relative residual tolerance.
       *
       * The relative tolerance must be finite and non-negative.
       *
       * @param rtol Relative residual tolerance.
       * @returns Reference to this solver.
       */
      NewtonSolver& setRelativeTolerance(Real rtol)
      {
        if (!std::isfinite(rtol) || rtol < 0.0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Relative tolerance must be finite and non-negative."
            << Alert::Raise;
        }
        m_rtol = rtol;
        return *this;
      }

      /**
       * @brief Gets the relative residual tolerance.
       *
       * Default value: 1e-8.
       *
       * @returns Current relative tolerance.
       */
      Real getRelativeTolerance() const noexcept
      {
        return m_rtol;
      }

      /**
       * @brief Sets the step tolerance.
       *
       * The step tolerance must be finite and non-negative. A value of zero
       * disables the step-based convergence test.
       *
       * @param stol Step tolerance.
       * @returns Reference to this solver.
       */
      NewtonSolver& setStepTolerance(Real stol)
      {
        if (!std::isfinite(stol) || stol < 0.0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Step tolerance must be finite and non-negative."
            << Alert::Raise;
        }
        m_stol = stol;
        return *this;
      }

      /**
       * @brief Gets the step tolerance.
       *
       * Default value: 0.0, which disables the step criterion.
       *
       * @returns Current step tolerance.
       */
      Real getStepTolerance() const noexcept
      {
        return m_stol;
      }

      /**
       * @brief Sets the constant damping factor.
       *
       * The damping factor must be finite and strictly positive. A value of
       * 1.0 corresponds to the standard undamped Newton update.
       *
       * @param alpha Damping factor.
       * @returns Reference to this solver.
       */
      NewtonSolver& setDampingFactor(Real alpha)
      {
        if (!std::isfinite(alpha) || alpha <= 0.0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Damping factor must be finite and strictly positive."
            << Alert::Raise;
        }
        m_alpha = alpha;
        return *this;
      }

      /**
       * @brief Gets the damping factor.
       *
       * Default value: 1.0.
       *
       * @returns Current damping factor.
       */
      Real getDampingFactor() const noexcept
      {
        return m_alpha;
      }

      /**
       * @brief Sets the maximum number of Newton iterations.
       *
       * The iteration count must be strictly positive.
       *
       * @param maxIt Maximum number of Newton iterations.
       * @returns Reference to this solver.
       */
      NewtonSolver& setMaxIterations(size_t maxIt)
      {
        if (maxIt == 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Maximum iterations must be strictly positive."
            << Alert::Raise;
        }
        m_maxIt = maxIt;
        return *this;
      }

      /**
       * @brief Gets the maximum number of Newton iterations.
       *
       * Default value: 100.
       *
       * @returns Maximum number of iterations.
       */
      size_t getMaxIterations() const noexcept
      {
        return m_maxIt;
      }

      /**
       * @brief Sets the optional iteration monitor.
       *
       * Passing a callable installs a monitor invoked once per Newton
       * iteration. Passing @c std::nullopt disables monitoring.
       *
       * Default value: none.
       *
       * @param monitor Optional monitor callback.
       * @returns Reference to this solver.
       */
      NewtonSolver& setMonitor(Optional<Monitor> monitor)
      {
        m_monitor = std::move(monitor);
        return *this;
      }

      /**
       * @brief Gets the optional iteration monitor.
       *
       * @returns The currently installed monitor, if any.
       */
      const Optional<Monitor>& getMonitor() const noexcept
      {
        return m_monitor;
      }

      NewtonSolver& setStepPolicy(Optional<StepPolicy> policy)
      {
        m_stepPolicy = std::move(policy);
        return *this;
      }

      const Optional<StepPolicy>& getStepPolicy() const noexcept
      {
        return m_stepPolicy;
      }

      /**
       * @brief Gets the report of the most recent solve.
       *
       * @returns Diagnostic report for the last call to @c solve().
       */
      const Report& getReport() const noexcept
      {
        return m_report;
      }

      /**
       * @brief Returns whether the most recent solve converged.
       *
       * @returns @c true if the last solve terminated by a convergence
       * criterion, @c false otherwise.
       */
      bool converged() const noexcept
      {
        return m_report.converged;
      }

      /**
       * @brief Solves the nonlinear problem starting from @p x.
       *
       * The vector @p x is used both as the initial guess and as the storage for
       * the final nonlinear iterate. The associated tangential problem is
       * reassembled at each Newton iteration.
       *
       * The solve proceeds as follows:
       * - assemble the tangential problem,
       * - compute the residual norm from the assembled right-hand side,
       * - test convergence,
       * - solve the tangential linear system,
       * - compute the damped step norm,
       * - test the optional step criterion,
       * - update @p x with the damped Newton correction.
       *
       * The internal report is reset at the beginning of the solve.
       *
       * @param[in,out] x Nonlinear state vector updated in place.
       *
       * @throws Alert::MemberFunctionException
       *   If the residual norm or the step norm is not finite.
       */
      void solve(SolutionType& x) override
      {
        m_report = Report{};
        m_report.dampingFactor = m_alpha;

        Real r0 = 0.0;

        for (size_t it = 0; it < m_maxIt; ++it)
        {
          auto& pb = this->getProblem();
          pb.assemble();

          auto& linearSystem = pb.getLinearSystem();
          const Real r = linearSystem.getVector().norm();

          if (!std::isfinite(r))
          {
            m_report.iterations = it;
            m_report.initialResidual = r0;
            m_report.finalResidual = r;
            m_report.finalStepNorm = 0.0;
            m_report.reason = ConvergedReason::ResidualNormIsNotFinite;
            m_report.converged = false;
            notify();
            return;
          }

          if (it == 0)
            r0 = r;

          m_report.iterations = it;
          m_report.initialResidual = r0;
          m_report.finalResidual = r;
          m_report.finalStepNorm = 0.0;

          if (r <= m_atol)
          {
            m_report.reason = ConvergedReason::AbsoluteTolerance;
            m_report.converged = true;
            notify();
            return;
          }

          if (r0 > 0.0 && r <= m_rtol * r0)
          {
            m_report.reason = ConvergedReason::RelativeTolerance;
            m_report.converged = true;
            notify();
            return;
          }

          this->getLinearSolver().solve();

          if (m_stepPolicy)
          {
            const StepResult step =
              (*m_stepPolicy)(x, linearSystem, m_report);
            if (!std::isfinite(step.stepNorm))
            {
              m_report.reason = ConvergedReason::StepNormIsNotFinite;
              m_report.converged = false;
              notify();
              return;
            }

            m_report.finalStepNorm = step.stepNorm;
            if (!step.accepted)
            {
              m_report.reason = ConvergedReason::StepRejected;
              m_report.converged = false;
              notify();
              return;
            }

            if (step.converged)
            {
              m_report.reason = ConvergedReason::UserCriterion;
              m_report.converged = true;
              notify();
              return;
            }

            notify();
            continue;
          }

          const Real dxNorm = m_alpha * linearSystem.getSolution().norm();
          if (!std::isfinite(dxNorm))
          {
            m_report.reason = ConvergedReason::StepNormIsNotFinite;
            m_report.converged = false;
            notify();
            return;
          }

          m_report.finalStepNorm = dxNorm;

          if (m_stol > 0.0 && dxNorm <= m_stol)
          {
            m_report.reason = ConvergedReason::StepTolerance;
            m_report.converged = true;
            notify();
            x += m_alpha * linearSystem.getSolution();
            return;
          }

          notify();
          x += m_alpha * linearSystem.getSolution();
        }

        m_report.iterations = m_maxIt;
        m_report.reason = ConvergedReason::MaxIterations;
        m_report.converged = false;
      }

    private:
      /**
       * @brief Invokes the optional monitor with the current report.
       */
      void notify() const
      {
        if (m_monitor)
          (*m_monitor)(m_report);
      }

    private:
      /**
       * @brief Maximum number of Newton iterations.
       *
       * Default value: 100.
       */
      size_t m_maxIt;

      /**
       * @brief Absolute residual tolerance.
       *
       * Default value: 1e-12.
       */
      Real m_atol;

      /**
       * @brief Relative residual tolerance.
       *
       * Default value: 1e-8.
       */
      Real m_rtol;

      /**
       * @brief Step tolerance.
       *
       * Default value: 0.0. A zero value disables the step criterion.
       */
      Real m_stol;

      /**
       * @brief Constant damping factor.
       *
       * Default value: 1.0.
       */
      Real m_alpha;

      /**
       * @brief Optional iteration monitor.
       *
       * Default value: none.
       */
      Optional<Monitor> m_monitor;

      Optional<StepPolicy> m_stepPolicy;

      /**
       * @brief Diagnostic report of the most recent solve.
       */
      Report m_report;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for NewtonSolver.
   *
   * Allows writing:
   * @code{.cpp}
   * SparseLU solver(tangent);
   * NewtonSolver newton(solver);
   * @endcode
   */
  template <class LS>
  NewtonSolver(LS&) -> NewtonSolver<LS>;
}

#endif
