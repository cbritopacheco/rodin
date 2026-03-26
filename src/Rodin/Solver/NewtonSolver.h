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
      using LinearSystemType = typename FormLanguage::Traits<LinearSolver>::LinearSystemType;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using SolutionType = typename FormLanguage::Traits<LinearSystemType>::VectorType;

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
   * @brief Newton solver for nonlinear systems.
   *
    * The tangent linear system is assembled from the associated
    * Variational::ProblemBase as:
    * @f[
    *   J(x) dx = -F(x),
    * @f]
    * and the state is updated with:
   * @f[
   *   x^{k + 1} = x^k + dx,
   * @f]
   * with @f$ x^0 @f$ given by the initial guess.
   *
   * Convergence is checked using the residual norm @f$ \|F(x^k)\| @f$ with
   * absolute and relative tolerances:
   * @f[
   *   \|F(x^k)\| \le \text{atol}
   *   \quad\text{or}\quad
   *   \|F(x^k)\| \le \text{rtol} \, \|F(x^0)\|.
   * @f]
   *
   * @tparam LinearSolver Type of the linear solver used at each Newton step.
   *   Must have a FormLanguage::Traits specialization providing LinearSystemType.
   */
  template <class LinearSolver>
  class NewtonSolver final
    : public NewtonSolverBase<LinearSolver>
  {
    public:
      using Parent = NewtonSolverBase<LinearSolver>;

      using LinearSystemType = typename Parent::LinearSystemType;
      using ProblemBaseType = typename Parent::ProblemBaseType;
      using SolutionType = typename Parent::SolutionType;
      using LinearSolverType = LinearSolver;

      using Parent::solve;

      /**
       * @brief Constructs a NewtonSolver from a linear solver reference.
       * @param solver The linear solver used at each Newton step. The
       *   associated ProblemBase is obtained from the solver via
       *   LinearSolverBase::getProblem().
       */
      explicit NewtonSolver(LinearSolver& solver)
        : Parent(solver),
          m_maxIt(100),
          m_atol(1e-12),
          m_rtol(1e-8),
          m_alpha(1.0)
      {}

      ~NewtonSolver() override = default;

      NewtonSolver* copy() const noexcept override
      {
        return new NewtonSolver(*this);
      }

      NewtonSolver& setAbsoluteTolerance(Real atol)
      {
        m_atol = atol;
        return *this;
      }

      NewtonSolver& setRelativeTolerance(Real rtol)
      {
        m_rtol = rtol;
        return *this;
      }

      NewtonSolver& setDampingFactor(Real alpha)
      {
        m_alpha = alpha;
        return *this;
      }

      Real getDampingFactor() const
      {
        return m_alpha;
      }

      Real getAbsoluteTolerance() const
      {
        return m_atol;
      }

      Real getRelativeTolerance() const
      {
        return m_rtol;
      }

      NewtonSolver& setMaxIterations(size_t maxIt)
      {
        m_maxIt = maxIt;
        return *this;
      }

      size_t getMaxIterations() const
      {
        return m_maxIt;
      }

      void solve(SolutionType& x) override
      {
        if (m_atol < 0.0)
        {
            Alert::MemberFunctionException(*this, __func__)
              << "Absolute tolerance must be non-negative."
              << Alert::Raise;
        }

        if (m_rtol < 0.0)
        {
            Alert::MemberFunctionException(*this, __func__)
              << "Relative tolerance must be non-negative."
              << Alert::Raise;
        }

        Real r0 = 0.0;
        for (size_t it = 0; it < m_maxIt; ++it)
        {
          auto& pb = this->getProblem();
          pb.assemble();

          auto& linearSystem = pb.getLinearSystem();
          const Real r = linearSystem.getVector().norm();

          if (!std::isfinite(r))
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Newton residual is not finite."
              << Alert::Raise;
          }

          if (it == 0)
            r0 = r;

          if (converged(r, r0))
          {
            return;
          }

          this->getLinearSolver().solve();
          x += m_alpha * linearSystem.getSolution();
        }

        std::cout << "Reached max iterations\n";
      }

    private:
      bool converged(Real r, Real r0) const
      {
        return r <= m_atol || (r0 > 0.0 && r <= m_rtol * r0);
      }

    private:
      size_t m_maxIt;
      Real m_atol;
      Real m_rtol;
      Real m_alpha;
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
