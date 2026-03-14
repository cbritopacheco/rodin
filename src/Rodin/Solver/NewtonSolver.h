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

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Copyable.h"
#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Types.h"

namespace Rodin::Solver
{
  /**
   * @brief Base interface for Newton-type nonlinear solvers.
   *
   * Newton solver class for solving nonlinear systems of the form
   * @f[
   * F(x) = 0,
   * @f]
   * where @f$ F @f$ is a possibly non-linear, differentiable function.
   *
   * @tparam Solution Type of the nonlinear state.
   * @tparam Function Type of the residual vector.
   * @tparam Jacobian Type of the Jacobian operator.
   * @tparam LinearSolver Type of the linear solver used at each Newton step.
   */
  template <class Solution, class Function, class Jacobian, class LinearSolver>
  class NewtonSolverBase : public Copyable
  {
    public:
      using SolutionType = Solution;
      using FunctionType = Function;
      using JacobianType = Jacobian;
      using LinearSolverType = LinearSolver;

      using ResidualAssembly = std::function<void(Function&, const Solution&)>;
      using JacobianAssembly = std::function<void(Jacobian&, const Solution&)>;

      using Parent = Copyable;

      virtual ~NewtonSolverBase() = default;

      virtual const ResidualAssembly& getFunction() const = 0;

      virtual const JacobianAssembly& getJacobian() const = 0;

      virtual const LinearSolver& getLinearSolver() const = 0;

      virtual LinearSolver& getLinearSolver() = 0;

      /**
       * @brief Solve a nonlinear system starting from the initial guess stored in @p x.
       * @param[in,out] x On entry: initial guess. On exit: final Newton iterate.
       */
      virtual void solve(Solution& x) = 0;
  };

  /**
   * @brief Newton solver for nonlinear systems.
   *
   * The nonlinear residual is assembled as F(x), the Newton system is solved
   * as:
   * @f[
   *   J(x) dx = F(x),
   * @f]
   * and the state is updated with:
   * @f[
   *   x^{k + 1} = x^k - dx,
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
   * @tparam Solution Type of the nonlinear state.
   * @tparam Function Type of the residual vector.
   * @tparam Jacobian Type of the Jacobian operator.
   * @tparam LinearSolver Type of the linear solver used at each Newton step.
   */
  template <class Solution, class Function, class Jacobian, class LinearSolver>
  class NewtonSolver final
    : public NewtonSolverBase<Solution, Function, Jacobian, LinearSolver>
  {
    public:
      using Parent = NewtonSolverBase<Solution, Function, Jacobian, LinearSolver>;

      using SolutionType = Solution;
      using FunctionType = Function;
      using JacobianType = Jacobian;
      using LinearSolverType = LinearSolver;

      using ResidualAssembly = std::function<void(FunctionType&, const SolutionType&)>;
      using JacobianAssembly = std::function<void(JacobianType&, const SolutionType&)>;
      using ResidualNorm = std::function<double(const FunctionType&)>;

      /**
       * Assumes the linear system stores:
       *   - operator: Jacobian
       *   - rhs/vector: residual F(x)
       *   - solution: Newton correction dx
       *
       * If Rodin's LinearSystem uses a different convention, this alias should be
       * adjusted accordingly.
       */
      using LinearSystemType = Math::LinearSystem<JacobianType, FunctionType>;

      explicit NewtonSolver(const LinearSolver& linearSolver)
        : m_linearSolver(linearSolver),
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

      NewtonSolver& setFunction(const ResidualAssembly& f)
      {
        m_function = f;
        return *this;
      }

      NewtonSolver& setJacobian(const JacobianAssembly& J)
      {
        m_jacobian = J;
        return *this;
      }

      NewtonSolver& setLinearSolver(const LinearSolver& linearSolver)
      {
        m_linearSolver = linearSolver;
        return *this;
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

      const ResidualAssembly& getFunction() const override
      {
        return m_function;
      }

      const JacobianAssembly& getJacobian() const override
      {
        return m_jacobian;
      }

      const LinearSolver& getLinearSolver() const override
      {
        return m_linearSolver;
      }

      LinearSolver& getLinearSolver() override
      {
        return m_linearSolver;
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

      void solve(Solution& x) override
      {
        if (!m_function)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Residual assembly not set."
            << Alert::Raise;
        }

        if (!m_jacobian)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Jacobian assembly not set."
            << Alert::Raise;
        }

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

        SolutionType xCurr = x;
        LinearSystemType linearSystem;

        FunctionType F;
        m_function(F, xCurr);

        const Real r0 = F.norm();
        if (converged(r0, r0))
        {
          x = xCurr;
          return;
        }

        for (size_t it = 0; it < m_maxIt; ++it)
        {
          auto& J = linearSystem.getOperator();
          m_jacobian(J, xCurr);

          auto& rhs = linearSystem.getVector();
          rhs = F;

          m_linearSolver.solve(linearSystem);

          xCurr = xCurr - m_alpha * linearSystem.getSolution();

          m_function(F, xCurr);
          const Real r = F.norm();

          if (converged(r, r0))
          {
            x = xCurr;
            return;
          }
        }

        x = xCurr;
      }

    private:
      bool converged(Real r, Real r0) const
      {
        return r <= m_atol || (r0 > 0.0 && r <= m_rtol * r0);
      }

    private:
      ResidualAssembly m_function;
      JacobianAssembly m_jacobian;
      LinearSolver m_linearSolver;
      size_t m_maxIt;
      Real m_atol;
      Real m_rtol;
      Real m_alpha;
  };
}

#endif
