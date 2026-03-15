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
#include <utility>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Copyable.h"
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
   * @tparam LinearSystem Type of linear system assembled at each Newton step.
   * @tparam LinearSolver Type of the linear solver used at each Newton step.
   */
  template <class LinearSystem, class LinearSolver>
  class NewtonSolverBase : public Copyable
  {
    public:
      using LinearSystemType = LinearSystem;
      using LinearSolverType = LinearSolver;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using SolutionType = typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using Parent = Copyable;

      virtual ~NewtonSolverBase() = default;

      explicit NewtonSolverBase(ProblemBaseType& pb)
        : m_pb(pb)
      {}

      virtual const LinearSolver& getLinearSolver() const = 0;

      virtual LinearSolver& getLinearSolver() = 0;

      /**
       * @brief Solve a nonlinear system starting from the initial guess stored in @p x.
       * @param[in,out] x On entry: initial guess. On exit: final Newton iterate.
       */
      virtual void solve(SolutionType& x) = 0;

    protected:
      ProblemBaseType& getProblem() noexcept
      {
        return m_pb.get();
      }

      const ProblemBaseType& getProblem() const noexcept
      {
        return m_pb.get();
      }

    private:
      std::reference_wrapper<ProblemBaseType> m_pb;
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
   * @tparam LinearSystem Type of linear system assembled at each Newton step.
   * @tparam LinearSolver Type of the linear solver used at each Newton step.
   */
  template <class LinearSystem, class LinearSolver>
  class NewtonSolver final
    : public NewtonSolverBase<LinearSystem, LinearSolver>
  {
    public:
      using Parent = NewtonSolverBase<LinearSystem, LinearSolver>;

      using LinearSystemType = LinearSystem;
      using ProblemBaseType = typename Parent::ProblemBaseType;
      using SolutionType = typename Parent::SolutionType;
      using LinearSolverType = LinearSolver;

      explicit NewtonSolver(ProblemBaseType& pb, const LinearSolver& linearSolver)
        : Parent(pb),
          m_linearSolver(linearSolver),
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
          if (it == 0)
            r0 = r;

          if (converged(r, r0))
          {
            return;
          }

          m_linearSolver.solve();
          x += m_alpha * linearSystem.getSolution();
        }
      }

    private:
      bool converged(Real r, Real r0) const
      {
        return r <= m_atol || (r0 > 0.0 && r <= m_rtol * r0);
      }

    private:
      LinearSolver m_linearSolver;
      size_t m_maxIt;
      Real m_atol;
      Real m_rtol;
      Real m_alpha;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD guide for NewtonSolver from a Rodin problem and linear solver.
   *
   * Deduces:
   * - LinearSystemType = LinearSystem
   * - LinearSolverType = LinearSolver
   */
  template <class LinearSystem, class LinearSolver>
  NewtonSolver(Variational::ProblemBase<LinearSystem>&, const LinearSolver&)
    -> NewtonSolver<LinearSystem, LinearSolver>;
}

#endif
