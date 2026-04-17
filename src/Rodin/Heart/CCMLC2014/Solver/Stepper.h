/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_SOLVER_STEPPER_H
#define RODIN_HEART_CCMLC2014_SOLVER_STEPPER_H

#include <algorithm>
#include <cassert>
#include <cmath>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/ForwardDecls.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/PartialPivLU.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/Problem.h"
#include "Rodin/Heart/CCMLC2014/Model/History.h"
#include "Rodin/Heart/CCMLC2014/Model/State.h"
#include "Rodin/Heart/CCMLC2014/Numerics/Jacobian.h"
#include "Rodin/Heart/CCMLC2014/Numerics/Residual.h"

namespace Rodin::Heart::CCMLC2014::Solver
{
  template <class PassiveEnergyLaw, class PassiveLaw>
  class StepperT
  {
    public:
      using Scalar = Real;
      using DenseMatrix = Math::Matrix<Scalar>;
      using DenseVector = Math::Vector<Scalar>;
      using DenseLinearSystem = Math::LinearSystem<DenseMatrix, DenseVector>;

      using State = Model::StateT<Scalar>;
      using Input = Model::InputT<Scalar, PassiveEnergyLaw>;
      using Report = Model::ReportT<Scalar, DenseLinearSystem>;
      using EvalData = Model::EvalDataT<Scalar>;
      using History = Model::HistoryT<State>;

      class Problem final : public Variational::ProblemBase<DenseLinearSystem>
      {
        public:
          using Parent = Variational::ProblemBase<DenseLinearSystem>;
          using ProblemBodyType = typename Parent::ProblemBodyType;

          explicit Problem(const Input& in)
            : m_input(in)
          {
            m_system.getOperator().resize(Model::NVAR, Model::NVAR);
            m_system.getVector().resize(Model::NVAR);
            m_system.getSolution().resize(Model::NVAR);
          }

          Parent& operator=(const ProblemBodyType&) override
          {
            return *this;
          }

          Problem& assemble() override
          {
            assert(m_xCurrent);

            auto& A = m_system.getOperator();
            auto& b = m_system.getVector();
            auto& s = m_system.getSolution();
            s.setZero();

            EvalData d;
            Numerics::buildEvalData<PassiveLaw>(m_input, *m_xCurrent, m_history.n, m_history.nm1, m_time, m_dt, d);
            Numerics::evaluateDynamicResidual(m_input, d, b);
            b = -b;
            Numerics::evaluateDynamicJacobian(m_input, A, d, m_dt);
            return *this;
          }

          void solve(::Rodin::Solver::LinearSolverBase<DenseLinearSystem>& solver) override
          {
            solver.solve(m_system);
          }

          DenseLinearSystem& getLinearSystem() override { return m_system; }
          const DenseLinearSystem& getLinearSystem() const override { return m_system; }
          Problem* copy() const noexcept override { return new Problem(*this); }

          void setCurrent(DenseVector& xCurrent) { m_xCurrent = &xCurrent; }

          void setStepData(const History& history, Scalar time, Scalar dt)
          {
            m_history = history;
            m_time = time;
            m_dt = dt;
          }

        private:
          Input m_input;
          DenseVector* m_xCurrent = nullptr;
          DenseLinearSystem m_system;

          History m_history;
          Scalar m_time = 0.0;
          Scalar m_dt = 0.0;
      };

      explicit StepperT(const Input& input)
        : m_input(input), m_problem(input), m_solver(m_problem), m_newton(m_solver)
      {
        m_x = packUnknowns(m_state);
      }

      StepperT& setAbsoluteTolerance(const Scalar atol)
      {
        m_atol = atol;
        m_newton.setAbsoluteTolerance(atol);
        return *this;
      }

      StepperT& setRelativeTolerance(const Scalar rtol)
      {
        m_rtol = rtol;
        m_newton.setRelativeTolerance(rtol);
        return *this;
      }

      StepperT& setStepTolerance(const Scalar stol)
      {
        m_stol = stol;
        m_newton.setStepTolerance(stol);
        return *this;
      }

      StepperT& setMaxIterations(const size_t maxIt)
      {
        m_maxIt = maxIt;
        m_newton.setMaxIterations(maxIt);
        return *this;
      }

      StepperT& setDampingFactor(const Scalar alpha)
      {
        m_damping = alpha;
        m_newton.setDampingFactor(alpha);
        return *this;
      }

      void initialize(const State& initial)
      {
        m_state = initial;
        m_x = packUnknowns(m_state);
        m_report = {};

        if (initial.gamma > Scalar(0))
        {
          m_state.gamma = initial.gamma;
          m_state.beta = initial.beta;
          m_state.kc = initial.gamma * initial.gamma;
          m_state.tauc = initial.gamma * initial.beta;
        }
        else
        {
          if (initial.kc > Scalar(0))
          {
            m_state.gamma = std::sqrt(initial.kc);
            m_state.beta = (m_state.gamma > Scalar(0))
                         ? (initial.tauc / m_state.gamma)
                         : Scalar(0);
          }
          else
          {
            m_state.gamma = std::sqrt(std::max<Scalar>(m_input.initActiveStiffness, Scalar(0)));
            m_state.beta =
              (m_state.gamma > Scalar(0))
              ? (m_input.initActiveStress / m_state.gamma)
              : Scalar(0);
            m_state.kc = m_state.gamma * m_state.gamma;
            m_state.tauc = m_state.gamma * m_state.beta;
          }
        }

        if (std::abs(initial.ec) > Scalar(0))
          m_state.ec = initial.ec;
        else
          m_state.ec = m_input.initFibDef;

        m_history.nm1 = m_state;
        m_x = packUnknowns(m_state);
      }

      Report step(const Scalar dt)
      {
        assert(dt > Scalar(0));

        m_history.n = m_state;

        m_x = packUnknowns(m_history.n);

        m_problem.setCurrent(m_x);
        m_problem.setStepData(m_history, m_history.n.t + dt, dt);

        m_newton.solve(m_x);

        const auto& nr = m_newton.getReport();
        m_report.converged = nr.converged;
        m_report.iterations = nr.iterations;
        m_report.finalResidual = nr.final_residual;
        m_report.finalStepNorm = nr.final_step_norm;
        m_report.reason = nr.reason;

        if (nr.converged)
        {
          EvalData d;
          Numerics::buildEvalData<PassiveLaw>(m_input, m_x, m_history.n, m_history.nm1, m_history.n.t + dt, dt, d);

          m_history.nm2 = m_history.nm1;
          m_history.nm1 = m_history.n;
          m_state = unpackUnknownsIntoState(m_x, m_state, m_history.n.t + dt);
          m_state.v = (m_state.y - m_history.n.y) / dt;
          m_state.ec = d.active.fib1;
          m_state.gamma = d.active.gammaNew;
          m_state.beta = d.active.betaNew;
          m_state.kc = m_state.gamma * m_state.gamma;
          m_state.tauc = m_state.gamma * m_state.beta;
        }

        return m_report;
      }

      const State& getState() const noexcept { return m_state; }
      const History& getHistory() const noexcept { return m_history; }
      const Report& getReport() const noexcept { return m_report; }
      const DenseVector& getUnknowns() const noexcept { return m_x; }

    private:
      static DenseVector packUnknowns(const State& s)
      {
        DenseVector x(Model::NVAR);
        x[Model::DISP] = s.y;
        x[Model::PV]   = s.pv;
        x[Model::PAR]  = s.par;
        x[Model::PD]   = s.pd;
        return x;
      }

      static State unpackUnknownsIntoState(const DenseVector& x, const State& base, Scalar t)
      {
        State s = base;
        s.y = x[Model::DISP];
        s.pv = x[Model::PV];
        s.par = x[Model::PAR];
        s.pd = x[Model::PD];
        s.t = t;
        return s;
      }

    private:
      Input m_input;
      State m_state;
      History m_history;

      Scalar m_atol = Scalar(1e-10);
      Scalar m_rtol = Scalar(1e-10);
      Scalar m_stol = Scalar(1e-12);
      size_t m_maxIt = 50;
      Scalar m_damping = Scalar(1.0);

      Report m_report;
      DenseVector m_x;

      Problem m_problem;
      ::Rodin::Solver::PartialPivLU<DenseLinearSystem> m_solver;
      ::Rodin::Solver::NewtonSolver<::Rodin::Solver::PartialPivLU<DenseLinearSystem>> m_newton;
  };
}

#endif
