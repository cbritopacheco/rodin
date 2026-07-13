/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Stepper.h
 * @brief Nonlinear time-stepper for the CCMLC2014 reduced 0D ventricular model.
 *
 * Wraps the DynamicSystem assembler inside a Variational::ProblemBase
 * and drives the Newton solver, advancing the coupled system state
 * from one time step to the next (Caruel et al. 2014, §4).
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
#include "Rodin/Heart/CCMLC2014/Numerics/DynamicSystem.h"

namespace Rodin::Heart::CCMLC2014::Solver
{
  /**
   * @brief Time-stepper for the coupled 0D CCMLC2014 ventricular model.
   *
   * Manages the Newton iteration for one time step, delegating residual and
   * Jacobian assembly to the Numerics::DynamicSystem class, and updating
   * the model state on convergence.
   *
   * @tparam PassiveEnergyLaw Passive-energy law used by the constitutive response.
   * @tparam PassiveLaw Passive stress operator using passive-energy derivatives.
   */
  template <class PassiveEnergyLaw, class PassiveLaw>
  class StepperT
  {
    public:
      /// @brief Scalar type used by the reduced model.
      using Scalar = Real;
      /// @brief Dense matrix type used for Newton systems.
      using DenseMatrix = Math::Matrix<Scalar>;
      /// @brief Dense vector type used for states and residuals.
      using DenseVector = Math::Vector<Scalar>;
      /// @brief Dense linear system assembled by the reduced model.
      using DenseLinearSystem = Math::LinearSystem<DenseMatrix, DenseVector>;

      /// @brief Dynamic state type.
      using State = Model::StateT<Scalar>;
      /// @brief Model input and parameter type.
      using Input = Model::InputT<Scalar, PassiveEnergyLaw>;
      /// @brief Nonlinear solve report type.
      using Report = Model::ReportT<Scalar, DenseLinearSystem>;
      /// @brief Residual/Jacobian evaluation cache type.
      using EvalData = Model::EvalDataT<Scalar>;
      /// @brief Time-history type.
      using History = Model::HistoryT<State>;

      /**
       * @brief Variational problem adapter for Newton iteration.
       *
       * Bridges the Numerics::DynamicSystem assembler into the
       * Rodin::Variational::ProblemBase interface expected by the
       * Newton solver.
       */
      class Problem final : public Variational::ProblemBase<DenseLinearSystem>
      {
        public:
          /// @brief Parent class type.
          using Parent = Variational::ProblemBase<DenseLinearSystem>;
          /// @brief Problem body type accepted by the variational interface.
          using ProblemBodyType = typename Parent::ProblemBodyType;

          /**
           * @brief Constructs the problem adapter from model inputs.
           * @param in Model input parameters.
           */
          explicit Problem(const Input& in)
            : m_dynamicSystem(in)
          {
            m_system.getOperator().resize(
                Model::NumberOfVariables, Model::NumberOfVariables);
            m_system.getVector().resize(Model::NumberOfVariables);
            m_system.getSolution().resize(Model::NumberOfVariables);
          }

          Parent& operator=(const ProblemBodyType&) override
          {
            return *this;
          }

          Problem& assemble() override
          {
            assert(m_xCurrent);

            auto& operatorMatrix = m_system.getOperator();
            auto& rightHandSide = m_system.getVector();
            auto& solutionIncrement = m_system.getSolution();
            solutionIncrement.setZero();

            EvalData evalData;
            m_dynamicSystem.buildEvalData(
                *m_xCurrent, m_history.n, m_history.nm1, m_time, m_dt,
                evalData);
            m_dynamicSystem.evaluateResidual(evalData, rightHandSide);
            rightHandSide = -rightHandSide;
            m_dynamicSystem.evaluateJacobian(evalData, operatorMatrix, m_dt);
            return *this;
          }

          /**
           * @brief Solves the assembled dense linear system.
           * @param solver Linear solver used for the Newton correction.
           */
          void solve(
              ::Rodin::Solver::LinearSolverBase<DenseLinearSystem>& solver) override
          {
            solver.solve(m_system);
          }

          DenseLinearSystem& getLinearSystem() override { return m_system; }
          const DenseLinearSystem& getLinearSystem() const override
          {
            return m_system;
          }
          Problem* copy() const noexcept override { return new Problem(*this); }

          /**
           * @brief Sets the current nonlinear unknown vector.
           * @param xCurrent Current unknown vector.
           */
          void setCurrent(DenseVector& xCurrent) { m_xCurrent = &xCurrent; }

          /**
           * @brief Sets state history and time-step data for assembly.
           * @param history Previous state history.
           * @param time Target time for the step.
           * @param dt Time-step size.
           */
          void setStepData(const History& history, Scalar time, Scalar dt)
          {
            m_history = history;
            m_time = time;
            m_dt = dt;
          }

        private:
          Numerics::DynamicSystem<PassiveLaw, Input> m_dynamicSystem;
          DenseVector* m_xCurrent = nullptr;
          DenseLinearSystem m_system;

          History m_history;
          Scalar m_time = 0.0;
          Scalar m_dt = 0.0;
      };

      /**
       * @brief Constructs a stepper from model input parameters.
       * @param input Model input parameters.
       */
      explicit StepperT(const Input& input)
        : m_input(input), m_problem(input), m_solver(m_problem), m_newton(m_solver)
      {
        m_x = packUnknowns(m_state);
      }

      /**
       * @brief Sets the Newton absolute residual tolerance.
       * @param atol Absolute tolerance.
       * @return Reference to this stepper.
       */
      StepperT& setAbsoluteTolerance(const Scalar atol)
      {
        m_atol = atol;
        m_newton.setAbsoluteTolerance(atol);
        return *this;
      }

      /**
       * @brief Sets the Newton relative residual tolerance.
       * @param rtol Relative tolerance.
       * @return Reference to this stepper.
       */
      StepperT& setRelativeTolerance(const Scalar rtol)
      {
        m_rtol = rtol;
        m_newton.setRelativeTolerance(rtol);
        return *this;
      }

      /**
       * @brief Sets the Newton step tolerance.
       * @param stol Step tolerance.
       * @return Reference to this stepper.
       */
      StepperT& setStepTolerance(const Scalar stol)
      {
        m_stol = stol;
        m_newton.setStepTolerance(stol);
        return *this;
      }

      /**
       * @brief Sets the maximum number of Newton iterations.
       * @param maxIt Maximum iteration count.
       * @return Reference to this stepper.
       */
      StepperT& setMaxIterations(const size_t maxIt)
      {
        m_maxIt = maxIt;
        m_newton.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Sets the Newton damping factor.
       * @param alpha Damping factor.
       * @return Reference to this stepper.
       */
      StepperT& setDampingFactor(const Scalar alpha)
      {
        m_damping = alpha;
        m_newton.setDampingFactor(alpha);
        return *this;
      }

      /**
       * @brief Initializes the dynamic state and history.
       * @param initial Initial state.
       */
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

        if (initial.w > Scalar(0))
          m_state.w = initial.w;
        else
          m_state.w = std::max<Scalar>(m_input.m0(m_state.ec), Scalar(0));

        m_history.nm1 = m_state;
        m_x = packUnknowns(m_state);
      }

      /**
       * @brief Advances the model by one time step.
       * @param dt Time-step size.
       * @return Report for the nonlinear solve.
       */
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
        m_report.finalResidual = nr.finalResidual;
        m_report.finalStepNorm = nr.finalStepNorm;
        m_report.reason = nr.reason;

        if (nr.converged)
        {
          m_history.nm2 = m_history.nm1;
          m_history.nm1 = m_history.n;
          m_state = unpackUnknownsIntoState(m_x, m_state, m_history.n.t + dt);
        }

        return m_report;
      }

      /// @brief Returns the current state.
      const State& getState() const noexcept { return m_state; }
      /// @brief Returns the stored time history.
      const History& getHistory() const noexcept { return m_history; }
      /// @brief Returns the most recent nonlinear solve report.
      const Report& getReport() const noexcept { return m_report; }
      /// @brief Returns the packed unknown vector.
      const DenseVector& getUnknowns() const noexcept { return m_x; }

      /**
       * @brief Restores a previously saved stepper state.
       * @param state Current dynamic state.
       * @param history Stored time history.
       * @param unknowns Packed unknown vector.
       * @param report Last nonlinear solve report.
       */
      void restore(
          const State& state,
          const History& history,
          const DenseVector& unknowns,
          const Report& report)
      {
        m_state = state;
        m_history = history;
        m_x = unknowns;
        m_report = report;
      }

    private:
      static DenseVector packUnknowns(const State& s)
      {
        DenseVector x(Model::NumberOfVariables);
        x[Model::RadialDisplacement] = s.y;
        x[Model::RadialVelocity] = s.v;
        x[Model::VentricularPressure] = s.pv;
        x[Model::ArterialPressure] = s.par;
        x[Model::DistalPressure] = s.pd;
        x[Model::FiberDeformation] = s.ec;
        x[Model::ActiveStiffness] = s.kc;
        x[Model::ActiveStress] = s.tauc;
        x[Model::LoadDependentRelaxation] = s.w;
        return x;
      }

      static State unpackUnknownsIntoState(
          const DenseVector& x, const State& base, Scalar t)
      {
        State s = base;
        s.y = x[Model::RadialDisplacement];
        s.v = x[Model::RadialVelocity];
        s.pv = x[Model::VentricularPressure];
        s.par = x[Model::ArterialPressure];
        s.pd = x[Model::DistalPressure];
        s.ec = x[Model::FiberDeformation];
        s.kc = normalizeActiveStiffness(x[Model::ActiveStiffness]);
        s.tauc = x[Model::ActiveStress];
        s.w = x[Model::LoadDependentRelaxation];
        s.gamma = std::sqrt(std::max<Scalar>(s.kc, Scalar(0)));
        s.beta = (s.gamma > Scalar(0)) ? (s.tauc / s.gamma) : Scalar(0);
        s.t = t;
        return s;
      }

      static Scalar normalizeActiveStiffness(Scalar kc)
      {
        if (kc < Scalar(0) && std::abs(kc) < Scalar(1e-14))
          return Scalar(0);
        return kc;
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
