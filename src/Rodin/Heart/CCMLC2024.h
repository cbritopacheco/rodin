/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2024_H
#define RODIN_HEART_CCMLC2024_H

#include <cassert>
#include <cmath>
#include <functional>
#include <numbers>

#include "Rodin/Types.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/Problem.h"
#include "Rodin/Solver/FiniteDifferenceProbe.h"
#include "Rodin/Solver/LDLT.h"
#include "Rodin/Solver/NewtonSolver.h"

namespace Rodin::Heart
{
  /**
   * @brief Coupled 0D cardio-circulatory model with implicit step + Newton solve.
   */
  class CCMLC2024
  {
    public:
      using Scalar = Real;
      using DenseMatrix = Math::Matrix<Scalar>;
      using DenseVector = Math::Vector<Scalar>;
      using DenseLinearSystem = Math::LinearSystem<DenseMatrix, DenseVector>;

      enum Variable : size_t
      {
        Y = 0,
        V,
        PV,
        PAR,
        PD,
        EC,
        KC,
        TAUC,
        NVAR
      };

      struct State
      {
        Scalar y = 0.0;
        Scalar v = 0.0;
        Scalar pv = 0.0;
        Scalar par = 0.0;
        Scalar pd = 0.0;
        Scalar ec = 0.0;
        Scalar kc = 0.0;
        Scalar tauc = 0.0;
        Scalar t = 0.0;
      };

      struct Input
      {
        Scalar rho = 1.0;
        Scalar d0 = 1.0;
        Scalar R0 = 1.0;
        Scalar eta = 0.0;

        Scalar Es = 1.0;
        Scalar mu = 1.0;
        Scalar alpha = 0.0;
        Scalar n0 = 0.0;
        Scalar k0 = 0.0;
        Scalar sigma0 = 0.0;
        Scalar w = 1.0;

        Scalar Cp = 1.0;
        Scalar Cd = 1.0;
        Scalar Rp = 1.0;
        Scalar Rd = 1.0;
        Scalar Psv = 0.0;

        std::function<Scalar(Scalar)> e1D =
          [](Scalar C) { return Scalar(0.5) * (C * C - Scalar(1)); };

        std::function<Scalar(Scalar, Scalar, Scalar)> dWe_dJ1 =
          [](Scalar, Scalar, Scalar) { return Scalar(0); };
        std::function<Scalar(Scalar, Scalar, Scalar)> dWe_dJ2 =
          [](Scalar, Scalar, Scalar) { return Scalar(0); };
        std::function<Scalar(Scalar, Scalar, Scalar)> dWe_dJ4 =
          [](Scalar, Scalar, Scalar) { return Scalar(0); };

        std::function<Scalar(Scalar)> ubar =
          [](Scalar) { return Scalar(0); };

        std::function<Scalar(Scalar)> pAt =
          [](Scalar) { return Scalar(0); };

        std::function<Scalar(Scalar, Scalar, Scalar)> valveFlow =
          [](Scalar pv, Scalar par, Scalar) { return pv > par ? pv - par : Scalar(0); };
      };

      struct Report
      {
        bool converged = false;
        size_t iterations = 0;
        Scalar finalResidual = 0.0;
        Scalar finalStepNorm = 0.0;
        Solver::NewtonSolver<Solver::LDLT<DenseLinearSystem>>::ConvergenceReason reason =
          Solver::NewtonSolver<Solver::LDLT<DenseLinearSystem>>::ConvergenceReason::MaxIterations;
      };

    private:
      class Problem final : public Variational::ProblemBase<DenseLinearSystem>
      {
        public:
          using Parent = Variational::ProblemBase<DenseLinearSystem>;
          using ProblemBodyType = typename Parent::ProblemBodyType;

          explicit Problem(const Input& in)
            : m_input(in), m_dt(0.0), m_time(0.0)
          {
            m_system.getOperator().resize(NVAR, NVAR);
            m_system.getVector().resize(NVAR);
            m_system.getSolution().resize(NVAR);
            m_xPrev.resize(NVAR);
            m_xCurrent = nullptr;
            m_xPrev.setZero();
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

            DenseVector R(NVAR);
            computeResidual(*m_xCurrent, R);
            b = -R;

            A = Solver::FiniteDifferenceProbe::jacobian(
              *m_xCurrent,
              [&](const DenseVector& x, DenseVector& out) { computeResidual(x, out); });
            return *this;
          }

          void solve(Solver::LinearSolverBase<DenseLinearSystem>& solver) override
          {
            solver.solve(m_system);
          }

          DenseLinearSystem& getLinearSystem() override
          {
            return m_system;
          }

          const DenseLinearSystem& getLinearSystem() const override
          {
            return m_system;
          }

          Problem* copy() const noexcept override
          {
            return new Problem(*this);
          }

          void setPrevious(const DenseVector& xPrev)
          {
            m_xPrev = xPrev;
          }

          void setCurrent(DenseVector& xCurrent)
          {
            m_xCurrent = &xCurrent;
          }

          void setTime(const Scalar t)
          {
            m_time = t;
          }

          void setTimeStep(const Scalar dt)
          {
            m_dt = dt;
          }

        private:
          static Scalar posPart(const Scalar x)
          {
            return (x > 0) ? x : 0;
          }

          static Scalar negPart(const Scalar x)
          {
            return (x < 0) ? -x : 0;
          }

          void computeResidual(const DenseVector& x, DenseVector& R) const
          {
            const Scalar y = x[Y];
            const Scalar v = x[V];
            const Scalar pv = x[PV];
            const Scalar par = x[PAR];
            const Scalar pd = x[PD];
            const Scalar ec = x[EC];
            const Scalar kc = x[KC];
            const Scalar tauc = x[TAUC];

            const Scalar yPrev = m_xPrev[Y];
            const Scalar vPrev = m_xPrev[V];
            const Scalar parPrev = m_xPrev[PAR];
            const Scalar pdPrev = m_xPrev[PD];
            const Scalar ecPrev = m_xPrev[EC];
            const Scalar kcPrev = m_xPrev[KC];
            const Scalar taucPrev = m_xPrev[TAUC];

            const Scalar C = 1.0 + y / m_input.R0;
            const Scalar invC = 1.0 / C;
            const Scalar C2 = C * C;
            const Scalar C6inv = invC * invC * invC * invC * invC * invC;
            const Scalar cdot = v / m_input.R0;

            const Scalar j1 = C2 + 2.0 * invC;
            const Scalar j2 = 2.0 * C + invC * invC;
            const Scalar j4 = C2;

            const Scalar dwe1 = m_input.dWe_dJ1(j1, j2, j4);
            const Scalar dwe2 = m_input.dWe_dJ2(j1, j2, j4);
            const Scalar dwe4 = m_input.dWe_dJ4(j1, j2, j4);

            const Scalar e1d = m_input.e1D(C);
            const Scalar ecDot = (ec - ecPrev) / m_dt;

            const Scalar denom1 = 1.0 + 2.0 * ec;
            const Scalar sigma1D = m_input.Es * (e1d - ec) / (denom1 * denom1);

            const Scalar sigmaSph =
              sigma1D
              + 4.0 * (1.0 - invC * invC * invC) * (dwe1 + C * dwe2)
              + 2.0 * dwe4
              + 2.0 * m_input.eta * cdot * (1.0 - 2.0 * C6inv);

            const Scalar pat = m_input.pAt(m_time + m_dt);
            const Scalar Q = m_input.valveFlow(pv, par, pat);

            const Scalar ubar = m_input.ubar(m_time + m_dt);
            const Scalar A =
              posPart(ubar)
              + m_input.w * negPart(ubar)
              + m_input.alpha * std::abs(ecDot);

            const Scalar volumeRate =
              4.0 * std::numbers::pi_v<Scalar> * m_input.R0 * m_input.R0 * C2 * v;

            R[Y] = y - yPrev - m_dt * v;

            R[V] =
              m_input.rho * m_input.d0 * (v - vPrev) / m_dt
              + (m_input.d0 / m_input.R0) * C * sigmaSph
              - pv * C2;

            R[PV] =
              tauc + m_input.mu * ecDot
              - m_input.Es * ((e1d - ec) * (1.0 + 2.0 * e1d))
                / (denom1 * denom1 * denom1);

            R[PAR] = (kc - kcPrev) / m_dt + A * kc - m_input.n0 * m_input.k0 * posPart(ubar);

            R[PD] =
              (tauc - taucPrev) / m_dt
              + A * tauc
              - m_input.n0 * m_input.sigma0 * posPart(ubar)
              - kc * ecDot;

            R[EC] = volumeRate + Q;

            R[KC] =
              m_input.Cp * (par - parPrev) / m_dt
              + (par - pd) / m_input.Rp
              - Q;

            R[TAUC] =
              m_input.Cd * (pd - pdPrev) / m_dt
              + (pd - par) / m_input.Rp
              - (m_input.Psv - pd) / m_input.Rd;
          }

        private:
          DenseLinearSystem m_system;
          Input m_input;
          DenseVector m_xPrev;
          DenseVector* m_xCurrent;
          Scalar m_dt;
          Scalar m_time;
      };

    public:
      explicit CCMLC2024(const Input& input)
        : m_input(input),
          m_problem(input),
          m_solver(m_problem),
          m_newton(m_solver)
      {
        m_x.resize(NVAR);
        m_xPrev.resize(NVAR);
        m_x.setZero();
        m_xPrev.setZero();
      }

      CCMLC2024& setAbsoluteTolerance(const Scalar atol)
      {
        m_newton.setAbsoluteTolerance(atol);
        return *this;
      }

      CCMLC2024& setRelativeTolerance(const Scalar rtol)
      {
        m_newton.setRelativeTolerance(rtol);
        return *this;
      }

      CCMLC2024& setStepTolerance(const Scalar stol)
      {
        m_newton.setStepTolerance(stol);
        return *this;
      }

      CCMLC2024& setMaxIterations(const size_t maxIt)
      {
        m_newton.setMaxIterations(maxIt);
        return *this;
      }

      CCMLC2024& setDampingFactor(const Scalar alpha)
      {
        m_newton.setDampingFactor(alpha);
        return *this;
      }

      void initialize(const State& initial)
      {
        m_state = initial;
        m_x[Y] = initial.y;
        m_x[V] = initial.v;
        m_x[PV] = initial.pv;
        m_x[PAR] = initial.par;
        m_x[PD] = initial.pd;
        m_x[EC] = initial.ec;
        m_x[KC] = initial.kc;
        m_x[TAUC] = initial.tauc;
        m_xPrev = m_x;
      }

      Report step(const Scalar dt)
      {
        m_xPrev = m_x;

        m_problem.setPrevious(m_xPrev);
        m_problem.setCurrent(m_x);
        m_problem.setTime(m_state.t);
        m_problem.setTimeStep(dt);

        m_newton.solve(m_x);

        const auto& nr = m_newton.getReport();
        m_report.converged = nr.converged;
        m_report.iterations = nr.iterations;
        m_report.finalResidual = nr.final_residual;
        m_report.finalStepNorm = nr.final_step_norm;
        m_report.reason = nr.reason;

        if (nr.converged)
        {
          m_state.y = m_x[Y];
          m_state.v = m_x[V];
          m_state.pv = m_x[PV];
          m_state.par = m_x[PAR];
          m_state.pd = m_x[PD];
          m_state.ec = m_x[EC];
          m_state.kc = m_x[KC];
          m_state.tauc = m_x[TAUC];
          m_state.t += dt;
        }

        return m_report;
      }

      const State& getState() const noexcept
      {
        return m_state;
      }

      const Report& getReport() const noexcept
      {
        return m_report;
      }

      const DenseVector& getUnknowns() const noexcept
      {
        return m_x;
      }

    private:
      Input m_input;
      State m_state;
      Report m_report;

      DenseVector m_x;
      DenseVector m_xPrev;

      Problem m_problem;
      Solver::LDLT<DenseLinearSystem> m_solver;
      Solver::NewtonSolver<Solver::LDLT<DenseLinearSystem>> m_newton;
  };
}

#endif
