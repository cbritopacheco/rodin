/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_CCMLC2014_H
#define RODIN_HEART_CCMLC2014_CCMLC2014_H

#include <cassert>
#include <cmath>
#include <functional>
#include <numbers>

#include "Rodin/Types.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/Problem.h"
#include "Rodin/Solver/LDLT.h"
#include "Rodin/Solver/NewtonSolver.h"

#include "PassiveLaw.h"

namespace Rodin::Heart
{
  template <class PassiveLaw = CCMLC2014Laws::HolzapfelOgdenLaw>
  class CCMLC2014
  {
    public:
      using Scalar = Real;
      using DenseMatrix = Math::Matrix<Scalar>;
      using DenseVector = Math::Vector<Scalar>;
      using DenseLinearSystem = Math::LinearSystem<DenseMatrix, DenseVector>;

      enum Variable : size_t { Y = 0, V, PV, PAR, PD, EC, KC, TAUC, NVAR };

      struct State
      {
        Scalar y = 0.0, v = 0.0, pv = 0.0, par = 0.0, pd = 0.0, ec = 0.0, kc = 0.0, tauc = 0.0, t = 0.0;
      };

      struct Input
      {
        Scalar rho = 1.0, d0 = 1.0, R0 = 1.0, eta = 0.0;
        Scalar Es = 1.0, mu = 1.0, alpha = 0.0, n0 = 0.0, k0 = 0.0, sigma0 = 0.0, w = 1.0;
        Scalar Cp = 1.0, Cd = 1.0, Rp = 1.0, Rd = 1.0, Psv = 0.0;

        PassiveLaw passive;

        std::function<Scalar(Scalar)> e1D = [](Scalar C) { return Scalar(0.5) * (C * C - Scalar(1)); };
        std::function<Scalar(Scalar)> e1DPrime = [](Scalar C) { return C; };

        std::function<Scalar(Scalar)> ubar = [](Scalar) { return Scalar(0); };
        std::function<Scalar(Scalar)> pAt = [](Scalar) { return Scalar(0); };

        std::function<Scalar(Scalar, Scalar, Scalar)> valveFlow =
          [](Scalar pv, Scalar par, Scalar) { return pv > par ? pv - par : Scalar(0); };
        std::function<Scalar(Scalar, Scalar, Scalar)> dValveFlow_dPv =
          [](Scalar pv, Scalar par, Scalar) { return pv > par ? Scalar(1) : Scalar(0); };
        std::function<Scalar(Scalar, Scalar, Scalar)> dValveFlow_dPar =
          [](Scalar pv, Scalar par, Scalar) { return pv > par ? Scalar(-1) : Scalar(0); };
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

      static void evaluateResidual(const Input& input, const DenseVector& x, const DenseVector& xPrev, Scalar t, Scalar dt, DenseVector& R)
      {
        const Scalar y = x[Y], v = x[V], pv = x[PV], par = x[PAR], pd = x[PD], ec = x[EC], kc = x[KC], tauc = x[TAUC];
        const Scalar yPrev = xPrev[Y], vPrev = xPrev[V], parPrev = xPrev[PAR], pdPrev = xPrev[PD], ecPrev = xPrev[EC], kcPrev = xPrev[KC], taucPrev = xPrev[TAUC];

        const Scalar C = 1.0 + y / input.R0;
        const Scalar invC = 1.0 / C;
        const Scalar C2 = C * C;
        const Scalar C6inv = invC * invC * invC * invC * invC * invC;
        const Scalar cdot = v / input.R0;

        const Scalar j1 = C2 + 2.0 * invC;
        const Scalar j2 = 2.0 * C + invC * invC;
        const Scalar j4 = C2;

        const auto passive = input.passive(C, j1, j2, j4);
        const Scalar e1d = input.e1D(C);
        const Scalar ecDot = (ec - ecPrev) / dt;

        const Scalar denom1 = 1.0 + 2.0 * ec;
        const Scalar sigma1D = input.Es * (e1d - ec) / (denom1 * denom1);
        const Scalar sigmaSph =
          sigma1D
          + 4.0 * (1.0 - invC * invC * invC) * (passive.dW1 + C * passive.dW2)
          + 2.0 * passive.dW4
          + 2.0 * input.eta * cdot * (1.0 - 2.0 * C6inv);

        const Scalar pat = input.pAt(t + dt);
        const Scalar Q = input.valveFlow(pv, par, pat);
        const Scalar ubar = input.ubar(t + dt);
        const Scalar A = posPart(ubar) + input.w * negPart(ubar) + input.alpha * std::abs(ecDot);

        const Scalar volumeRate = 4.0 * std::numbers::pi_v<Scalar> * input.R0 * input.R0 * C2 * v;

        R.resize(NVAR);
        R[Y] = y - yPrev - dt * v;
        R[V] = input.rho * input.d0 * (v - vPrev) / dt + (input.d0 / input.R0) * C * sigmaSph - pv * C2;
        R[PV] = tauc + input.mu * ecDot - input.Es * ((e1d - ec) * (1.0 + 2.0 * e1d)) / (denom1 * denom1 * denom1);
        R[PAR] = (kc - kcPrev) / dt + A * kc - input.n0 * input.k0 * posPart(ubar);
        R[PD] = (tauc - taucPrev) / dt + A * tauc - input.n0 * input.sigma0 * posPart(ubar) - kc * ecDot;
        R[EC] = volumeRate + Q;
        R[KC] = input.Cp * (par - parPrev) / dt + (par - pd) / input.Rp - Q;
        R[TAUC] = input.Cd * (pd - pdPrev) / dt + (pd - par) / input.Rp - (input.Psv - pd) / input.Rd;
      }

      static void evaluateJacobian(const Input& input, const DenseVector& x, const DenseVector& xPrev, Scalar t, Scalar dt, DenseMatrix& J)
      {
        const Scalar y = x[Y], v = x[V], pv = x[PV], par = x[PAR], ec = x[EC], kc = x[KC], tauc = x[TAUC];
        const Scalar ecPrev = xPrev[EC];

        const Scalar C = 1.0 + y / input.R0;
        const Scalar C2 = C * C;
        const Scalar C6inv = std::pow(C, -6);
        const Scalar C7inv = std::pow(C, -7);
        const Scalar cdot = v / input.R0;
        const Scalar dCdy = 1.0 / input.R0;

        const Scalar j1 = C2 + 2.0 / C;
        const Scalar j2 = 2.0 * C + 1.0 / (C * C);
        const Scalar j4 = C2;
        const auto passive = input.passive(C, j1, j2, j4);

        const Scalar e1d = input.e1D(C);
        const Scalar de1d_dy = input.e1DPrime(C) * dCdy;
        const Scalar ecDot = (ec - ecPrev) / dt;
        const Scalar signEcDot = (ecDot > 0) ? 1.0 : ((ecDot < 0) ? -1.0 : 0.0);

        const Scalar denom1 = 1.0 + 2.0 * ec;
        const Scalar denom2 = denom1 * denom1;
        const Scalar denom3 = denom2 * denom1;
        const Scalar denom4 = denom3 * denom1;

        const Scalar sigma1D = input.Es * (e1d - ec) / denom2;
        const Scalar dsigma1D_dy = input.Es * de1d_dy / denom2;
        const Scalar dsigma1D_dec = input.Es * (-1.0 / denom2 - 4.0 * (e1d - ec) / denom3);

        const Scalar sigmaPassive = 4.0 * (1.0 - std::pow(C, -3)) * (passive.dW1 + C * passive.dW2) + 2.0 * passive.dW4;
        const Scalar dsigmaPassive_dy = passive.dSigmaPassive_dC * dCdy;

        const Scalar sigmaViscous = 2.0 * input.eta * cdot * (1.0 - 2.0 * C6inv);
        const Scalar dsigmaViscous_dv = 2.0 * input.eta / input.R0 * (1.0 - 2.0 * C6inv);
        const Scalar dsigmaViscous_dy = 24.0 * input.eta * cdot * C7inv * dCdy;

        const Scalar sigmaSph = sigma1D + sigmaPassive + sigmaViscous;
        const Scalar dsigma_dy = dsigma1D_dy + dsigmaPassive_dy + dsigmaViscous_dy;
        const Scalar dsigma_dv = dsigmaViscous_dv;
        const Scalar dsigma_dec = dsigma1D_dec;

        const Scalar pat = input.pAt(t + dt);
        const Scalar Q_pv = input.dValveFlow_dPv(pv, par, pat);
        const Scalar Q_par = input.dValveFlow_dPar(pv, par, pat);

        const Scalar ubar = input.ubar(t + dt);
        const Scalar A = posPart(ubar) + input.w * negPart(ubar) + input.alpha * std::abs(ecDot);
        const Scalar dA_dec = input.alpha * signEcDot / dt;

        const Scalar g = (e1d - ec) * (1.0 + 2.0 * e1d);
        const Scalar dg_dec = -(1.0 + 2.0 * e1d);
        const Scalar dg_de1d = 1.0 + 4.0 * e1d - 2.0 * ec;

        J.resize(NVAR, NVAR);
        J.setZero();

        J(Y, Y) = 1.0;
        J(Y, V) = -dt;

        J(V, Y) = (input.d0 / input.R0) * (dCdy * sigmaSph + C * dsigma_dy) - 2.0 * pv * C * dCdy;
        J(V, V) = input.rho * input.d0 / dt + (input.d0 / input.R0) * C * dsigma_dv;
        J(V, PV) = -C2;
        J(V, EC) = (input.d0 / input.R0) * C * dsigma_dec;

        J(PV, Y) = -input.Es * (dg_de1d * de1d_dy) / denom3;
        J(PV, EC) = input.mu / dt - input.Es * (dg_dec / denom3 - 6.0 * g / denom4);
        J(PV, TAUC) = 1.0;

        J(PAR, EC) = dA_dec * kc;
        J(PAR, KC) = 1.0 / dt + A;

        J(PD, EC) = dA_dec * tauc - kc / dt;
        J(PD, KC) = -ecDot;
        J(PD, TAUC) = 1.0 / dt + A;

        const Scalar K = 4.0 * std::numbers::pi_v<Scalar> * input.R0 * input.R0;
        J(EC, Y) = K * 2.0 * C * dCdy * v;
        J(EC, V) = K * C2;
        J(EC, PV) = Q_pv;
        J(EC, PAR) = Q_par;

        J(KC, PV) = -Q_pv;
        J(KC, PAR) = input.Cp / dt + 1.0 / input.Rp - Q_par;
        J(KC, PD) = -1.0 / input.Rp;

        J(TAUC, PAR) = -1.0 / input.Rp;
        J(TAUC, PD) = input.Cd / dt + 1.0 / input.Rp + 1.0 / input.Rd;
      }

    private:
      static Scalar posPart(Scalar x) { return (x > 0) ? x : 0; }
      static Scalar negPart(Scalar x) { return (x < 0) ? -x : 0; }

      class Problem final : public Variational::ProblemBase<DenseLinearSystem>
      {
        public:
          using Parent = Variational::ProblemBase<DenseLinearSystem>;
          using ProblemBodyType = typename Parent::ProblemBodyType;

          explicit Problem(const Input& in) : m_input(in)
          {
            m_system.getOperator().resize(NVAR, NVAR);
            m_system.getVector().resize(NVAR);
            m_system.getSolution().resize(NVAR);
          }

          Parent& operator=(const ProblemBodyType&) override { return *this; }

          Problem& assemble() override
          {
            assert(m_xCurrent);
            auto& A = m_system.getOperator();
            auto& b = m_system.getVector();
            m_system.getSolution().setZero();

            DenseVector R;
            CCMLC2014::evaluateResidual(m_input, *m_xCurrent, m_xPrev, m_time, m_dt, R);
            CCMLC2014::evaluateJacobian(m_input, *m_xCurrent, m_xPrev, m_time, m_dt, A);
            b = -R;
            return *this;
          }

          void solve(Solver::LinearSolverBase<DenseLinearSystem>& solver) override { solver.solve(m_system); }
          DenseLinearSystem& getLinearSystem() override { return m_system; }
          const DenseLinearSystem& getLinearSystem() const override { return m_system; }
          Problem* copy() const noexcept override { return new Problem(*this); }

          void setPrevious(const DenseVector& xPrev) { m_xPrev = xPrev; }
          void setCurrent(DenseVector& xCurrent) { m_xCurrent = &xCurrent; }
          void setTime(Scalar t) { m_time = t; }
          void setTimeStep(Scalar dt) { m_dt = dt; }

        private:
          Input m_input;
          DenseVector m_xPrev;
          DenseVector* m_xCurrent = nullptr;
          DenseLinearSystem m_system;
          Scalar m_time = 0.0;
          Scalar m_dt = 0.0;
      };

    public:
      explicit CCMLC2014(const Input& input)
        : m_input(input), m_problem(input), m_solver(m_problem), m_newton(m_solver)
      {
        m_x.resize(NVAR);
        m_xPrev.resize(NVAR);
        m_x.setZero();
        m_xPrev.setZero();
      }

      CCMLC2014& setAbsoluteTolerance(Scalar v){ m_newton.setAbsoluteTolerance(v); return *this; }
      CCMLC2014& setRelativeTolerance(Scalar v){ m_newton.setRelativeTolerance(v); return *this; }
      CCMLC2014& setStepTolerance(Scalar v){ m_newton.setStepTolerance(v); return *this; }
      CCMLC2014& setMaxIterations(size_t v){ m_newton.setMaxIterations(v); return *this; }
      CCMLC2014& setDampingFactor(Scalar v){ m_newton.setDampingFactor(v); return *this; }

      void initialize(const State& initial)
      {
        m_state = initial;
        m_x[Y]=initial.y; m_x[V]=initial.v; m_x[PV]=initial.pv; m_x[PAR]=initial.par;
        m_x[PD]=initial.pd; m_x[EC]=initial.ec; m_x[KC]=initial.kc; m_x[TAUC]=initial.tauc;
        m_xPrev = m_x;
      }

      Report step(Scalar dt)
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
          m_state.y=m_x[Y]; m_state.v=m_x[V]; m_state.pv=m_x[PV]; m_state.par=m_x[PAR];
          m_state.pd=m_x[PD]; m_state.ec=m_x[EC]; m_state.kc=m_x[KC]; m_state.tauc=m_x[TAUC];
          m_state.t += dt;
        }
        return m_report;
      }

      const State& getState() const noexcept { return m_state; }
      const Report& getReport() const noexcept { return m_report; }
      const DenseVector& getUnknowns() const noexcept { return m_x; }

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
