/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_CCMLC2014_H
#define RODIN_HEART_CCMLC2014_CCMLC2014_H

#include <algorithm>
#include <array>
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
#include "Rodin/Heart/CCMLC2014/HolzapfelReducedLaw.h"
#include "Rodin/Heart/CCMLC2014/PassiveLaw.h"

namespace Rodin::Heart
{
  template <
    class PassiveEnergyLaw = HolzapfelReducedLaw<Real>,
    class PassiveLaw = CCMLC2014PassiveLaw<Real>>
  class CCMLC2014T
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
        W,
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
        Scalar w = 1.0;
        Scalar t = 0.0;
      };

      struct Input
      {
        Scalar rho = 1.0;
        Scalar d0 = 1.0;
        Scalar R0 = 1.0;

        Scalar Es = 1.0;
        Scalar mu = 0.0;
        Scalar alpha = 0.0;
        Scalar alphaR = 1.0;
        Scalar k0 = 0.0;
        Scalar sigma0 = 0.0;

        Scalar eta = 0.0;

        Scalar Cp = 1.0;
        Scalar Cd = 1.0;
        Scalar Rp = 1.0;
        Scalar Rd = 1.0;

        Scalar Kat = 0.0;
        Scalar Kp = 0.0;
        Scalar Kar = 0.0;

        std::function<Scalar(Scalar)> u =
          [](Scalar) { return Scalar(0); };

        std::function<Scalar(Scalar)> pAt =
          [](Scalar) { return Scalar(0); };

        std::function<Scalar(Scalar)> pSv =
          [](Scalar) { return Scalar(0); };

        std::function<Scalar(Scalar)> n0 =
          [](Scalar) { return Scalar(1); };

        std::function<Scalar(Scalar)> m0 =
          [](Scalar) { return Scalar(1); };

        std::function<Scalar(Scalar)> m0Prime =
          [](Scalar) { return Scalar(0); };

        PassiveEnergyLaw passiveEnergy;
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
      struct EvalData
      {
        Scalar t = 0.0;
        Scalar dt = 0.0;

        Scalar y = 0.0;
        Scalar v = 0.0;
        Scalar pv = 0.0;
        Scalar par = 0.0;
        Scalar pd = 0.0;
        Scalar ec = 0.0;
        Scalar kc = 0.0;
        Scalar tauc = 0.0;
        Scalar w = 0.0;

        Scalar yPrev = 0.0;
        Scalar vPrev = 0.0;
        Scalar pvPrev = 0.0;
        Scalar parPrev = 0.0;
        Scalar pdPrev = 0.0;
        Scalar ecPrev = 0.0;
        Scalar kcPrev = 0.0;
        Scalar taucPrev = 0.0;
        Scalar wPrev = 0.0;

        Scalar yMid = 0.0;
        Scalar vMid = 0.0;
        Scalar pvMid = 0.0;
        Scalar parMid = 0.0;
        Scalar pdMid = 0.0;
        Scalar ecMid = 0.0;
        Scalar kcMid = 0.0;
        Scalar taucMid = 0.0;
        Scalar wMid = 0.0;

        Scalar uMidInput = 0.0;
        Scalar pAtMid = 0.0;
        Scalar pSvMid = 0.0;
        Scalar n0Frozen = 0.0;
        Scalar m0Mid = 0.0;

        Scalar uPlus = 0.0;
        Scalar uMinus = 0.0;

        Scalar sqrtC = 0.0;
        Scalar C = 0.0;
        Scalar dC_dy = 0.0;
        Scalar e1D = 0.0;
        Scalar de1D_dy = 0.0;
        Scalar Cdot = 0.0;
        Scalar dCdot_dy = 0.0;
        Scalar dCdot_dv = 0.0;

        Scalar sigmaPassive = 0.0;
        Scalar dsigmaPassive_dy = 0.0;

        Scalar sigma1D = 0.0;
        Scalar dsigma1D_dy = 0.0;
        Scalar dsigma1D_dec = 0.0;

        Scalar sigmaSph = 0.0;
        Scalar dsigmaSph_dy = 0.0;
        Scalar dsigmaSph_dv = 0.0;
        Scalar dsigmaSph_dec = 0.0;

        Scalar ecdot = 0.0;
        Scalar absEcdot = 0.0;
        Scalar signEcdot = 0.0;

        Scalar A = 0.0;
        Scalar dA_dec = 0.0;
        Scalar dA_dw = 0.0;

        Scalar Q = 0.0;
        Scalar dQ_dpv = 0.0;
        Scalar dQ_dpar = 0.0;
      };

      static Scalar posPart(Scalar x)
      {
        return x > Scalar(0) ? x : Scalar(0);
      }

      static Scalar negPart(Scalar x)
      {
        return x < Scalar(0) ? -x : Scalar(0);
      }

      static void unpack(const DenseVector& x, EvalData& d)
      {
        d.y = x[Y];
        d.v = x[V];
        d.pv = x[PV];
        d.par = x[PAR];
        d.pd = x[PD];
        d.ec = x[EC];
        d.kc = x[KC];
        d.tauc = x[TAUC];
        d.w = x[W];
      }

      static DenseVector pack(const State& s)
      {
        DenseVector x(NVAR);
        x[Y] = s.y;
        x[V] = s.v;
        x[PV] = s.pv;
        x[PAR] = s.par;
        x[PD] = s.pd;
        x[EC] = s.ec;
        x[KC] = s.kc;
        x[TAUC] = s.tauc;
        x[W] = s.w;
        return x;
      }

      static void prepareEvalData(
          const Input& input,
          const DenseVector& x,
          const State& prev,
          const Scalar t,
          const Scalar dt,
          EvalData& d)
      {
        assert(dt > Scalar(0));
        unpack(x, d);

        d.t = t;
        d.dt = dt;

        d.yPrev = prev.y;
        d.vPrev = prev.v;
        d.pvPrev = prev.pv;
        d.parPrev = prev.par;
        d.pdPrev = prev.pd;
        d.ecPrev = prev.ec;
        d.kcPrev = prev.kc;
        d.taucPrev = prev.tauc;
        d.wPrev = prev.w;

        d.yMid = Scalar(0.5) * (d.y + d.yPrev);
        d.vMid = Scalar(0.5) * (d.v + d.vPrev);
        d.pvMid = Scalar(0.5) * (d.pv + d.pvPrev);
        d.parMid = Scalar(0.5) * (d.par + d.parPrev);
        d.pdMid = Scalar(0.5) * (d.pd + d.pdPrev);
        d.ecMid = Scalar(0.5) * (d.ec + d.ecPrev);
        d.kcMid = Scalar(0.5) * (d.kc + d.kcPrev);
        d.taucMid = Scalar(0.5) * (d.tauc + d.taucPrev);
        d.wMid = Scalar(0.5) * (d.w + d.wPrev);

        const Scalar tMid = t + Scalar(0.5) * dt;
        d.uMidInput = input.u(tMid);
        d.pAtMid = input.pAt(tMid);
        d.pSvMid = input.pSv(tMid);
        d.n0Frozen = input.n0(d.ecPrev);
        d.m0Mid = input.m0(d.ecMid);

        d.uPlus = posPart(d.uMidInput);
        d.uMinus = negPart(d.uMidInput);

        d.sqrtC = Scalar(1) + d.yMid / input.R0;
        assert(d.sqrtC > Scalar(0));
        d.C = d.sqrtC * d.sqrtC;

        d.dC_dy = d.sqrtC / input.R0;
        d.e1D = Scalar(0.5) * (d.C - Scalar(1));
        d.de1D_dy = Scalar(0.5) * d.dC_dy;

        d.Cdot = Scalar(2) * d.sqrtC / input.R0 * d.vMid;
        d.dCdot_dy = d.vMid / (input.R0 * input.R0);
        d.dCdot_dv = d.sqrtC / input.R0;

        PassiveLaw passiveLaw;
        passiveLaw(input.passiveEnergy, d.C, d.dC_dy, d.sigmaPassive, d.dsigmaPassive_dy);

        d.ecdot = (d.ec - d.ecPrev) / dt;
        d.absEcdot = std::abs(d.ecdot);
        d.signEcdot =
          d.ecdot > Scalar(0) ? Scalar(1)
          : (d.ecdot < Scalar(0) ? Scalar(-1) : Scalar(0));

        {
          const Scalar den = Scalar(1) + Scalar(2) * d.ecMid;
          const Scalar den2 = den * den;
          const Scalar den3 = den2 * den;

          d.sigma1D = input.Es * (d.e1D - d.ecMid) / den2;
          d.dsigma1D_dy = input.Es * d.de1D_dy / den2;
          d.dsigma1D_dec =
            Scalar(0.5) * input.Es *
            ( -Scalar(1) / den2 - Scalar(4) * (d.e1D - d.ecMid) / den3 );
        }

        const Scalar viscFactor = Scalar(1) - Scalar(2) * std::pow(d.C, -6);
        const Scalar dViscFactor_dy =
          Scalar(12) * std::pow(d.C, -7) * d.dC_dy;

        const Scalar sigmaViscous = Scalar(2) * input.eta * d.Cdot * viscFactor;
        const Scalar dsigmaViscous_dy =
          Scalar(2) * input.eta *
            (d.dCdot_dy * viscFactor + d.Cdot * dViscFactor_dy);
        const Scalar dsigmaViscous_dv =
          Scalar(2) * input.eta * d.dCdot_dv * viscFactor;

        d.sigmaSph = d.sigma1D + d.sigmaPassive + sigmaViscous;
        d.dsigmaSph_dy = d.dsigma1D_dy + d.dsigmaPassive_dy + dsigmaViscous_dy;
        d.dsigmaSph_dv = dsigmaViscous_dv;
        d.dsigmaSph_dec = d.dsigma1D_dec;

        d.A = d.uPlus + d.wMid * d.uMinus + input.alpha * d.absEcdot;
        d.dA_dec = input.alpha * d.signEcdot / dt;
        d.dA_dw = Scalar(0.5) * d.uMinus;

        if (d.pvMid <= d.pAtMid)
        {
          d.Q = input.Kat * (d.pvMid - d.pAtMid);
          d.dQ_dpv = Scalar(0.5) * input.Kat;
          d.dQ_dpar = Scalar(0);
        }
        else if (d.pvMid <= d.parMid)
        {
          d.Q = input.Kp * (d.pvMid - d.pAtMid);
          d.dQ_dpv = Scalar(0.5) * input.Kp;
          d.dQ_dpar = Scalar(0);
        }
        else
        {
          d.Q = input.Kar * (d.pvMid - d.parMid)
              + input.Kp * (d.parMid - d.pAtMid);
          d.dQ_dpv = Scalar(0.5) * input.Kar;
          d.dQ_dpar = Scalar(0.5) * (-input.Kar + input.Kp);
        }
      }

    public:
      static void evaluateResidual(
          const Input& input,
          const DenseVector& x,
          const State& prev,
          const Scalar t,
          const Scalar dt,
          DenseVector& R)
      {
        EvalData d;
        prepareEvalData(input, x, prev, t, dt, d);

        R.resize(NVAR);
        R.setZero();

        const Scalar geom1 = Scalar(1) + d.yMid / input.R0;
        const Scalar geom2 = geom1 * geom1;

        R[Y] = d.y - d.yPrev - dt * d.vMid;

        R[V] =
          input.rho * input.d0 * (d.v - d.vPrev) / dt
          + (input.d0 / input.R0) * geom1 * d.sigmaSph
          - d.pvMid * geom2;

        R[PV] =
          Scalar(4) * std::numbers::pi_v<Scalar> * input.R0 * input.R0 * geom2 * d.vMid
          - d.Q;

        R[PAR] =
          input.Cp * (d.par - d.parPrev) / dt
          + (d.parMid - d.pdMid) / input.Rp
          - d.Q;

        R[PD] =
          input.Cd * (d.pd - d.pdPrev) / dt
          + (d.pdMid - d.parMid) / input.Rp
          - (d.pSvMid - d.pdMid) / input.Rd;

        {
          const Scalar den = Scalar(1) + Scalar(2) * d.ecMid;
          const Scalar den3 = den * den * den;
          const Scalar rhs =
            input.Es * (d.e1D - d.ecMid) * (Scalar(1) + Scalar(2) * d.e1D) / den3;

          R[EC] = d.taucMid + input.mu * d.ecdot - rhs;
        }

        R[KC] =
          (d.kc - d.kcPrev) / dt
          + d.A * d.kcMid
          - d.n0Frozen * input.k0 * d.uPlus;

        R[TAUC] =
          (d.tauc - d.taucPrev) / dt
          + d.A * d.taucMid
          - d.n0Frozen * input.sigma0 * d.uPlus
          - d.kcMid * d.ecdot;

        R[W] =
          input.alphaR * (d.w - d.wPrev) / dt
          - (d.m0Mid - d.wMid);
      }

      static void evaluateJacobian(
          const Input& input,
          const DenseVector& x,
          const State& prev,
          const Scalar t,
          const Scalar dt,
          DenseMatrix& J)
      {
        EvalData d;
        prepareEvalData(input, x, prev, t, dt, d);

        J.resize(NVAR, NVAR);
        J.setZero();

        const Scalar geom1 = Scalar(1) + d.yMid / input.R0;
        const Scalar geom2 = geom1 * geom1;
        const Scalar dGeom1_dy = Scalar(1) / (Scalar(2) * input.R0);
        const Scalar dGeom2_dy = geom1 / input.R0;

        J(Y, Y) = Scalar(1);
        J(Y, V) = -Scalar(0.5) * dt;

        J(V, Y) =
          (input.d0 / input.R0) * (dGeom1_dy * d.sigmaSph + geom1 * d.dsigmaSph_dy)
          - d.pvMid * dGeom2_dy;
        J(V, V) =
          input.rho * input.d0 / dt
          + (input.d0 / input.R0) * geom1 * d.dsigmaSph_dv;
        J(V, PV) = -Scalar(0.5) * geom2;
        J(V, EC) = (input.d0 / input.R0) * geom1 * d.dsigmaSph_dec;

        J(PV, Y) =
          Scalar(4) * std::numbers::pi_v<Scalar> * input.R0 * input.R0 * dGeom2_dy * d.vMid;
        J(PV, V) =
          Scalar(2) * std::numbers::pi_v<Scalar> * input.R0 * input.R0 * geom2;
        J(PV, PV) = -d.dQ_dpv;
        J(PV, PAR) = -d.dQ_dpar;

        J(PAR, PV) = -d.dQ_dpv;
        J(PAR, PAR) = input.Cp / dt + Scalar(0.5) / input.Rp - d.dQ_dpar;
        J(PAR, PD) = -Scalar(0.5) / input.Rp;

        J(PD, PAR) = -Scalar(0.5) / input.Rp;
        J(PD, PD) = input.Cd / dt + Scalar(0.5) / input.Rp + Scalar(0.5) / input.Rd;

        {
          const Scalar den = Scalar(1) + Scalar(2) * d.ecMid;
          const Scalar den3 = den * den * den;
          const Scalar den4 = den3 * den;
          const Scalar g = (d.e1D - d.ecMid) * (Scalar(1) + Scalar(2) * d.e1D);
          const Scalar dg_dy =
            d.de1D_dy * (Scalar(1) + Scalar(2) * d.e1D)
            + (d.e1D - d.ecMid) * (Scalar(2) * d.de1D_dy);
          const Scalar dg_dec = -Scalar(0.5) * (Scalar(1) + Scalar(2) * d.e1D);

          J(EC, Y) =
            -input.Es * (dg_dy / den3);

          J(EC, EC) =
            input.mu / dt
            - input.Es * (dg_dec / den3 - Scalar(3) * g / den4);

          J(EC, TAUC) = Scalar(0.5);
        }

        J(KC, EC) = d.dA_dec * d.kcMid;
        J(KC, KC) = Scalar(1) / dt + Scalar(0.5) * d.A;
        J(KC, W) = d.dA_dw * d.kcMid;

        J(TAUC, EC) = d.dA_dec * d.taucMid - d.kcMid / dt;
        J(TAUC, KC) = -Scalar(0.5) * d.ecdot;
        J(TAUC, TAUC) = Scalar(1) / dt + Scalar(0.5) * d.A;
        J(TAUC, W) = d.dA_dw * d.taucMid;

        J(W, EC) = -Scalar(0.5) * input.m0Prime(d.ecMid);
        J(W, W) = input.alphaR / dt + Scalar(0.5);
      }

      class Problem final : public Variational::ProblemBase<DenseLinearSystem>
      {
        public:
          using Parent = Variational::ProblemBase<DenseLinearSystem>;
          using ProblemBodyType = typename Parent::ProblemBodyType;

          explicit Problem(const Input& in)
            : m_input(in)
          {
            m_system.getOperator().resize(NVAR, NVAR);
            m_system.getVector().resize(NVAR);
            m_system.getSolution().resize(NVAR);
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

            DenseVector R;
            CCMLC2014T::evaluateResidual(m_input, *m_xCurrent, m_prev, m_time, m_dt, R);
            CCMLC2014T::evaluateJacobian(m_input, *m_xCurrent, m_prev, m_time, m_dt, A);
            b = -R;
            return *this;
          }

          void solve(Solver::LinearSolverBase<DenseLinearSystem>& solver) override
          {
            solver.solve(m_system);
          }

          DenseLinearSystem& getLinearSystem() override { return m_system; }
          const DenseLinearSystem& getLinearSystem() const override { return m_system; }

          Problem* copy() const noexcept override { return new Problem(*this); }

          void setPrevious(const State& prev) { m_prev = prev; }
          void setCurrent(DenseVector& xCurrent) { m_xCurrent = &xCurrent; }
          void setTime(const Scalar t) { m_time = t; }
          void setTimeStep(const Scalar dt) { m_dt = dt; }

        private:
          Input m_input;
          State m_prev;
          DenseVector* m_xCurrent = nullptr;
          DenseLinearSystem m_system;
          Scalar m_time = 0.0;
          Scalar m_dt = 0.0;
      };

    public:
      explicit CCMLC2014T(const Input& input)
        : m_input(input), m_problem(input), m_solver(m_problem), m_newton(m_solver)
      {
        m_x = pack(m_state);
      }

      CCMLC2014T& setAbsoluteTolerance(const Scalar atol)
      {
        m_newton.setAbsoluteTolerance(atol);
        return *this;
      }

      CCMLC2014T& setRelativeTolerance(const Scalar rtol)
      {
        m_newton.setRelativeTolerance(rtol);
        return *this;
      }

      CCMLC2014T& setStepTolerance(const Scalar stol)
      {
        m_newton.setStepTolerance(stol);
        return *this;
      }

      CCMLC2014T& setMaxIterations(const size_t maxIt)
      {
        m_newton.setMaxIterations(maxIt);
        return *this;
      }

      CCMLC2014T& setDampingFactor(const Scalar alpha)
      {
        m_newton.setDampingFactor(alpha);
        return *this;
      }

      void initialize(const State& initial)
      {
        m_state = initial;
        m_x = pack(m_state);
      }

      Report step(const Scalar dt)
      {
        m_x = pack(m_state);

        m_problem.setPrevious(m_state);
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
          m_state.w = m_x[W];
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
      Problem m_problem;
      Solver::LDLT<DenseLinearSystem> m_solver;
      Solver::NewtonSolver<Solver::LDLT<DenseLinearSystem>> m_newton;
  };

  using CCMLC2014 = CCMLC2014T<>;
}

#endif
