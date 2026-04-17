/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_CCMLC2014_H
#define RODIN_HEART_CCMLC2014_CCMLC2014_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>

#include "Rodin/Solver/ForwardDecls.h"
#include "Rodin/Types.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/Problem.h"
#include "Rodin/Solver/PartialPivLU.h"
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
        DISP = 0,
        PV,
        PAR,
        PD,
        NVAR
      };

      struct State
      {
        Scalar y = 0.0;
        Scalar v = 0.0;
        Scalar pv = 0.0;
        Scalar par = 0.0;
        Scalar pd = 0.0;

        // Local active state
        Scalar ec = 0.0;
        Scalar kc = 0.0;
        Scalar tauc = 0.0;
        Scalar gamma = 0.0;
        Scalar beta = 0.0;

        Scalar t = 0.0;
      };

      struct Input
      {
        // Geometry / inertia
        Scalar rho = 1.0;
        Scalar d0 = 1.0;
        Scalar R0 = 1.0;

        // Passive / active material
        Scalar Es = 1.0;
        Scalar eta = 0.0;                // spherical viscous stress coefficient
        Scalar mu = 0.0;                 // DampingParallel
        Scalar alpha = 0.0;              // DestructionRate
        Scalar k0 = 0.0;                 // CrossBridgeStiffness
        Scalar sigma0 = 0.0;             // Contractility

        // Windkessel
        Scalar Cp = 1.0;
        Scalar Cd = 1.0;
        Scalar Rp = 1.0;
        Scalar Rd = 1.0;

        // Valve law
        Scalar Kat = 0.0;                // K_cavities
        Scalar Kp = 0.0;                 // K_closed
        Scalar Kar = 0.0;                // K_artery

        // Tiny p_cav capacity regularization used in valve block
        Scalar cavityCapacity = Scalar(5e-12);

        // Local active solve options
        Scalar localTolerance = Scalar(1e-12);
        size_t localMaxIterations = 50;
        Scalar localDamping = Scalar(1.0);
        Scalar absRegularization = Scalar(1e-14);

        // Initial local active state
        Scalar initFibDef = 0.0;
        Scalar initActiveStiffness = 0.0; // kc(0)
        Scalar initActiveStress = 0.0;    // tauc(0)

        std::function<Scalar(Scalar)> u =
          [](Scalar) { return Scalar(0); };

        std::function<Scalar(Scalar)> pAt =
          [](Scalar) { return Scalar(0); };

        std::function<Scalar(Scalar)> pSv =
          [](Scalar) { return Scalar(0); };

        PassiveEnergyLaw passiveEnergy;
      };

      struct Report
      {
        bool converged = false;
        size_t iterations = 0;
        Scalar finalResidual = 0.0;
        Scalar finalStepNorm = 0.0;
        Solver::NewtonSolver<Solver::PartialPivLU<DenseLinearSystem>>::ConvergenceReason reason =
          Solver::NewtonSolver<Solver::PartialPivLU<DenseLinearSystem>>::ConvergenceReason::MaxIterations;
      };

      struct LocalActiveData
      {
        Scalar fib0 = 0.0;
        Scalar fib1 = 0.0;
        Scalar fib12 = 0.0;

        Scalar gammaOld = 0.0;
        Scalar betaOld = 0.0;
        Scalar gammaNew = 0.0;
        Scalar betaNew = 0.0;

        Scalar u1 = 0.0;
        Scalar u1Plus = 0.0;
        Scalar n0 = 0.0;

        Scalar k21 = 0.0;
        Scalar k22 = 0.0;
        Scalar krc_k22 = 0.0;

        Scalar sigma1d = 0.0;
        Scalar partialSigma1dWrtDisp = 0.0;
        Scalar partialSigma1dWrtEc = 0.0;

        Scalar stressActive = 0.0;
        Scalar diffStressActive = 0.0;

        bool converged = false;
        size_t iterations = 0;
      };

      struct EvalData
      {
        State sn;
        State snm1;

        Scalar tnp1 = 0.0;
        Scalar dt = 0.0;

        // Current Newton iterate (global unknowns)
        Scalar y = 0.0;
        Scalar pv = 0.0;
        Scalar par = 0.0;
        Scalar pd = 0.0;

        // Previous step
        Scalar yPrev = 0.0;
        Scalar pvPrev = 0.0;
        Scalar parPrev = 0.0;
        Scalar pdPrev = 0.0;

        // Previous-previous step (for inertia only)
        Scalar yPrevPrev = 0.0;

        // Midpoint quantities
        Scalar yMid = 0.0;
        Scalar pvMid = 0.0;
        Scalar parMid = 0.0;
        Scalar pdMid = 0.0;
        Scalar vel = 0.0;

        // Geometry at midpoint
        Scalar sqrtC = 0.0;
        Scalar C = 0.0;
        Scalar strain1D = 0.0;
        Scalar diffGreen = 0.0; // d strain1D / d(mid disp)

        // Passive / viscous
        Scalar stressPassive = 0.0;
        Scalar diffStressPassive = 0.0;
        Scalar stressViscous = 0.0;
        Scalar diffStressViscous = 0.0;

        // Flow branches
        Scalar pAtCur = 0.0;
        Scalar pAtPrev = 0.0;
        Scalar pSvMid = 0.0;

        Scalar cavityFluxCur = 0.0;
        Scalar cavityFluxPrev = 0.0;
        Scalar dCavityFluxCur_dPv = 0.0;
        Scalar dCavityFluxCur_dPar = 0.0;

        Scalar windkesselOutflow = 0.0;
        Scalar dWindkesselOutflow_dPv = 0.0;
        Scalar dWindkesselOutflow_dPar = 0.0;

        LocalActiveData active;
      };

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
            CCMLC2014T::evaluateDynamicResidual(m_input, *m_xCurrent, m_sn, m_snm1, m_time, m_dt, R);
            CCMLC2014T::evaluateDynamicJacobian(m_input, *m_xCurrent, m_sn, m_snm1, m_time, m_dt, A);
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

          void setCurrent(DenseVector& xCurrent) { m_xCurrent = &xCurrent; }

          void setStepData(const State& sn, const State& snm1, Scalar time, Scalar dt)
          {
            m_sn = sn;
            m_snm1 = snm1;
            m_time = time;
            m_dt = dt;
          }

        private:
          Input m_input;
          DenseVector* m_xCurrent = nullptr;
          DenseLinearSystem m_system;

          State m_sn;
          State m_snm1;
          Scalar m_time = 0.0;
          Scalar m_dt = 0.0;
      };

    private:
      static DenseVector packUnknowns(const State& s)
      {
        DenseVector x(NVAR);
        x[DISP] = s.y;
        x[PV]   = s.pv;
        x[PAR]  = s.par;
        x[PD]   = s.pd;
        return x;
      }

      static State unpackUnknownsIntoState(const DenseVector& x, const State& base, Scalar t)
      {
        State s = base;
        s.y = x[DISP];
        s.pv = x[PV];
        s.par = x[PAR];
        s.pd = x[PD];
        s.t = t;
        return s;
      }

      static Scalar n0_piecewise(Scalar fib0)
      {
        const Scalar x1 = Scalar(-0.4);
        const Scalar y1 = Scalar(0.0);
        const Scalar x2 = Scalar(0.3);
        const Scalar y2 = Scalar(0.38);
        const Scalar x3 = Scalar(0.73);
        const Scalar y3 = Scalar(0.74);
        const Scalar x4 = Scalar(1.0);
        const Scalar y4 = Scalar(1.0);
        const Scalar x5 = Scalar(1.3);
        const Scalar y5 = Scalar(1.0);
        const Scalar x6 = Scalar(2.4);
        const Scalar y6 = Scalar(0.0);

        Scalar n0 = Scalar(0.0);

        if (fib0 < x2)
          n0 = ((y2 - y1) / (x2 - x1)) * (fib0 - x2) + y2;
        else if (fib0 < x3)
          n0 = ((y3 - y2) / (x3 - x2)) * (fib0 - x3) + y3;
        else if (fib0 < x4)
          n0 = ((y4 - y3) / (x4 - x3)) * (fib0 - x4) + y4;
        else if (fib0 < x5)
          n0 = y4;
        else if (fib0 < x6)
          n0 = ((y6 - y5) / (x6 - x5)) * (fib0 - x6) + y6;

        return std::max<Scalar>(n0, Scalar(0));
      }

      static void computeMidpointKinematics(const Input& in, EvalData& d)
      {
        d.yMid = Scalar(0.5) * (d.y + d.yPrev);
        d.pvMid = Scalar(0.5) * (d.pv + d.pvPrev);
        d.parMid = Scalar(0.5) * (d.par + d.parPrev);
        d.pdMid = Scalar(0.5) * (d.pd + d.pdPrev);
        d.vel = (d.y - d.yPrev) / d.dt;

        d.sqrtC = Scalar(1) + d.yMid / in.R0;
        assert(d.sqrtC > Scalar(0));
        d.C = d.sqrtC * d.sqrtC;

        d.strain1D = Scalar(0.5) * (d.C - Scalar(1));
        d.diffGreen = d.sqrtC / in.R0;
      }

      static void computePassiveContribution(const Input& in, EvalData& d)
      {
        Scalar dC_dyMid = Scalar(2) * d.sqrtC / in.R0;

        Scalar sigmaPassive = 0.0;
        Scalar dsigmaPassive_dyMid = 0.0;
        PassiveLaw passiveLaw;
        passiveLaw(in.passiveEnergy, d.C, dC_dyMid, sigmaPassive, dsigmaPassive_dyMid);

        d.stressPassive = sigmaPassive;
        d.diffStressPassive = Scalar(0.5) * dsigmaPassive_dyMid;
      }

      static void computeViscousContribution(const Input& in, EvalData& d)
      {
        const Scalar R0 = in.R0;
        const Scalar nu = in.eta;
        const Scalar vel = d.vel;
        const Scalar sqrtC = d.sqrtC;
        const Scalar C = d.C;

        const Scalar diffGreen_rr = Scalar(-2) / R0 * std::pow(sqrtC, -5);
        const Scalar diffGreen_pp = d.diffGreen;

        const Scalar stress_rr = nu * diffGreen_rr * vel;
        const Scalar stress_pp = nu * diffGreen_pp * vel;
        const Scalar stress_tt = stress_pp;

        d.stressViscous =
          stress_pp + stress_tt - Scalar(2) * std::pow(C, -3) * stress_rr;

        const Scalar d_dotC_rr_dy =
          Scalar(10) / (R0 * R0) * nu * std::pow(sqrtC, -6) * vel
          + Scalar(2) * nu * diffGreen_rr / d.dt;

        const Scalar d_dotC_pp_dy =
          nu / (R0 * R0) * vel
          + Scalar(2) * nu * diffGreen_pp / d.dt;

        const Scalar diffStress =
          d_dotC_pp_dy + d_dotC_pp_dy
          - Scalar(2) * (std::pow(C, -3) * d_dotC_rr_dy
              - Scalar(6) / R0 * std::pow(sqrtC, -7) * stress_rr);

        d.diffStressViscous = diffStress;
      }

      static void updateInternalVariables0D(
          const Input& in,
          Scalar dt,
          Scalar fib0,
          Scalar fib1,
          Scalar gammaOld,
          Scalar betaOld,
          Scalar u1,
          Scalar& gammaNew,
          Scalar& betaNew,
          Scalar& n0)
      {
        const Scalar alpha = in.alpha;
        const Scalar k0 = in.k0;
        const Scalar sigma0 = in.sigma0;

        const Scalar u1Plus = std::max<Scalar>(u1, Scalar(0));

        n0 = n0_piecewise(fib0);

        const Scalar denominatorGamma =
          Scalar(1) + dt * (std::abs(u1) + alpha * std::abs(fib1 - fib0) / dt);

        Scalar gammaSquare =
          (gammaOld * gammaOld + dt * n0 * k0 * u1Plus) / denominatorGamma;
        gammaSquare = std::max<Scalar>(Scalar(1e-16), gammaSquare);
        gammaNew = std::sqrt(gammaSquare);

        const Scalar denominatorBeta =
          Scalar(1)
          + Scalar(0.5) * dt * n0 * k0 * u1Plus / gammaSquare
          + Scalar(0.5) * dt * (std::abs(u1) + alpha * std::abs(fib1 - fib0) / dt);

        betaNew =
          (betaOld + gammaNew * (fib1 - fib0) + dt * n0 * sigma0 * u1Plus / gammaNew)
          / denominatorBeta;
      }

      static void getPartialDerivativesInternalVariablesWrtFibDeformation(
          const Input& in,
          Scalar dt,
          Scalar gammaOld,
          Scalar betaOld,
          Scalar u1,
          Scalar fib0,
          Scalar fib1,
          Scalar& derBetaGammaWrtFib)
      {
        const Scalar alpha = in.alpha;
        const Scalar k0 = in.k0;
        const Scalar sigma0 = in.sigma0;

        const Scalar deltaFib = fib1 - fib0;
        const Scalar absDelta = std::abs(deltaFib);
        const Scalar sign =
          (deltaFib > Scalar(0)) ? Scalar(1) :
          (deltaFib < Scalar(0)) ? Scalar(-1) : Scalar(0);

        const Scalar u1Plus = std::max<Scalar>(u1, Scalar(0));
        const Scalar n0 = n0_piecewise(fib0);

        const Scalar Dg =
          Scalar(1) + dt * std::abs(u1) + alpha * absDelta;

        const Scalar Ng =
          gammaOld * gammaOld + dt * n0 * k0 * u1Plus;

        const Scalar gammaNewSq = std::max<Scalar>(Scalar(1e-16), Ng / Dg);
        const Scalar gammaNew = std::sqrt(gammaNewSq);

        const Scalar dGammaSq =
          -Ng * alpha * sign / (Dg * Dg);

        const Scalar dGamma =
          Scalar(0.5) * dGammaSq / gammaNew;

        const Scalar Nb =
          betaOld + gammaNew * deltaFib + dt * n0 * sigma0 * u1Plus / gammaNew;

        const Scalar Db =
          Scalar(1)
          + Scalar(0.5) * dt * n0 * k0 * u1Plus / gammaNewSq
          + Scalar(0.5) * dt * std::abs(u1)
          + Scalar(0.5) * alpha * absDelta;

        const Scalar dNb =
          dGamma * deltaFib
          + gammaNew
          - dt * n0 * sigma0 * u1Plus * dGamma / (gammaNew * gammaNew);

        const Scalar dDb =
          -Scalar(0.5) * dt * n0 * k0 * u1Plus * dGammaSq / (gammaNewSq * gammaNewSq)
          + Scalar(0.5) * alpha * sign;

        const Scalar dBeta =
          (dNb * Db - Nb * dDb) / (Db * Db);

        derBetaGammaWrtFib =
          dGamma * (Nb / Db) + gammaNew * dBeta;
      }

      static bool solveLocalDynamicActive(
          const Input& in,
          const EvalData& d,
          LocalActiveData& a)
      {
        a.fib0 = d.sn.ec;
        a.gammaOld = d.sn.gamma;
        a.betaOld = d.sn.beta;
        a.u1 = in.u(d.tnp1);
        a.u1Plus = std::max<Scalar>(a.u1, Scalar(0));
        a.n0 = n0_piecewise(a.fib0);

        a.fib1 = a.fib0;

        bool ok = false;
        for (size_t it = 0; it < in.localMaxIterations; ++it)
        {
          a.iterations = it + 1;
          a.fib12 = Scalar(0.5) * (a.fib1 + a.fib0);

          updateInternalVariables0D(
              in, d.dt, a.fib0, a.fib1, a.gammaOld, a.betaOld, a.u1,
              a.gammaNew, a.betaNew, a.n0);

          a.k21 = -in.Es * (Scalar(1) + Scalar(4) * d.strain1D - Scalar(2) * a.fib12);

          {
            Scalar derBetaGammaWrtFib = 0.0;
            getPartialDerivativesInternalVariablesWrtFibDeformation(
                in, d.dt, a.gammaOld, a.betaOld, a.u1, a.fib0, a.fib1, derBetaGammaWrtFib);

            a.k22 =
                Scalar(3) * std::pow(Scalar(1) + Scalar(2) * a.fib12, 2) *
                    (a.gammaNew * a.betaNew + in.mu * (a.fib1 - a.fib0) / d.dt)
              + std::pow(Scalar(1) + Scalar(2) * a.fib12, 3) *
                    (derBetaGammaWrtFib + in.mu / d.dt)
              + Scalar(0.5) * in.Es * (Scalar(1) + Scalar(2) * d.strain1D);
          }

          const Scalar Rraw =
            (a.gammaNew * a.betaNew + in.mu * (a.fib1 - a.fib0) / d.dt)
              * std::pow(Scalar(1) + Scalar(2) * a.fib12, 3)
            - in.Es * (d.strain1D - a.fib12) * (Scalar(1) + Scalar(2) * d.strain1D);

          a.krc_k22 = Rraw / a.k22;

          if (std::abs(Rraw) < in.localTolerance)
          {
            ok = true;
            break;
          }

          if (std::abs(a.k22) < std::numeric_limits<Scalar>::epsilon())
            break;

          a.fib1 += -in.localDamping * a.krc_k22;
        }

        a.fib12 = Scalar(0.5) * (a.fib1 + a.fib0);
        updateInternalVariables0D(
            in, d.dt, a.fib0, a.fib1, a.gammaOld, a.betaOld, a.u1,
            a.gammaNew, a.betaNew, a.n0);

        a.sigma1d =
          in.Es / std::pow(Scalar(1) + Scalar(2) * a.fib12, 2)
          * (d.strain1D - a.fib12);

        a.partialSigma1dWrtDisp =
          in.Es / std::pow(Scalar(1) + Scalar(2) * a.fib12, 2);

        a.partialSigma1dWrtEc =
          in.Es / std::pow(Scalar(1) + Scalar(2) * a.fib12, 3)
          * (Scalar(2) * a.fib12 - Scalar(4) * d.strain1D - Scalar(1));

        const Scalar coefSchurD2W =
          Scalar(0.5) * a.partialSigma1dWrtEc / a.k22 * a.k21;

        const Scalar tangentCorrection =
          a.partialSigma1dWrtDisp - coefSchurD2W;

        const Scalar rhsCorrection =
          Scalar(0.5) * a.partialSigma1dWrtEc * a.krc_k22;

        a.stressActive = a.sigma1d - rhsCorrection;
        a.diffStressActive = tangentCorrection * d.diffGreen;
        a.converged = ok;
        return ok;
      }

      static void buildEvalData(
          const Input& in,
          const DenseVector& x,
          const State& sn,
          const State& snm1,
          Scalar tnp1,
          Scalar dt,
          EvalData& d)
      {
        d.sn = sn;
        d.snm1 = snm1;

        d.tnp1 = tnp1;
        d.dt = dt;

        d.y = x[DISP];
        d.pv = x[PV];
        d.par = x[PAR];
        d.pd = x[PD];

        d.yPrev = sn.y;
        d.pvPrev = sn.pv;
        d.parPrev = sn.par;
        d.pdPrev = sn.pd;
        d.yPrevPrev = snm1.y;

        d.pAtCur = in.pAt(tnp1);
        d.pAtPrev = in.pAt(sn.t);
        d.pSvMid = in.pSv(sn.t + Scalar(0.5) * dt);

        computeMidpointKinematics(in, d);
        computePassiveContribution(in, d);
        computeViscousContribution(in, d);
        solveLocalDynamicActive(in, d, d.active);

        {
          const bool mitralOpenCur = d.pv <= d.pAtCur;
          const bool bothClosedCur = d.pAtCur <= d.pv && d.pv <= d.par;

          if (mitralOpenCur)
          {
            d.cavityFluxCur = in.Kat * (d.pv - d.pAtCur);
            d.dCavityFluxCur_dPv = in.Kat;
            d.dCavityFluxCur_dPar = Scalar(0);
          }
          else if (bothClosedCur)
          {
            d.cavityFluxCur = in.Kp * (d.pv - d.pAtCur);
            d.dCavityFluxCur_dPv = in.Kp;
            d.dCavityFluxCur_dPar = Scalar(0);
          }
          else
          {
            d.cavityFluxCur = in.Kar * (d.pv - d.par) + in.Kp * (d.par - d.pAtCur);
            d.dCavityFluxCur_dPv = in.Kar;
            d.dCavityFluxCur_dPar = -in.Kar + in.Kp;
          }
        }

        {
          const bool mitralOpenPrev = d.pvPrev <= d.pAtPrev;
          const bool bothClosedPrev = d.pAtPrev <= d.pvPrev && d.pvPrev <= d.parPrev;

          if (mitralOpenPrev)
          {
            d.cavityFluxPrev = in.Kat * (d.pvPrev - d.pAtPrev);
          }
          else if (bothClosedPrev)
          {
            d.cavityFluxPrev = in.Kp * (d.pvPrev - d.pAtPrev);
          }
          else
          {
            d.cavityFluxPrev = in.Kar * (d.pvPrev - d.parPrev) + in.Kp * (d.parPrev - d.pAtPrev);
          }
        }

        {
          const Scalar flowAr = in.Kar * (d.pvMid - d.parMid);
          d.windkesselOutflow = (flowAr > Scalar(0)) ? flowAr : Scalar(0);

          const Scalar difFlow = (flowAr > Scalar(0)) ? (Scalar(0.5) * in.Kar) : Scalar(0);
          d.dWindkesselOutflow_dPv = difFlow;
          d.dWindkesselOutflow_dPar = -difFlow;
        }
      }

    public:
      static void evaluateDynamicResidual(
          const Input& in,
          const DenseVector& x,
          const State& sn,
          const State& snm1,
          Scalar tnp1,
          Scalar dt,
          DenseVector& R)
      {
        EvalData d;
        buildEvalData(in, x, sn, snm1, tnp1, dt, d);

        R.resize(NVAR);
        R.setZero();

        const Scalar coeff = in.d0 * in.rho;
        const Scalar geom = Scalar(1) + d.yMid / in.R0;
        const Scalar geom2 = geom * geom;

        const Scalar totalStress =
          d.stressPassive + d.stressViscous + d.active.stressActive;

        R[DISP] =
          coeff / (dt * dt) * (d.y - Scalar(2) * d.yPrev + d.yPrevPrev)
          + in.d0 / in.R0 * geom * totalStress
          - d.pvMid * geom2;

        {
          const Scalar capacityCur = in.cavityCapacity / dt * d.pv;
          const Scalar capacityPrev = in.cavityCapacity / dt * d.pvPrev;
          const Scalar volumeTerm =
            Scalar(4) * std::numbers::pi_v<Scalar> * in.R0 * in.R0 * geom2 * d.vel;

          R[PV] =
            (capacityCur - capacityPrev)
            + Scalar(0.5) * (d.cavityFluxCur + d.cavityFluxPrev)
            + volumeTerm;
        }

        {
          R[PAR] =
            in.Cp / dt * (d.par - d.parPrev)
            + (d.parMid - d.pdMid) / in.Rp
            - d.windkesselOutflow;
        }

        {
          R[PD] =
            in.Cd / dt * (d.pd - d.pdPrev)
            + (d.pdMid - d.parMid) / in.Rp
            - (d.pSvMid - d.pdMid) / in.Rd;
        }
      }

      static void evaluateDynamicJacobian(
          const Input& in,
          const DenseVector& x,
          const State& sn,
          const State& snm1,
          Scalar tnp1,
          Scalar dt,
          DenseMatrix& J)
      {
        EvalData d;
        buildEvalData(in, x, sn, snm1, tnp1, dt, d);

        J.resize(NVAR, NVAR);
        J.setZero();

        const Scalar coeff = in.d0 * in.rho;
        const Scalar geom = Scalar(1) + d.yMid / in.R0;
        const Scalar geom2 = geom * geom;

        const Scalar totalStress =
          d.stressPassive + d.stressViscous + d.active.stressActive;
        const Scalar totalDiffStress =
          d.diffStressPassive + d.diffStressViscous + d.active.diffStressActive;

        J(DISP, DISP) +=
          coeff / (dt * dt)
          + Scalar(0.5) * in.d0 / (in.R0 * in.R0) * totalStress
          + Scalar(0.5) * in.d0 / in.R0 * geom * totalDiffStress
          - Scalar(1) / in.R0 * d.pvMid * geom;
        J(DISP, PV) += -Scalar(0.5) * geom2;

        J(PV, DISP) +=
          Scalar(4) * std::numbers::pi_v<Scalar> * in.R0 * geom * d.vel
          + Scalar(4) * std::numbers::pi_v<Scalar> * in.R0 * in.R0 * geom2 * (Scalar(1) / dt);

        J(PV, PV) += in.cavityCapacity / dt + Scalar(0.5) * d.dCavityFluxCur_dPv;
        J(PV, PAR) += Scalar(0.5) * d.dCavityFluxCur_dPar;

        J(PAR, PV) += -d.dWindkesselOutflow_dPv;
        J(PAR, PAR) += in.Cp / dt + Scalar(1) / (Scalar(2) * in.Rp) - d.dWindkesselOutflow_dPar;
        J(PAR, PD) += -Scalar(1) / (Scalar(2) * in.Rp);

        J(PD, PAR) += -Scalar(1) / (Scalar(2) * in.Rp);
        J(PD, PD) += in.Cd / dt + Scalar(1) / (Scalar(2) * in.Rp) + Scalar(1) / (Scalar(2) * in.Rd);
      }

      explicit CCMLC2014T(const Input& input)
        : m_input(input), m_problem(input), m_solver(m_problem), m_newton(m_solver)
      {
        m_x = packUnknowns(m_state);
      }

      CCMLC2014T& setAbsoluteTolerance(const Scalar atol)
      {
        m_atol = atol;
        m_newton.setAbsoluteTolerance(atol);
        return *this;
      }

      CCMLC2014T& setRelativeTolerance(const Scalar rtol)
      {
        m_rtol = rtol;
        m_newton.setRelativeTolerance(rtol);
        return *this;
      }

      CCMLC2014T& setStepTolerance(const Scalar stol)
      {
        m_stol = stol;
        m_newton.setStepTolerance(stol);
        return *this;
      }

      CCMLC2014T& setMaxIterations(const size_t maxIt)
      {
        m_maxIt = maxIt;
        m_newton.setMaxIterations(maxIt);
        return *this;
      }

      CCMLC2014T& setDampingFactor(const Scalar alpha)
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

        m_prevState = m_state;
        m_x = packUnknowns(m_state);
      }

      Report step(const Scalar dt)
      {
        assert(dt > Scalar(0));

        const State sn = m_state;
        const State snm1 = m_prevState;

        m_x = packUnknowns(sn);

        m_problem.setCurrent(m_x);
        m_problem.setStepData(sn, snm1, sn.t + dt, dt);

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
          buildEvalData(m_input, m_x, sn, snm1, sn.t + dt, dt, d);

          m_prevState = m_state;
          m_state = unpackUnknownsIntoState(m_x, m_state, sn.t + dt);
          m_state.v = (m_state.y - sn.y) / dt;
          m_state.ec = d.active.fib1;
          m_state.gamma = d.active.gammaNew;
          m_state.beta = d.active.betaNew;
          m_state.kc = m_state.gamma * m_state.gamma;
          m_state.tauc = m_state.gamma * m_state.beta;
        }

        return m_report;
      }

      const State& getState() const noexcept { return m_state; }
      const Report& getReport() const noexcept { return m_report; }
      const DenseVector& getUnknowns() const noexcept { return m_x; }

    private:
      Input m_input;
      State m_state;
      State m_prevState;

      Scalar m_atol = Scalar(1e-10);
      Scalar m_rtol = Scalar(1e-10);
      Scalar m_stol = Scalar(1e-12);
      size_t m_maxIt = 50;
      Scalar m_damping = Scalar(1.0);

      Report m_report;
      DenseVector m_x;

      Problem m_problem;
      Solver::PartialPivLU<DenseLinearSystem> m_solver;
      Solver::NewtonSolver<Solver::PartialPivLU<DenseLinearSystem>> m_newton;
  };

  using CCMLC2014 = CCMLC2014T<>;
}

#endif
