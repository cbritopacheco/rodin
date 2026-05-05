/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file DynamicSystem.h
 * @brief Residual and analytic Jacobian assembly for the CCMLC2014 0D system.
 *
 * The active variables are solved as part of the global Newton system:
 * @f[
 *   (y, v, p_v, p_{ar}, p_d, e_c, k_c, \tau_c, w).
 * @f]
 * This avoids the previous local condensation of @f$ e_c @f$ and keeps the
 * 0D wall, active, pressure, and relaxation equations on the same implicit
 * time grid.
 */
#ifndef RODIN_HEART_CCMLC2014_NUMERICS_DYNAMICSYSTEM_H
#define RODIN_HEART_CCMLC2014_NUMERICS_DYNAMICSYSTEM_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

#include "Rodin/Heart/CCMLC2014/Model/State.h"
#include "Rodin/Heart/CCMLC2014/Physics/Windkessel.h"

namespace Rodin::Heart::CCMLC2014::Numerics
{
  /**
   * @brief Assembles the coupled full-state residual and exact tangent.
   *
   * @tparam PassiveLaw Passive stress operator type.
   * @tparam Input Model input parameter type.
   */
  template <class PassiveLaw, class Input>
  class DynamicSystem
  {
    public:
      /**
       * @brief Construct the dynamic system assembler.
       */
      explicit DynamicSystem(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Build intermediate evaluation data from candidate unknowns.
       */
      template <class DenseVector, class StateType, class EvalData>
      void buildEvalData(
          const DenseVector& candidateUnknowns,
          const StateType& currentState,
          const StateType& previousState,
          typename DenseVector::Scalar tnp1,
          typename DenseVector::Scalar dt,
          EvalData& evalData) const
      {
        evalData.sn = currentState;
        evalData.snm1 = previousState;

        evalData.tnp1 = tnp1;
        evalData.dt = dt;

        evalData.y = candidateUnknowns[Model::RadialDisplacement];
        evalData.v = candidateUnknowns[Model::RadialVelocity];
        evalData.pv = candidateUnknowns[Model::VentricularPressure];
        evalData.par = candidateUnknowns[Model::ArterialPressure];
        evalData.pd = candidateUnknowns[Model::DistalPressure];
        evalData.ec = candidateUnknowns[Model::FiberDeformation];
        evalData.kc = candidateUnknowns[Model::ActiveStiffness];
        evalData.tauc = candidateUnknowns[Model::ActiveStress];
        evalData.w = candidateUnknowns[Model::LoadDependentRelaxation];

        evalData.yPrev = currentState.y;
        evalData.vPrev = currentState.v;
        evalData.pvPrev = currentState.pv;
        evalData.parPrev = currentState.par;
        evalData.pdPrev = currentState.pd;
        evalData.yPrevPrev = previousState.y;

        evalData.pAtCur = m_input.pAt(tnp1);
        evalData.pAtPrev = m_input.pAt(currentState.t);
        evalData.pSvMid = m_input.pSv(tnp1);

        evaluateKinematicsAndStresses(evalData);
        evaluateActiveDiagnostics(evalData);
        evaluateValveDiagnostics(evalData);
      }

      /**
       * @brief Evaluate the coupled 0D residual vector.
       */
      template <class DenseVector, class EvalData>
      void evaluateResidual(
          const EvalData& evalData,
          DenseVector& residualVector) const
      {
        using Scalar = typename DenseVector::Scalar;
        residualVector.resize(Model::NumberOfVariables);
        residualVector.setZero();

        const auto tc = getTimeCoefficients(evalData);
        const Scalar yDot =
          timeDerivative(evalData.y, evalData.sn.y, evalData.snm1.y, tc);
        const Scalar vDot =
          timeDerivative(evalData.v, evalData.sn.v, evalData.snm1.v, tc);
        const Scalar pvDot =
          timeDerivative(evalData.pv, evalData.sn.pv, evalData.snm1.pv, tc);
        const Scalar parDot =
          timeDerivative(evalData.par, evalData.sn.par, evalData.snm1.par, tc);
        const Scalar pdDot =
          timeDerivative(evalData.pd, evalData.sn.pd, evalData.snm1.pd, tc);
        const Scalar ecDot =
          timeDerivative(evalData.ec, evalData.sn.ec, evalData.snm1.ec, tc);
        const Scalar kcDot =
          timeDerivative(evalData.kc, evalData.sn.kc, evalData.snm1.kc, tc);
        const Scalar taucDot =
          timeDerivative(evalData.tauc, evalData.sn.tauc, evalData.snm1.tauc, tc);
        const Scalar wDot =
          timeDerivative(evalData.w, evalData.sn.w, evalData.snm1.w, tc);

        const Scalar radius = m_input.R0 + evalData.y;
        const Scalar geometricStretch = evalData.sqrtC;
        const Scalar totalStress =
          evalData.stressPassive + evalData.stressViscous
          + evalData.active.activeStress;

        residualVector[Model::RadialDisplacement] =
          yDot - evalData.v;

        residualVector[Model::RadialVelocity] =
          m_input.d0 * m_input.rho * vDot
          + m_input.d0 / m_input.R0 * geometricStretch * totalStress
          - evalData.pv * geometricStretch * geometricStretch;

        residualVector[Model::VentricularPressure] =
          m_input.cavityCapacity * pvDot
          + evalData.cavityFluxCur
          + Scalar(4) * std::numbers::pi_v<Scalar>
            * radius * radius * evalData.v;

        residualVector[Model::ArterialPressure] =
          m_input.Cp * parDot
          + evalData.windkesselflowP;

        residualVector[Model::DistalPressure] =
          m_input.Cd * pdDot
          + evalData.windkesselflowD;

        {
          const Scalar h = activeStretch(evalData.ec);
          const Scalar h3 = h * h * h;
          residualVector[Model::FiberDeformation] =
            (evalData.tauc + m_input.mu * ecDot) * h3
            - m_input.Es
              * (evalData.strain1D - evalData.ec)
              * (Scalar(1) + Scalar(2) * evalData.strain1D);
        }

        const auto active = evaluateActiveRates(evalData, ecDot);
        residualVector[Model::ActiveStiffness] =
          kcDot + active.rate * evalData.kc
          - active.recruitment * m_input.k0 * active.uPositive;

        residualVector[Model::ActiveStress] =
          taucDot - evalData.kc * ecDot
          - active.recruitment * m_input.sigma0 * active.uPositive
          + active.rate * evalData.tauc;

        if (m_input.alphaR <= std::numeric_limits<Scalar>::epsilon())
        {
          residualVector[Model::LoadDependentRelaxation] =
            evalData.w - m_input.m0(evalData.ec);
        }
        else
        {
          residualVector[Model::LoadDependentRelaxation] =
            wDot - (m_input.m0(evalData.ec) - evalData.w) / m_input.alphaR;
        }
      }

      /**
       * @brief Assemble the exact Jacobian of the coupled residual.
       */
      template <class DenseMatrix, class EvalData>
      void evaluateJacobian(
          const EvalData& evalData,
          DenseMatrix& jacobianMatrix,
          typename DenseMatrix::Scalar) const
      {
        using Scalar = typename DenseMatrix::Scalar;
        jacobianMatrix.resize(Model::NumberOfVariables, Model::NumberOfVariables);
        jacobianMatrix.setZero();

        const auto tc = getTimeCoefficients(evalData);
        const Scalar a0 = tc.current;

        const Scalar radius = m_input.R0 + evalData.y;
        const Scalar s = evalData.sqrtC;
        const Scalar dsdy = Scalar(1) / m_input.R0;
        const Scalar dEdy = evalData.diffGreen;
        const Scalar F = Scalar(1) + Scalar(2) * evalData.strain1D;

        const Scalar totalStress =
          evalData.stressPassive + evalData.stressViscous
          + evalData.active.activeStress;

        // Kinematic equation: D(y) - v = 0
        jacobianMatrix(Model::RadialDisplacement, Model::RadialDisplacement) = a0;
        jacobianMatrix(Model::RadialDisplacement, Model::RadialVelocity) = -Scalar(1);

        // Wall momentum equation.
        jacobianMatrix(Model::RadialVelocity, Model::RadialDisplacement) =
          m_input.d0 / m_input.R0
          * (dsdy * totalStress
             + s * (evalData.diffStressPassive
                    + evalData.diffStressViscous
                    + evalData.active.dActiveStressWrtDisplacement))
          - Scalar(2) * evalData.pv * s * dsdy;

        jacobianMatrix(Model::RadialVelocity, Model::RadialVelocity) =
          m_input.d0 * m_input.rho * a0
          + m_input.d0 / m_input.R0
            * s * evalData.diffStressViscousWrtVelocity;

        jacobianMatrix(Model::RadialVelocity, Model::VentricularPressure) =
          -s * s;

        jacobianMatrix(Model::RadialVelocity, Model::FiberDeformation) =
          m_input.d0 / m_input.R0
          * s * evalData.active.partialActiveStressWrtFiberDeformation;

        // Cavity pressure balance.
        jacobianMatrix(Model::VentricularPressure, Model::RadialDisplacement) =
          Scalar(8) * std::numbers::pi_v<Scalar> * radius * evalData.v;

        jacobianMatrix(Model::VentricularPressure, Model::RadialVelocity) =
          Scalar(4) * std::numbers::pi_v<Scalar> * radius * radius;

        jacobianMatrix(Model::VentricularPressure, Model::VentricularPressure) =
          m_input.cavityCapacity * a0 + evalData.dCavityFluxCur_dPv;

        jacobianMatrix(Model::VentricularPressure, Model::ArterialPressure) =
          evalData.dCavityFluxCur_dPar;

        // Proximal Windkessel balance.
        jacobianMatrix(Model::ArterialPressure, Model::VentricularPressure) =
          -evalData.dWindkesselOutflow_dPv;

        jacobianMatrix(Model::ArterialPressure, Model::ArterialPressure) +=
          m_input.Cp * a0 + evalData.dWindkesselflowP_dPar;

        jacobianMatrix(Model::ArterialPressure, Model::DistalPressure) +=
          evalData.dWindkesselflowP_dPd;

        // --- Row: DistalPressure ---
        jacobianMatrix(Model::DistalPressure, Model::ArterialPressure) +=
          evalData.dWindkesselflowD_dPar;

        jacobianMatrix(Model::DistalPressure, Model::DistalPressure) +=
          m_input.Cd * a0
          + evalData.dWindkesselflowD_dPd;

        // Fiber-deformation equilibrium.
        const Scalar ecDot =
          timeDerivative(evalData.ec, evalData.sn.ec, evalData.snm1.ec, tc);
        const Scalar h = activeStretch(evalData.ec);
        const Scalar h2 = h * h;
        const Scalar h3 = h2 * h;
        const Scalar activeRateStress = evalData.tauc + m_input.mu * ecDot;

        jacobianMatrix(Model::FiberDeformation, Model::RadialDisplacement) =
          -m_input.Es * dEdy * (F + Scalar(2) * (evalData.strain1D - evalData.ec));

        jacobianMatrix(Model::FiberDeformation, Model::FiberDeformation) =
          m_input.mu * a0 * h3
          + Scalar(6) * activeRateStress * h2
          + m_input.Es * F;

        jacobianMatrix(Model::FiberDeformation, Model::ActiveStress) = h3;

        const auto active = evaluateActiveRates(evalData, ecDot);
        const Scalar dRateDec = active.dRateDEcDot * a0;
        const Scalar dn0Dec = active.dRecruitmentDEc;

        // Active stiffness evolution.
        jacobianMatrix(Model::ActiveStiffness, Model::FiberDeformation) =
          evalData.kc * dRateDec
          - dn0Dec * m_input.k0 * active.uPositive;

        jacobianMatrix(Model::ActiveStiffness, Model::ActiveStiffness) =
          a0 + active.rate;

        jacobianMatrix(Model::ActiveStiffness, Model::LoadDependentRelaxation) =
          evalData.kc * active.uNegative;

        // Active stress evolution.
        jacobianMatrix(Model::ActiveStress, Model::FiberDeformation) =
          -evalData.kc * a0
          - dn0Dec * m_input.sigma0 * active.uPositive
          + evalData.tauc * dRateDec;

        jacobianMatrix(Model::ActiveStress, Model::ActiveStiffness) =
          -ecDot;

        jacobianMatrix(Model::ActiveStress, Model::ActiveStress) =
          a0 + active.rate;

        jacobianMatrix(Model::ActiveStress, Model::LoadDependentRelaxation) =
          evalData.tauc * active.uNegative;

        // Load-dependent relaxation.
        if (m_input.alphaR <= std::numeric_limits<Scalar>::epsilon())
        {
          jacobianMatrix(
              Model::LoadDependentRelaxation,
              Model::FiberDeformation) = -m_input.dm0(evalData.ec);
          jacobianMatrix(
              Model::LoadDependentRelaxation,
              Model::LoadDependentRelaxation) = Scalar(1);
        }
        else
        {
          jacobianMatrix(
              Model::LoadDependentRelaxation,
              Model::FiberDeformation) =
            -m_input.dm0(evalData.ec) / m_input.alphaR;
          jacobianMatrix(
              Model::LoadDependentRelaxation,
              Model::LoadDependentRelaxation) =
            a0 + Scalar(1) / m_input.alphaR;
        }
      }

    private:
      template <class Scalar>
      struct TimeCoefficients
      {
        Scalar current = 0.0;
        Scalar previous = 0.0;
        Scalar previousPrevious = 0.0;
      };

      template <class Scalar>
      struct RecruitmentData
      {
        Scalar value = 0.0;
        Scalar derivative = 0.0;
      };

      template <class Scalar>
      struct ActiveRateData
      {
        Scalar uPositive = 0.0;
        Scalar uNegative = 0.0;
        Scalar recruitment = 0.0;
        Scalar dRecruitmentDEc = 0.0;
        Scalar rate = 0.0;
        Scalar dRateDEcDot = 0.0;
      };

      template <class EvalData>
      TimeCoefficients<decltype(std::declval<EvalData>().y)>
      getTimeCoefficients(const EvalData& data) const
      {
        using Scalar = decltype(data.y);
        if (m_input.timeScheme == Model::TimeScheme::BDF2
            && data.sn.t > data.snm1.t + Scalar(0.5) * data.dt)
        {
          return {
            Scalar(1.5) / data.dt,
            -Scalar(2) / data.dt,
            Scalar(0.5) / data.dt
          };
        }
        return {
          Scalar(1) / data.dt,
          -Scalar(1) / data.dt,
          Scalar(0)
        };
      }

      template <class Scalar>
      static Scalar timeDerivative(
          Scalar current,
          Scalar previous,
          Scalar previousPrevious,
          const TimeCoefficients<Scalar>& tc)
      {
        return tc.current * current
             + tc.previous * previous
             + tc.previousPrevious * previousPrevious;
      }

      template <class Scalar>
      static Scalar activeStretch(Scalar ec)
      {
        return Scalar(1) + Scalar(2) * ec;
      }

      template <class Scalar>
      Scalar regularizedAbs(Scalar x) const
      {
        const Scalar eps = m_input.absRegularization;
        return std::sqrt(x * x + eps * eps);
      }

      template <class Scalar>
      Scalar regularizedAbsDerivative(Scalar x) const
      {
        const Scalar ax = regularizedAbs(x);
        if (ax <= std::numeric_limits<Scalar>::epsilon())
          return Scalar(0);
        return x / ax;
      }

      template <class Scalar>
      static RecruitmentData<Scalar> computeRecruitment(Scalar ec)
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

        RecruitmentData<Scalar> result;
        if (ec < x2)
        {
          result.derivative = (y2 - y1) / (x2 - x1);
          result.value = result.derivative * (ec - x2) + y2;
        }
        else if (ec < x3)
        {
          result.derivative = (y3 - y2) / (x3 - x2);
          result.value = result.derivative * (ec - x3) + y3;
        }
        else if (ec < x4)
        {
          result.derivative = (y4 - y3) / (x4 - x3);
          result.value = result.derivative * (ec - x4) + y4;
        }
        else if (ec < x5)
        {
          result.value = y4;
          result.derivative = Scalar(0);
        }
        else if (ec < x6)
        {
          result.derivative = (y6 - y5) / (x6 - x5);
          result.value = result.derivative * (ec - x6) + y6;
        }

        if (result.value <= Scalar(0))
        {
          result.value = Scalar(0);
          result.derivative = Scalar(0);
        }
        return result;
      }

      template <class EvalData>
      ActiveRateData<decltype(std::declval<EvalData>().y)>
      evaluateActiveRates(const EvalData& data, decltype(std::declval<EvalData>().y) ecDot) const
      {
        using Scalar = decltype(data.y);
        ActiveRateData<Scalar> result;
        const Scalar u = m_input.u(data.tnp1);
        result.uPositive = std::max<Scalar>(u, Scalar(0));
        result.uNegative = std::max<Scalar>(-u, Scalar(0));

        const auto recruitment = computeRecruitment(data.ec);
        result.recruitment = recruitment.value;
        result.dRecruitmentDEc = recruitment.derivative;

        result.rate =
          result.uPositive + data.w * result.uNegative
          + m_input.alpha * regularizedAbs(ecDot);
        result.dRateDEcDot =
          m_input.alpha * regularizedAbsDerivative(ecDot);
        return result;
      }

      template <class EvalData>
      void evaluateKinematicsAndStresses(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        data.yMid = data.y;
        data.pvMid = data.pv;
        data.parMid = data.par;
        data.pdMid = data.pd;
        data.vel = data.v;

        data.sqrtC = Scalar(1) + data.y / m_input.R0;
        data.C = data.sqrtC * data.sqrtC;
        data.strain1D = Scalar(0.5) * (data.C - Scalar(1));
        data.diffGreen = data.sqrtC / m_input.R0;

        Scalar passiveStress = 0.0;
        Scalar passiveStressDerivativeWrtDisplacement = 0.0;
        PassiveLaw passiveLaw;
        passiveLaw(
            m_input.passiveEnergy,
            data.C,
            Scalar(2) * data.sqrtC / m_input.R0,
            passiveStress,
            passiveStressDerivativeWrtDisplacement);

        data.stressPassive = passiveStress;
        data.diffStressPassive = passiveStressDerivativeWrtDisplacement;

        const Scalar s = data.sqrtC;
        const Scalar eta = m_input.eta;
        data.stressViscous =
          eta * data.v / m_input.R0
          * (Scalar(2) * s + Scalar(4) * std::pow(s, Scalar(-11)));
        data.diffStressViscous =
          eta * data.v / (m_input.R0 * m_input.R0)
          * (Scalar(2) - Scalar(44) * std::pow(s, Scalar(-12)));
        data.diffStressViscousWrtVelocity =
          eta / m_input.R0
          * (Scalar(2) * s + Scalar(4) * std::pow(s, Scalar(-11)));

        const Scalar h = activeStretch(data.ec);
        const Scalar h2 = h * h;
        const Scalar h3 = h2 * h;
        data.active.activeStressOneDimensional =
          m_input.Es / h2 * (data.strain1D - data.ec);
        data.active.partialActiveStressWrtDisplacement =
          m_input.Es / h2 * data.diffGreen;
        data.active.partialActiveStressWrtFiberDeformation =
          m_input.Es
          * (-Scalar(1) / h2
             - Scalar(4) * (data.strain1D - data.ec) / h3);

        data.active.activeStress = data.active.activeStressOneDimensional;
        data.active.dActiveStressWrtDisplacement =
          data.active.partialActiveStressWrtDisplacement;
      }

      template <class EvalData>
      void evaluateActiveDiagnostics(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        const auto tc = getTimeCoefficients(data);
        const Scalar ecDot =
          timeDerivative(data.ec, data.sn.ec, data.snm1.ec, tc);
        const auto active = evaluateActiveRates(data, ecDot);

        data.active.fiberDeformationPrevious = data.sn.ec;
        data.active.fiberDeformationCurrent = data.ec;
        data.active.fiberDeformationMidpoint = data.ec;
        data.active.gammaPrevious = data.sn.gamma;
        data.active.betaPrevious = data.sn.beta;
        data.active.wPrevious = data.sn.w;
        data.active.gammaCurrent =
          std::sqrt(std::max<Scalar>(data.kc, Scalar(0)));
        data.active.betaCurrent =
          (data.active.gammaCurrent > Scalar(0))
          ? data.tauc / data.active.gammaCurrent
          : Scalar(0);
        data.active.wCurrent = data.w;
        data.active.activationDrive = m_input.u(data.tnp1);
        data.active.activationDrivePositivePart = active.uPositive;
        data.active.activationDriveNegativePart = active.uNegative;
        data.active.relaxationTarget = m_input.m0(data.ec);
        data.active.relaxationDrive = active.rate;
        data.active.recruitmentFraction = active.recruitment;
        data.active.converged = true;
        data.active.iterations = 0;
      }

      template <class EvalData>
      void evaluateValveDiagnostics(EvalData& data) const
      {
        using Scalar = decltype(data.y);

        const bool mitralOpenCurrent = data.pv <= data.pAtCur;
        const bool bothClosedCurrent =
          data.pAtCur <= data.pv && data.pv <= data.par;

        if (mitralOpenCurrent)
        {
          data.cavityFluxCur = m_input.Kat * (data.pv - data.pAtCur);
          data.dCavityFluxCur_dPv = m_input.Kat;
          data.dCavityFluxCur_dPar = Scalar(0);
        }
        else if (bothClosedCurrent)
        {
          data.cavityFluxCur = m_input.Kp * (data.pv - data.pAtCur);
          data.dCavityFluxCur_dPv = m_input.Kp;
          data.dCavityFluxCur_dPar = Scalar(0);
        }
        else
        {
          data.cavityFluxCur =
            m_input.Kar * (data.pv - data.par)
            + m_input.Kp * (data.par - data.pAtCur);
          data.dCavityFluxCur_dPv = m_input.Kar;
          data.dCavityFluxCur_dPar = -m_input.Kar + m_input.Kp;
        }

        data.cavityFluxPrev = Scalar(0);

        if (m_input.windkesselRheology == Model::WindkesselRheology::CarreauYasuda)
        {
          Physics::WindkesselOutflowEvaluator<
            Input,
            Physics::Rheology::NonNewtonian::CY> windkessel(m_input);
          windkessel.evaluate(data);
        }
        else
        {
          Physics::WindkesselOutflowEvaluator<
            Input,
            Physics::Rheology::Newtonian> windkessel(m_input);
          windkessel.evaluate(data);
        }
      }

      const Input& m_input;
  };
}

#endif
