/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file DynamicSystem.h
 * @brief Residual and analytic Jacobian assembly for the 0D poroelastic
 * sphere.
 *
 * The unknowns are solved as one global Newton system:
 * @f[
 *   (y, \Phi, p_v, p_{ar}, p_d, e_c, k_c, \tau_c, w).
 * @f]
 * The wall is quasi-static — the pressure-volume law of the @f$ v_1 @f$
 * projection is an algebraic equation in @f$ y @f$ — and the fluid mass
 * balance @f$ V_{w0}\dot\Phi = Q @f$ with the lumped two-resistor perfusion
 * source @f$ Q = \gamma_{ar}(p_{ar} - \tilde p) - \gamma_{ven}(\tilde p - p_{sv}) @f$
 * is the differential equation of the porosity. The active, valve and
 * Windkessel equations are those of CCMLC2014 on the same implicit time grid.
 * The coronary inflow @f$ \gamma_{ar}(p_{ar} - \tilde p) @f$ is drawn from the
 * proximal Windkessel compartment so that fluid is conserved.
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_NUMERICS_DYNAMICSYSTEM_H
#define RODIN_HEART_POROELASTICSPHERE_NUMERICS_DYNAMICSYSTEM_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

#include "Rodin/Heart/CCMLC2014/Physics/Windkessel.h"
#include "Rodin/Heart/PoroelasticSphere/Model/State.h"
#include "Rodin/Heart/PoroelasticSphere/Physics/Kinematics.h"
#include "Rodin/Heart/PoroelasticSphere/Physics/WallStress.h"

namespace Rodin::Heart::PoroelasticSphere::Numerics
{
  /**
   * @brief Assembles the coupled full-state residual and exact tangent.
   *
   * @tparam PassiveLaw Passive stress operator type (radial Piola components).
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
        : m_input(input), m_kinematics(input), m_wall(input)
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
        evalData.phi = candidateUnknowns[Model::Porosity];
        evalData.pv = candidateUnknowns[Model::VentricularPressure];
        evalData.par = candidateUnknowns[Model::ArterialPressure];
        evalData.pd = candidateUnknowns[Model::DistalPressure];
        evalData.ec = candidateUnknowns[Model::FiberDeformation];
        evalData.kc = candidateUnknowns[Model::ActiveStiffness];
        evalData.tauc = candidateUnknowns[Model::ActiveStress];
        evalData.w = candidateUnknowns[Model::LoadDependentRelaxation];

        evalData.pAtCur = m_input.pAt(tnp1);
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

        // Wall equilibrium (v1 projection): P_v(r_in, Phi) - p_v = 0.
        residualVector[Model::RadialDisplacement] = evalData.wallPressure - evalData.pv;

        // Fluid mass balance: V_w0 dPhi/dt = Q_in - Q_out.
        residualVector[Model::Porosity] = evalData.phiDot -
          (evalData.perfusionInflow - evalData.perfusionOutflow) / evalData.Vw0;

        residualVector[Model::VentricularPressure] = m_input.cavityCapacity * pvDot +
          evalData.cavityFluxCur +
          Scalar(4) * std::numbers::pi_v<Scalar> * evalData.rin * evalData.rin * evalData.yDot;

        residualVector[Model::ArterialPressure] =
          m_input.Cp * parDot + evalData.windkesselflowP + evalData.perfusionInflow;

        residualVector[Model::DistalPressure] =
          m_input.Cd * pdDot + evalData.windkesselflowD;

        {
          const Scalar h = activeStretch(evalData.ec);
          const Scalar h3 = h * h * h;
          residualVector[Model::FiberDeformation] =
            (evalData.tauc + m_input.mu * ecDot) * h3 -
            m_input.Es * (evalData.strain1D - evalData.ec) *
              (Scalar(1) + Scalar(2) * evalData.strain1D);
        }

        const auto active = evaluateActiveRates(evalData, ecDot);
        residualVector[Model::ActiveStiffness] = kcDot + active.rate * evalData.kc -
          active.recruitment * m_input.k0 * active.uPositive;

        residualVector[Model::ActiveStress] = taucDot - evalData.kc * ecDot -
          active.recruitment * m_input.sigma0 * active.uPositive +
          active.rate * evalData.tauc;

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
      void evaluateJacobian(const EvalData& evalData, DenseMatrix& jacobianMatrix,
        typename DenseMatrix::Scalar) const
      {
        using Scalar = typename DenseMatrix::Scalar;
        jacobianMatrix.resize(Model::NumberOfVariables, Model::NumberOfVariables);
        jacobianMatrix.setZero();

        const auto tc = getTimeCoefficients(evalData);
        const Scalar a0 = tc.current;

        const Scalar dEdy = evalData.diffGreen;
        const Scalar dEdphi = evalData.diffGreenWrtPorosity;
        const Scalar F = Scalar(1) + Scalar(2) * evalData.strain1D;

        // Wall equilibrium.
        jacobianMatrix(Model::RadialDisplacement, Model::RadialDisplacement) =
          evalData.dWallPressure_dy;
        jacobianMatrix(Model::RadialDisplacement, Model::Porosity) =
          evalData.dWallPressure_dphi;
        jacobianMatrix(Model::RadialDisplacement, Model::VentricularPressure) = -Scalar(1);
        jacobianMatrix(Model::RadialDisplacement, Model::FiberDeformation) =
          evalData.dWallPressure_dec;

        // Fluid mass balance: d(Q_in - Q_out)/dpf = -(gamma_ar + gamma_ven).
        const Scalar gSum = (m_input.gammaAr + m_input.gammaVen) / evalData.Vw0;
        jacobianMatrix(Model::Porosity, Model::RadialDisplacement) = gSum * evalData.dPf_dy;
        jacobianMatrix(Model::Porosity, Model::Porosity) = a0 + gSum * evalData.dPf_dphi;
        jacobianMatrix(Model::Porosity, Model::ArterialPressure) =
          -m_input.gammaAr / evalData.Vw0;
        jacobianMatrix(Model::Porosity, Model::FiberDeformation) = gSum * evalData.dPf_dec;

        // Cavity pressure balance.
        jacobianMatrix(Model::VentricularPressure, Model::RadialDisplacement) =
          Scalar(4) * std::numbers::pi_v<Scalar> * evalData.rin *
          (Scalar(2) * evalData.yDot + evalData.rin * a0);

        jacobianMatrix(Model::VentricularPressure, Model::VentricularPressure) =
          m_input.cavityCapacity * a0 + evalData.dCavityFluxCur_dPv;

        jacobianMatrix(Model::VentricularPressure, Model::ArterialPressure) =
          evalData.dCavityFluxCur_dPar;

        // Proximal Windkessel balance (with the coronary inflow drawn from it).
        jacobianMatrix(Model::ArterialPressure, Model::VentricularPressure) =
          -evalData.dWindkesselOutflow_dPv;

        jacobianMatrix(Model::ArterialPressure, Model::ArterialPressure) +=
          m_input.Cp * a0 + evalData.dWindkesselflowP_dPar + m_input.gammaAr;

        jacobianMatrix(Model::ArterialPressure, Model::DistalPressure) +=
          evalData.dWindkesselflowP_dPd;

        jacobianMatrix(Model::ArterialPressure, Model::RadialDisplacement) =
          -m_input.gammaAr * evalData.dPf_dy;
        jacobianMatrix(Model::ArterialPressure, Model::Porosity) =
          -m_input.gammaAr * evalData.dPf_dphi;
        jacobianMatrix(Model::ArterialPressure, Model::FiberDeformation) =
          -m_input.gammaAr * evalData.dPf_dec;

        // --- Row: DistalPressure ---
        jacobianMatrix(Model::DistalPressure, Model::ArterialPressure) +=
          evalData.dWindkesselflowD_dPar;

        jacobianMatrix(Model::DistalPressure, Model::DistalPressure) +=
          m_input.Cd * a0 + evalData.dWindkesselflowD_dPd;

        // Fiber-deformation equilibrium.
        const Scalar ecDot =
          timeDerivative(evalData.ec, evalData.sn.ec, evalData.snm1.ec, tc);
        const Scalar h = activeStretch(evalData.ec);
        const Scalar h2 = h * h;
        const Scalar h3 = h2 * h;
        const Scalar activeRateStress = evalData.tauc + m_input.mu * ecDot;

        jacobianMatrix(Model::FiberDeformation, Model::RadialDisplacement) =
          -m_input.Es * dEdy * (F + Scalar(2) * (evalData.strain1D - evalData.ec));

        jacobianMatrix(Model::FiberDeformation, Model::Porosity) =
          -m_input.Es * dEdphi * (F + Scalar(2) * (evalData.strain1D - evalData.ec));

        jacobianMatrix(Model::FiberDeformation, Model::FiberDeformation) =
          m_input.mu * a0 * h3 + Scalar(6) * activeRateStress * h2 + m_input.Es * F;

        jacobianMatrix(Model::FiberDeformation, Model::ActiveStress) = h3;

        const auto active = evaluateActiveRates(evalData, ecDot);
        const Scalar dRateDec = active.dRateDEcDot * a0;
        const Scalar dn0Dec = active.dRecruitmentDEc;

        // Active stiffness evolution.
        jacobianMatrix(Model::ActiveStiffness, Model::FiberDeformation) =
          evalData.kc * dRateDec - dn0Dec * m_input.k0 * active.uPositive;

        jacobianMatrix(Model::ActiveStiffness, Model::ActiveStiffness) = a0 + active.rate;

        jacobianMatrix(Model::ActiveStiffness, Model::LoadDependentRelaxation) =
          evalData.kc * active.uNegative;

        // Active stress evolution.
        jacobianMatrix(Model::ActiveStress, Model::FiberDeformation) = -evalData.kc * a0 -
          dn0Dec * m_input.sigma0 * active.uPositive + evalData.tauc * dRateDec;

        jacobianMatrix(Model::ActiveStress, Model::ActiveStiffness) = -ecDot;

        jacobianMatrix(Model::ActiveStress, Model::ActiveStress) = a0 + active.rate;

        jacobianMatrix(Model::ActiveStress, Model::LoadDependentRelaxation) =
          evalData.tauc * active.uNegative;

        // Load-dependent relaxation.
        if (m_input.alphaR <= std::numeric_limits<Scalar>::epsilon())
        {
          jacobianMatrix(Model::LoadDependentRelaxation, Model::FiberDeformation) =
            -m_input.dm0(evalData.ec);
          jacobianMatrix(Model::LoadDependentRelaxation, Model::LoadDependentRelaxation) =
            Scalar(1);
        }
        else
        {
          jacobianMatrix(Model::LoadDependentRelaxation, Model::FiberDeformation) =
            -m_input.dm0(evalData.ec) / m_input.alphaR;
          jacobianMatrix(Model::LoadDependentRelaxation, Model::LoadDependentRelaxation) =
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
      TimeCoefficients<decltype(std::declval<EvalData>().y)> getTimeCoefficients(
        const EvalData& data) const
      {
        using Scalar = decltype(data.y);
        if (m_input.timeScheme == Model::TimeScheme::BDF2 &&
          data.sn.t > data.snm1.t + Scalar(0.5) * data.dt)
        {
          return {Scalar(1.5) / data.dt, -Scalar(2) / data.dt, Scalar(0.5) / data.dt};
        }
        return {Scalar(1) / data.dt, -Scalar(1) / data.dt, Scalar(0)};
      }

      template <class Scalar>
      static Scalar timeDerivative(Scalar current, Scalar previous,
        Scalar previousPrevious, const TimeCoefficients<Scalar>& tc)
      {
        return tc.current * current + tc.previous * previous +
          tc.previousPrevious * previousPrevious;
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

      // Recruitment fraction n_0(e_c): piecewise linear, identical to CCMLC2014.
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
      ActiveRateData<decltype(std::declval<EvalData>().y)> evaluateActiveRates(
        const EvalData& data, decltype(std::declval<EvalData>().y) ecDot) const
      {
        using Scalar = decltype(data.y);
        ActiveRateData<Scalar> result;
        const Scalar u = m_input.u(data.tnp1);
        result.uPositive = std::max<Scalar>(u, Scalar(0));
        result.uNegative = std::max<Scalar>(-u, Scalar(0));

        const auto recruitment = computeRecruitment(data.ec);
        result.recruitment = recruitment.value;
        result.dRecruitmentDEc = recruitment.derivative;

        result.rate = result.uPositive + data.w * result.uNegative +
          m_input.alpha * regularizedAbs(ecDot);
        result.dRateDEcDot = m_input.alpha * regularizedAbsDerivative(ecDot);
        return result;
      }

      template <class EvalData>
      void evaluateKinematicsAndStresses(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        const auto tc = getTimeCoefficients(data);
        data.a0 = tc.current;
        data.yDot = timeDerivative(data.y, data.sn.y, data.snm1.y, tc);
        data.phiDot = timeDerivative(data.phi, data.sn.phi, data.snm1.phi, tc);
        data.pvMid = data.pv;
        data.parMid = data.par;
        data.pdMid = data.pd;

        m_kinematics.evaluate(data);

        // Series-element fiber stress of the single contractile unit,
        // evaluated at the mid-wall hoop strain.
        const Scalar h = activeStretch(data.ec);
        const Scalar h2 = h * h;
        const Scalar h3 = h2 * h;
        data.active.activeStressOneDimensional =
          m_input.Es / h2 * (data.strain1D - data.ec);
        data.active.partialActiveStressWrtDisplacement = m_input.Es / h2 * data.diffGreen;
        data.active.partialActiveStressWrtPorosity =
          m_input.Es / h2 * data.diffGreenWrtPorosity;
        data.active.partialActiveStressWrtFiberDeformation =
          m_input.Es * (-Scalar(1) / h2 - Scalar(4) * (data.strain1D - data.ec) / h3);

        m_wall.evaluate(data);

        // Interstitial pressure and perfusion source.
        data.pf = m_input.KPhi * (data.phi - m_input.phi0) - data.lambdaBar;
        data.dPf_dy = -data.dLambdaBar_dy;
        data.dPf_dphi = m_input.KPhi - data.dLambdaBar_dphi;
        data.dPf_dec = -data.dLambdaBar_dec;
        data.perfusionInflow = m_input.gammaAr * (data.par - data.pf);
        data.perfusionOutflow = m_input.gammaVen * (data.pf - data.pSvMid);
      }

      template <class EvalData>
      void evaluateActiveDiagnostics(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        const auto tc = getTimeCoefficients(data);
        const Scalar ecDot = timeDerivative(data.ec, data.sn.ec, data.snm1.ec, tc);
        const auto active = evaluateActiveRates(data, ecDot);

        data.active.activationDrive = m_input.u(data.tnp1);
        data.active.activationDrivePositivePart = active.uPositive;
        data.active.activationDriveNegativePart = active.uNegative;
        data.active.relaxationTarget = m_input.m0(data.ec);
        data.active.relaxationDrive = active.rate;
        data.active.recruitmentFraction = active.recruitment;
      }

      template <class EvalData>
      void evaluateValveDiagnostics(EvalData& data) const
      {
        using Scalar = decltype(data.y);

        const bool mitralOpenCurrent = data.pv <= data.pAtCur;
        const bool bothClosedCurrent = data.pAtCur <= data.pv && data.pv <= data.par;

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
            m_input.Kar * (data.pv - data.par) + m_input.Kp * (data.par - data.pAtCur);
          data.dCavityFluxCur_dPv = m_input.Kar;
          data.dCavityFluxCur_dPar = -m_input.Kar + m_input.Kp;
        }

        using namespace CCMLC2014::Physics;
        if (m_input.windkesselRheology == Model::WindkesselRheology::CarreauYasuda)
        {
          WindkesselOutflowEvaluator<Input, Rheology::CarreauYasuda> windkessel(m_input);
          windkessel.evaluate(data);
        }
        else if (m_input.windkesselRheology == Model::WindkesselRheology::Cross)
        {
          WindkesselOutflowEvaluator<Input, Rheology::Cross> windkessel(m_input);
          windkessel.evaluate(data);
        }
        else if (m_input.windkesselRheology == Model::WindkesselRheology::PowerLaw)
        {
          WindkesselOutflowEvaluator<Input, Rheology::PowerLaw> windkessel(m_input);
          windkessel.evaluate(data);
        }
        else if (m_input.windkesselRheology == Model::WindkesselRheology::Quemada)
        {
          WindkesselOutflowEvaluator<Input, Rheology::Quemada> windkessel(m_input);
          windkessel.evaluate(data);
        }
        else
        {
          WindkesselOutflowEvaluator<Input, Rheology::Newtonian> windkessel(m_input);
          windkessel.evaluate(data);
        }
      }

      const Input& m_input;
      Physics::WallKinematics<Input> m_kinematics;
      Physics::WallStressEvaluator<PassiveLaw, Input> m_wall;
  };
}

#endif
