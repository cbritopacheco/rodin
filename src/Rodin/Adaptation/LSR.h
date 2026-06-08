/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_LSR_H
#define RODIN_ADAPTATION_LSR_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Rodin/Assembly.h"
#include "Rodin/Solid/Linear/LinearElasticityIntegral.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Variational.h"

#include "CellGeomCache.h"
#include "JacobianAdmissibilityBarrierSampled.h"
#include "LSRAdmissibility.h"
#include "LSRParameters.h"
#include "LSRReport.h"
#include "LSRRegistration.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Level-set registration solver.
   *
   * `LSR` owns no mesh data. It is bound to the displacement field `u`,
   * infers the carrier finite-element space and mesh from that field, and
   * overwrites `u` with the accepted registration displacement.
   *
   * ## Dimensionless energy formulation
   *
   * All user weights γ_* are dimensionless O(1) numbers. The internal
   * dimensional factor that brings every term to E_push's energy scale
   * is `lengthScaling = normalizer · h_ref²`, where
   *   normalizer = 1 / (M_w · h_ref²)        (auto)
   *   M_w        = ∫_Ω W(ψ) dX               (auto)
   *
   * Total objective (push branch):
   *
   *   E_total(u) =   γ_push   · normalizer · ½ ∫ W(ψ)(φ(X+u) − ψ(X))² dX
   *               + γ_Γ      · interfaceNormalizer · ½ ∫_Γ (φ(X+u)/|∇φ|)² ds
   *               + γ_shape  · normalizer · h_ref² · ∫ B_shape(I + ∇u) dX
   *               + γ_jBar   · normalizer · h_ref² · ∫ B_j(j) dX
   *               + γ_vol    · normalizer · h_ref² · ½ ∫ (log j)² dX
   *               + γ_damp   · normalizer        · ½ ∫ (1−W(ψ))² |u|² dX
   *               + γ_H1     · normalizer · h_ref² · ½ ∫ |∇u|² dX
   *
   * The Hilbert lift used by the initial guess can be stiffened outside
   * the band via the dimensionless `bandHilbertStiffness` s_H, with
   *   a(u, v) = ∫ (1 + s_H (1−W)²) ∇u : ∇v   (Harmonic metric only).
   *
   * ## Robust defaults from the (P_k, n, lobes) sweep
   *
   * Tested on the wavy-circle reconstruction with ψ = projected-phi-h1.
   * The interface term is assembled only by the push branch.
   *
   *   γ_push   = 1
   *   γ_Γ      = min(cap_k, 1 / h_ref)      # cap_1 = 10, cap_2 = 100
   *   γ_shape  = sqrt(h_ref)                # per-FE-order auto-scaling
   *   γ_damp   = 0                          # off by default
   *   γ_jBar   = 1                          # active near j < jBarSafe
   *   jBarSafe = 0.5
   *   γ_vol    = 0.01
   *   γ_H1     = 0
   *   s_H      = 0
   *
   * These values are empirical defaults for the tested benchmark family,
   * not theoretical constants. The optional safety-net terms should be
   * enabled separately, according to the observed failure mode:
   *
   *   γ_damp  > 0   : cell distortion in deep interior
   *                   (Hilbert lift's harmonic extension extremum).
   *   γ_H1    > 0   : extra C¹ smoothness on u (rarely needed).
   *   s_H     > 0   : Hilbert lift backs to α=0 (inadmissible Riesz lift,
   *                   typically from a sharp non-smooth ψ).
   *
   * γ_shape can NOT be set to zero: without the Q_rel barrier Newton
   * loses cell-quality control and stops converging (sweep showed
   * fit ≈ 9e-3 NC at γ_shape = 0).
   */
  template <class Displacement>
  class LSR
  {
    public:
      explicit LSR(Displacement& u)
        : m_u(u)
      {}

      LSR& setParameters(LSRParameters params)
      {
        m_params = std::move(params);
        return *this;
      }

      const LSRParameters& getParameters() const noexcept
      {
        return m_params;
      }

      const LSRReport& getReport() const noexcept
      {
        return m_report;
      }

      template <class Psi, class PhiDerived, class GradDerived, class HessDerived>
      LSRReport solve(
          const Psi& psi,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& gradPhi,
          const Variational::MatrixFunctionBase<Real, HessDerived>& hessPhi)
      {
        using Variational::DirichletBC;
        using Variational::Integral;
        using Variational::Jacobian;
        using Variational::Problem;
        using Variational::RealFunction;
        using Variational::TestFunction;
        using Variational::TrialFunction;
        using Variational::VectorFunction;
        using Variational::Zero;

        m_report = {};
        m_prevAcceptedAlpha = Real(0);

        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        LSRParameters params = m_params;
        completeParameters(params, psi);
        const Real jLineSearchRatio =
          std::max(params.jMinRatio,
                   params.lineSearchSafetyMargin * params.jSafeRatio);
        m_report.jLineSearchRatio = jLineSearchRatio;
        const bool barrierEnabled =
             params.shapeWeight != Real(0)
          || params.qBarrierWeight != Real(0)
          || params.jBarrierWeight != Real(0)
          || params.jVolumeTetherWeight != Real(0);

        TrialFunction du(fes);
        TestFunction  v(fes);
        auto zero = VectorFunction{ Zero(), Zero() };

        // Scale-aware dimensionless weights.
        //
        //   E_push  ~ γ_push  · normalizer · ∫ W(φ−ψ)² dX                 (existing)
        //   E_damp  ~ γ_damp  · normalizer · ∫ (1−W)² |u|² dX             (existing)
        //   E_H1    ~ γ_H1    · normalizer · h_ref² · ∫ |∇u|² dX          (NEW)
        //   E_shape ~ γ_shape · normalizer · h_ref² · ∫ B_shape(F) dX     (NEW)
        //   E_barrier (jacobian/volume) ~ γ_* · normalizer · h_ref²        (NEW)
        //
        // The `normalizer · h_ref²` factor turns per-qpt B integrands
        // with units of [energy/length^d] into dimensionless E_push-scale
        // quantities, so every γ_* is O(1) and h-independent.
        //
        // We mutate the LOCAL `params` copy in place so BOTH the tangent
        // assembly and the line-search objective evaluation see the same
        // scaled weights.
        const Real lengthScaling = params.normalizer * params.hRef * params.hRef;
        params.shapeWeight *= lengthScaling;
        params.h1RegularizationWeight *= lengthScaling;
        params.jBarrierWeight *= lengthScaling;
        params.jVolumeTetherWeight *= lengthScaling;
        params.qBarrierWeight *= lengthScaling;
	        RealFunction<Real> shapeWeight(params.shapeWeight);
	        RealFunction<Real> h1RegularizationWeight(
	            params.h1RegularizationWeight);
        // Outside-band L^2 damping weight.
        //   E_damp = 0.5 * (w_damp * normalizer) * (1 - W(psi))^2 * |u|^2.
        // Multiplying by `normalizer` (= 1 / (M_w h_ref^2)) makes
        // w_damp dimensionless and *scale-aware*: w_damp = 1 puts
        // E_damp on the same scale as E_push when |u| ~ h. Reads psi
        // by reference (captured by the lambda).
        const Real outsideBandWeightValue =
          params.outsideBandTikhonovWeight * params.normalizer;
        const Real deltaWLocal = params.deltaW;
        RealFunction outsideBandWeight(
            [&, outsideBandWeightValue, deltaWLocal](
                const Geometry::Point& p) -> Real
            {
              if (outsideBandWeightValue == Real(0))
                return Real(0);
              const Real s = psi.getValue(p);
              const Real W = std::exp(
                  -s * s / (Real(2) * deltaWLocal * deltaWLocal));
              const Real complement = Real(1) - W;
              return outsideBandWeightValue * complement * complement;
            });

        BarrierParameters barrierParams;
        barrierParams.jMin = params.jMinRatio;
        barrierParams.domainMeasure = computeDomainMeasure();
        barrierParams.qBarrierWeight = params.qBarrierWeight;
        barrierParams.qBarrierAct    = params.qBarrierAct;
        barrierParams.qBarrierMax    = params.qBarrierMax;
        barrierParams.jBarrierWeight = params.jBarrierWeight;
        barrierParams.jBarrierSafeRatio = params.jBarrierSafeRatio;
        barrierParams.jVolumeTetherWeight = params.jVolumeTetherWeight;

        LSRRegistration lsrTerm(du, v, u);
        JacobianAdmissibilityBarrierSampled barrier(
            du, v, u, params.quadratureOrder);

        if (params.initialGuess == LSRInitialGuess::Zero)
          u.getData().setZero();
        else if (params.initialGuess == LSRInitialGuess::Hilbert)
          applyHilbertInitialGuess(
              psi, phi, gradPhi, shapeWeight,
              barrierParams, params);

        // Adaptive σ via MAD on active-band residuals at the initial
        // guess. Computed once per solve(); overrides params.lossSigma.
        if (params.useAdaptiveLossSigma)
        {
          const Real sigmaEstimate =
            estimateLossSigmaMAD(psi, phi, params);
          if (sigmaEstimate > Real(0))
            params.lossSigma = sigmaEstimate;
        }

        Problem newton(du, v);
        switch (params.tangent)
        {
          case LSRTangent::GaussNewton:
            if (barrierEnabled)
	              newton =
	                  lsrTerm.Tangent(phi, gradPhi, psi, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, psi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceTangent(phi, gradPhi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceResidual(phi, gradPhi, makeLSRIntegratorParameters(params))
		                + barrier.TangentPSDProjected(shapeWeight, barrierParams)
	                + barrier.Residual(shapeWeight, barrierParams)
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
	                + DirichletBC(du, zero);
	            else
	              newton =
	                  lsrTerm.Tangent(phi, gradPhi, psi, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, psi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceTangent(phi, gradPhi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceResidual(phi, gradPhi, makeLSRIntegratorParameters(params))
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
	                + DirichletBC(du, zero);
            break;
          case LSRTangent::Newton:
            if (barrierEnabled)
	              newton =
	                  lsrTerm.Tangent(phi, gradPhi, hessPhi, psi, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, psi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceTangent(phi, gradPhi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceResidual(phi, gradPhi, makeLSRIntegratorParameters(params))
		                + barrier.TangentPSDProjected(shapeWeight, barrierParams)
	                + barrier.Residual(shapeWeight, barrierParams)
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
	                + DirichletBC(du, zero);
	            else
	              newton =
	                  lsrTerm.Tangent(phi, gradPhi, hessPhi, psi, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, psi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceTangent(phi, gradPhi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceResidual(phi, gradPhi, makeLSRIntegratorParameters(params))
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
	                + DirichletBC(du, zero);
            break;
          case LSRTangent::PSDProjectedNewton:
            if (barrierEnabled)
	              newton =
	                  lsrTerm.TangentPSDProjected(
	                      phi, gradPhi, hessPhi, psi, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, psi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceTangent(phi, gradPhi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceResidual(phi, gradPhi, makeLSRIntegratorParameters(params))
	                + barrier.TangentPSDProjected(shapeWeight, barrierParams)
	                + barrier.Residual(shapeWeight, barrierParams)
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
                  + Integral(outsideBandWeight * du, v)
                  + Integral(outsideBandWeight * u, v)
	                + DirichletBC(du, zero);
	            else
	              newton =
	                  lsrTerm.TangentPSDProjected(
	                      phi, gradPhi, hessPhi, psi, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, psi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceTangent(phi, gradPhi, makeLSRIntegratorParameters(params))
                  + lsrTerm.InterfaceResidual(phi, gradPhi, makeLSRIntegratorParameters(params))
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
                  + Integral(outsideBandWeight * du, v)
                  + Integral(outsideBandWeight * u, v)
	                + DirichletBC(du, zero);
            break;
        }

        auto objective =
          [&](const Math::Vector<Real>& uTry) -> Real
          {
            u.getData() = uTry;
            return computeObjective(
                psi, phi, gradPhi, params, barrierParams, barrierEnabled);
          };

        return driveNewton(
            newton, params, objective,
            barrierParams, jLineSearchRatio);
      }

      /**
       * @brief Pull-back LSR solve.
       *
       * E_pull(u) = (rhoS/2) integral W(psi(X)) ( phi(X) - psi(X - u(X)) )^2 dX.
       *
       * The data residual evaluates phi at the original X and psi at the
       * pulled-back point X - u(X). The caller supplies:
       *
       *   psi          : the band-weight psi field (used in W(psi(X))).
       *   phi, gradPhi : at original X. gradPhi is only needed for the
       *                  Hilbert initial guess; the data term itself uses
       *                  grad psi at the pulled-back point.
       *   psiDisp      : RealFunctionBase returning psi at X - u(X).
       *   gradPsiDisp  : VectorFunctionBase returning grad psi at X - u(X).
       *   hessPsiDisp  : MatrixFunctionBase returning Hess psi at X - u(X).
       *
       * Tangent. The same LSRTangent selector used by the push form is
       * honoured here. The exact pull Newton correction is
       * `-r * hess(psi)(X-u)`, with `r = phi(X) - psi(X-u)`.
       */
      template <class Psi, class PhiDerived, class GradPhiDerived,
                class PsiDispDerived, class GradPsiDispDerived,
                class HessPsiDispDerived>
      LSRReport solvePull(
          const Psi& psi,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradPhiDerived>& gradPhi,
          const Variational::RealFunctionBase<PsiDispDerived>& psiDisp,
          const Variational::VectorFunctionBase<Real, GradPsiDispDerived>& gradPsiDisp,
          const Variational::MatrixFunctionBase<Real, HessPsiDispDerived>& hessPsiDisp)
      {
        using Variational::DirichletBC;
        using Variational::Integral;
        using Variational::Jacobian;
        using Variational::Problem;
        using Variational::RealFunction;
        using Variational::TestFunction;
        using Variational::TrialFunction;
        using Variational::VectorFunction;
        using Variational::Zero;

        m_report = {};
        m_prevAcceptedAlpha = Real(0);

        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();

        LSRParameters params = m_params;
        completeMeshParameters(params);
        const Real jLineSearchRatio =
          std::max(params.jMinRatio,
                   params.lineSearchSafetyMargin * params.jSafeRatio);
        m_report.jLineSearchRatio = jLineSearchRatio;
        const bool barrierEnabled =
             params.shapeWeight != Real(0)
          || params.qBarrierWeight != Real(0);

        TrialFunction du(fes);
        TestFunction  v(fes);
        auto zero = VectorFunction{ Zero(), Zero() };

        RealFunction<Real> shapeWeight(params.shapeWeight);
        RealFunction<Real> h1RegularizationWeight(
            params.h1RegularizationWeight);

        BarrierParameters barrierParams;
        barrierParams.jMin = params.jMinRatio;
        barrierParams.domainMeasure = computeDomainMeasure();
        barrierParams.qBarrierWeight = params.qBarrierWeight;
        barrierParams.qBarrierAct    = params.qBarrierAct;
        barrierParams.qBarrierMax    = params.qBarrierMax;
        barrierParams.jBarrierWeight = params.jBarrierWeight;
        barrierParams.jBarrierSafeRatio = params.jBarrierSafeRatio;
        barrierParams.jVolumeTetherWeight = params.jVolumeTetherWeight;

        LSRRegistration lsrTerm(du, v, u);
        JacobianAdmissibilityBarrierSampled barrier(
            du, v, u, params.quadratureOrder);

        if (params.initialGuess == LSRInitialGuess::Zero)
          u.getData().setZero();
        else if (params.initialGuess == LSRInitialGuess::Hilbert)
          applyHilbertInitialGuess(
              psi, phi, gradPhi, shapeWeight,
              barrierParams, params);

        Problem newton(du, v);
        switch (params.tangent)
        {
          case LSRTangent::GaussNewton:
            if (barrierEnabled)
              newton =
                  lsrTerm.PullTangent(psi, gradPsiDisp, makeLSRIntegratorParameters(params))
                + lsrTerm.PullResidual(phi, psi, psiDisp, gradPsiDisp,
                                       makeLSRIntegratorParameters(params))
                + barrier.TangentPSDProjected(shapeWeight, barrierParams)
                + barrier.Residual(shapeWeight, barrierParams)
                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
                + DirichletBC(du, zero);
            else
              newton =
                  lsrTerm.PullTangent(psi, gradPsiDisp, makeLSRIntegratorParameters(params))
                + lsrTerm.PullResidual(phi, psi, psiDisp, gradPsiDisp,
                                       makeLSRIntegratorParameters(params))
                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
                + DirichletBC(du, zero);
            break;
          case LSRTangent::Newton:
            if (barrierEnabled)
              newton =
                  lsrTerm.PullTangent(
                      phi, psi, psiDisp, gradPsiDisp, hessPsiDisp,
                      makeLSRIntegratorParameters(params))
                + lsrTerm.PullResidual(phi, psi, psiDisp, gradPsiDisp,
                                       makeLSRIntegratorParameters(params))
                + barrier.TangentPSDProjected(shapeWeight, barrierParams)
                + barrier.Residual(shapeWeight, barrierParams)
                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
                + DirichletBC(du, zero);
            else
              newton =
                  lsrTerm.PullTangent(
                      phi, psi, psiDisp, gradPsiDisp, hessPsiDisp,
                      makeLSRIntegratorParameters(params))
                + lsrTerm.PullResidual(phi, psi, psiDisp, gradPsiDisp,
                                       makeLSRIntegratorParameters(params))
                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
                + DirichletBC(du, zero);
            break;
          case LSRTangent::PSDProjectedNewton:
            if (barrierEnabled)
              newton =
                  lsrTerm.PullTangentPSDProjected(
                      phi, psi, psiDisp, gradPsiDisp, hessPsiDisp,
                      makeLSRIntegratorParameters(params))
                + lsrTerm.PullResidual(phi, psi, psiDisp, gradPsiDisp,
                                       makeLSRIntegratorParameters(params))
                + barrier.TangentPSDProjected(shapeWeight, barrierParams)
                + barrier.Residual(shapeWeight, barrierParams)
                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
                + DirichletBC(du, zero);
            else
              newton =
                  lsrTerm.PullTangentPSDProjected(
                      phi, psi, psiDisp, gradPsiDisp, hessPsiDisp,
                      makeLSRIntegratorParameters(params))
                + lsrTerm.PullResidual(phi, psi, psiDisp, gradPsiDisp,
                                       makeLSRIntegratorParameters(params))
                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
                + DirichletBC(du, zero);
            break;
        }

        auto objective =
          [&](const Math::Vector<Real>& uTry) -> Real
          {
            u.getData() = uTry;
            return computeObjectivePull(
                psi, phi, psiDisp, params, barrierParams, barrierEnabled);
          };

        return driveNewton(
            newton, params, objective,
            barrierParams, jLineSearchRatio);
      }

    private:
      // Newton driver for the standard E_LSR solve.
      // Takes the already-assembled Problem; sets up SparseLU + NewtonSolver
      // with the same line-search + admissibility machinery either way.
      template <class Problem_, class Objective>
      LSRReport driveNewton(
          Problem_& newton,
          LSRParameters& params,
          Objective&& objective,
          const BarrierParameters& barrierParams,
          Real jLineSearchRatio)
      {
        auto& u = m_u.get();

        Solver::SparseLU linearSolver(newton);
        Solver::NewtonSolver newtonSolver(linearSolver);

        auto evaluator = [&](const Math::Vector<Real>& uTry)
        {
          return evaluateLSRAdmissibilitySampled(
              u, uTry, params.jMinRatio, params.quadratureOrder);
        };

        auto residualNorm = [&](const Math::Vector<Real>& uTry) -> Real
        {
          u.getData() = uTry;
          newton.assemble();
          return newton.getLinearSystem().getVector().norm();
        };

        Math::Vector<Real> uBest = u.getData();
        Real energyBest = std::numeric_limits<Real>::infinity();
        Real residualBest = std::numeric_limits<Real>::infinity();
        std::size_t stallCount = 0;
        bool initialized = false;

        // --- Trust region state (A3 + A4) ---
        // trustRadius caps ‖du‖_∞ per Newton iteration. It grows on
        // successful steps with little backtracking and shrinks on
        // step rejections. Set trustRadiusFactorInit ≤ 0 to disable.
        const Real h_ref = std::max(params.hRef, Real(1e-30));
        const bool trustRegionEnabled = params.trustRadiusFactorInit > Real(0);
        Real trustRadius = trustRegionEnabled
          ? params.trustRadiusFactorInit * h_ref
          : std::numeric_limits<Real>::infinity();
        const Real trustRadiusMin =
          std::max(params.trustRadiusFactorMin * h_ref, params.alphaMin * h_ref);
        const Real trustRadiusMax = params.trustRadiusFactorMax * h_ref;
        std::size_t trustRadiusShrinks = 0;
        std::size_t trustRadiusExpansions = 0;

        newtonSolver
          .setMaxIterations(params.maxNewtonIterations)
          .setAbsoluteTolerance(params.absoluteTolerance)
          .setRelativeTolerance(params.relativeTolerance)
          .setStepTolerance(Real(0))
          .setDampingFactor(Real(1));

        newtonSolver.setStepPolicy(
          typename decltype(newtonSolver)::StepPolicy(
          [&](auto& x, auto& linearSystem, auto& solverReport)
            -> typename decltype(newtonSolver)::StepResult
        {
          const Math::Vector<Real> uOld = x;
          const Real residual = solverReport.final_residual;
          const Real energy = objective(uOld);

          if (!std::isfinite(residual) || !std::isfinite(energy))
          {
            m_report.converged = false;
            x = uBest;
            return {false, false, Real(0)};
          }

          if (!initialized)
          {
            m_report.initialResidual = residual;
            m_report.initialEnergy = energy;
            uBest = uOld;
            energyBest = energy;
            residualBest = residual;
            initialized = true;
          }

          Math::Vector<Real> newtonDirection = linearSystem.getSolution();
          Math::Vector<Real> residualDirection = linearSystem.getVector();
          const Real newtonDirectionNorm = newtonDirection.norm();
          if (!std::isfinite(newtonDirectionNorm))
          {
            m_report.converged = false;
            x = uBest;
            return {false, false, Real(0)};
          }

          // --- A4: step-norm cap via trust radius ---
          // Bound ‖du‖_∞ ≤ trustRadius. Direction is scaled by a
          // single scalar; sign and shape preserved.
          bool stepWasCapped = false;
          if (trustRegionEnabled && std::isfinite(trustRadius))
          {
            const Real maxCoeff =
              newtonDirection.cwiseAbs().maxCoeff();
            if (maxCoeff > trustRadius && maxCoeff > Real(0))
            {
              const Real scaleFactor = trustRadius / maxCoeff;
              newtonDirection *= scaleFactor;
              stepWasCapped = true;
            }
          }

          LSRLineSearchResult ls;
          Real acceptedResidual = std::numeric_limits<Real>::quiet_NaN();
          Real acceptedEnergy = std::numeric_limits<Real>::quiet_NaN();
          Math::Vector<Real> acceptedU = uOld;
          Real acceptedStepNorm = Real(0);
          LSRAdmissibilityReport lastAdm;
          const Real energyTol =
            params.energyDecreaseTolerance
            * std::max<Real>(Real(1), std::abs(energy));

          auto tryLineSearchDirection =
            [&](const Math::Vector<Real>& direction) -> bool
          {
            const Real directionNorm = direction.norm();
            if (!std::isfinite(directionNorm) || directionNorm == Real(0))
              return false;

            Real alphaStart = params.alphaInit;
            if (params.useWarmStartAlpha
                && m_prevAcceptedAlpha > Real(0))
            {
              alphaStart = std::min<Real>(
                  Real(1),
                  params.alphaWarmStartGrowth * m_prevAcceptedAlpha);
              alphaStart = std::max<Real>(alphaStart, params.alphaMin);
            }
            if (params.useSampledQuadraticAlphaPredictor)
            {
              const auto predicted =
                predictSampledQuadraticAlpha(
                    u, uOld, direction, jLineSearchRatio,
                    params.alphaPredictorSafety,
                    params.quadratureOrder);
              if (std::isfinite(predicted.alphaMax))
                alphaStart = std::min(alphaStart, predicted.alphaMax);
            }

            Real alpha = alphaStart;
            while (alpha >= params.alphaMin)
            {
              const Math::Vector<Real> uTrial = uOld + alpha * direction;
              const LSRAdmissibilityReport adm = evaluator(uTrial);
              lastAdm = adm;

              const bool jOK =
                   adm.minJRatio > jLineSearchRatio
                && adm.inadmissibleCount == 0;
              const bool qRelOK =
                std::isfinite(params.qRelMax)
                  ? adm.maxQRel <= params.qRelMax
                  : true;
              bool energyOK = false;
              Real trialEnergy = std::numeric_limits<Real>::infinity();
              if (jOK && qRelOK)
              {
                trialEnergy = objective(uTrial);
                energyOK =
                     std::isfinite(trialEnergy)
                  && trialEnergy <= energy + energyTol;
              }

              if (jOK && qRelOK && energyOK)
              {
                acceptedResidual = residualNorm(uTrial);
                if (std::isfinite(acceptedResidual))
                {
                  acceptedEnergy = trialEnergy;
                  acceptedU = uTrial;
                  acceptedStepNorm = alpha * directionNorm;
                  ls.alphaAccepted = alpha;
                  ls.minJRatioAccepted = adm.minJRatio;
                  ls.inadmissibleCountAccepted = adm.inadmissibleCount;
                  ls.maxQRelAccepted = adm.maxQRel;
                  ls.succeeded = true;
                  return true;
                }
              }

              alpha *= params.alphaReduction;
              ++ls.backtracks;
            }
            return false;
          };

          const bool acceptedNewton =
            tryLineSearchDirection(newtonDirection);
          if (!acceptedNewton)
          {
            const bool acceptedResidualDirection =
              tryLineSearchDirection(residualDirection);
            if (!acceptedResidualDirection)
              tryLineSearchDirection(-residualDirection);
          }

          if (!ls.succeeded)
          {
            x = uOld;
            ls.minJRatioAccepted = lastAdm.minJRatio;
            ls.inadmissibleCountAccepted = lastAdm.inadmissibleCount;
            ls.maxQRelAccepted = lastAdm.maxQRel;
          }
          else
          {
            x = acceptedU;
          }

          m_report.iterations = solverReport.iterations + 1;
          m_report.finalResidual =
            ls.succeeded ? acceptedResidual : residualBest;
          m_report.finalEnergy =
            ls.succeeded ? acceptedEnergy : energyBest;
          m_report.finalStepNorm = acceptedStepNorm;
          m_report.lastAcceptedAlpha = ls.alphaAccepted;
          m_report.totalBacktracks += ls.backtracks;
          if (ls.succeeded)
            m_prevAcceptedAlpha = ls.alphaAccepted;
          m_report.minJRatio = ls.minJRatioAccepted;
          m_report.inadmissibleCount = ls.inadmissibleCountAccepted;

          // --- A3: trust-region radius update ---
          // Heuristic based on accepted α as a proxy for model quality.
          //   α ≥ growExpandThreshold (≈ 0.5): step matched the model
          //      well — expand Δ.
          //   α ≤ shrinkThreshold (≈ 0.1): step needed heavy line
          //      search — shrink Δ.
          //   Line search failed: shrink Δ aggressively.
          // Δ only adapts when the cap was active or the step needed
          // significant backtracking; healthy uncapped Newton steps
          // with α = 1 leave Δ alone. Otherwise default Δ_init can be
          // smaller than the natural step and Δ collapses on iter 0.
          if (trustRegionEnabled)
          {
            if (ls.succeeded)
            {
              if (stepWasCapped
                  && ls.alphaAccepted >= params.trustRadiusGrowExpandThresh
                  && ls.backtracks == 0)
              {
                trustRadius = std::min(
                  trustRadius * params.trustRadiusGrowRate,
                  trustRadiusMax);
                ++trustRadiusExpansions;
              }
              else if (ls.alphaAccepted <= params.trustRadiusShrinkThresh
                       || ls.backtracks >= 4)
              {
                trustRadius = std::max(
                  trustRadius * params.trustRadiusShrinkRate,
                  trustRadiusMin);
                ++trustRadiusShrinks;
              }
            }
            else
            {
              trustRadius = std::max(
                trustRadius * params.trustRadiusFailShrinkRate,
                trustRadiusMin);
              ++trustRadiusShrinks;
            }
            m_report.lastTrustRadius = trustRadius;
            m_report.trustRadiusShrinks = trustRadiusShrinks;
            m_report.trustRadiusExpansions = trustRadiusExpansions;
          }

          if (!ls.succeeded)
          {
            m_report.lineSearchFailed = true;
            m_report.converged = false;
            x = uBest;
            solverReport.final_residual = m_report.finalResidual;
            solverReport.final_step_norm = m_report.finalStepNorm;
            return {false, false, m_report.finalStepNorm};
          }

          const bool energyImproved =
            acceptedEnergy < energyBest
              - params.energyDecreaseTolerance
                * std::max<Real>(Real(1), std::abs(energyBest));
          if (energyImproved)
          {
            uBest = u.getData();
            energyBest = acceptedEnergy;
            residualBest = acceptedResidual;
            stallCount = 0;
          }
          else
          {
            ++stallCount;
          }

          if (params.acceptedStateConvergenceTest
              && params.acceptedStateConvergenceTest(m_report))
          {
            m_report.converged = true;
            solverReport.final_residual = m_report.finalResidual;
            solverReport.final_step_norm = m_report.finalStepNorm;
            return {true, true, m_report.finalStepNorm};
          }

          if (acceptedResidual <= params.absoluteTolerance
              || (m_report.initialResidual > 0
                  && acceptedResidual <= params.relativeTolerance
                    * m_report.initialResidual))
          {
            m_report.converged = true;
            solverReport.final_residual = m_report.finalResidual;
            solverReport.final_step_norm = m_report.finalStepNorm;
            return {true, true, m_report.finalStepNorm};
          }

          if (params.stepTolerance > Real(0)
              && m_report.finalStepNorm <= params.stepTolerance)
          {
            m_report.converged = true;
            solverReport.final_residual = m_report.finalResidual;
            solverReport.final_step_norm = m_report.finalStepNorm;
            return {true, true, m_report.finalStepNorm};
          }

          if (params.stallPatience > 0 && stallCount >= params.stallPatience)
          {
            m_report.converged = true;
            x = uBest;
            m_report.finalResidual = residualBest;
            m_report.finalEnergy = energyBest;
            solverReport.final_residual = m_report.finalResidual;
            solverReport.final_step_norm = m_report.finalStepNorm;
            return {true, true, m_report.finalStepNorm};
          }

          solverReport.final_residual = m_report.finalResidual;
          solverReport.final_step_norm = m_report.finalStepNorm;
          return {true, false, m_report.finalStepNorm};
        }));

        newtonSolver.solve(u.getData());
        const auto& solverReport = newtonSolver.getReport();

        if (!initialized)
        {
          m_report.iterations = solverReport.iterations;
          m_report.initialResidual = solverReport.initial_residual;
          m_report.finalResidual = solverReport.final_residual;
          m_report.converged = solverReport.converged;
          return m_report;
        }

        if (!solverReport.converged && !m_report.lineSearchFailed)
        {
          u.getData() = uBest;
          m_report.finalResidual = residualBest;
          m_report.finalEnergy = energyBest;
        }
        m_report.converged = solverReport.converged;
        return m_report;
      }

      void completeMeshParameters(LSRParameters& params) const
      {
        const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
        const Real domainMeasure = computeDomainMeasure();

        if (params.hRef <= 0)
          params.hRef =
            std::sqrt(domainMeasure / static_cast<Real>(mesh.getCellCount()));
        if (params.deltaW <= 0)
          params.deltaW = params.hRef;
        if (params.normalizer <= 0)
          params.normalizer = Real(1);
      }

      /// Estimate Welsch σ via MAD on active-band residuals at the
      /// current `u`. Returns 0 on empty band (no estimate possible).
      ///
      ///   σ ≈ multiplier · 1.4826 · median |r|,
      ///
      /// where the median is over quadrature points with W(ψ) > band-
      /// Threshold. For an approximately symmetric inlier residual
      /// distribution `median(r) ≈ 0`, so MAD = median(|r − median(r)|)
      /// reduces to median(|r|). The factor 1.4826 makes σ the
      /// consistent estimator of the inlier standard deviation; the
      /// multiplier scales the saturation cutoff (residuals beyond
      /// multiplier·σ_inlier get downweighted).
      template <class Psi, class PhiDerived>
      Real estimateLossSigmaMAD(
          const Psi& psi,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const LSRParameters& params) const
      {
        const auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto lsrParams = makeLSRIntegratorParameters(params);

        std::vector<Real> absResiduals;
        absResiduals.reserve(mesh.getCellCount() * 4);

        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        {
          const auto& cell = *cellIt;
          const auto& fe = fes.getFiniteElement(
              cell.getDimension(), cell.getIndex());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(
                lsrQuadOrderFor(fe.getOrder(), lsrParams),
                cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& pt = quadrature.getPoint(q);
            const Variational::IntegrationPoint ip(pt, &qf, q);
            const Real s = psi.getValue(ip);
            const Real W =
              std::exp(-s * s / (2 * params.deltaW * params.deltaW));
            if (W < params.adaptiveLossSigmaBandThreshold)
              continue;
            const auto uq = u.getValue(ip);
            Math::SpatialVector<Real> displacement(cell.getDimension());
            displacement.setZero();
            for (std::size_t c = 0; c < fes.getVectorDimension(); ++c)
              displacement(c) = uq(c);
            const auto traced =
              detail::traceMovedPoint(
                  pt, displacement, phi, lsrParams.fieldEvaluation);
            if (traced.exited)
              continue;
            const Real r =
              phi.getValue(traced.point) + traced.correction - s;
            absResiduals.push_back(std::abs(r));
          }
        }

        if (absResiduals.empty())
          return Real(0);

        const std::size_t mid = absResiduals.size() / 2;
        std::nth_element(
            absResiduals.begin(),
            absResiduals.begin() + mid,
            absResiduals.end());
        const Real medianAbsR = absResiduals[mid];

        return params.adaptiveLossSigmaMultiplier
             * Real(1.4826)
             * medianAbsR;
      }

      template <class Psi, class PhiDerived, class GradDerived>
      Real computeObjective(
          const Psi& psi,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& gradPhi,
          const LSRParameters& params,
          const BarrierParameters& barrierParams,
          bool barrierEnabled) const
      {
        const auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto lsrParams = makeLSRIntegratorParameters(params);

        Real energy = 0;
        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        {
          const auto& cell = *cellIt;
          const auto& fe = fes.getFiniteElement(
              cell.getDimension(), cell.getIndex());
	          const auto& qf =
	            QF::PolytopeQuadratureFormula::get(
	                lsrQuadOrderFor(fe.getOrder(), lsrParams),
	                cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& pt = quadrature.getPoint(q);
            const Variational::IntegrationPoint ip(pt, &qf, q);
            const Real s = psi.getValue(ip);
            const auto uq = u.getValue(ip);
            Math::SpatialVector<Real> displacement(cell.getDimension());
            displacement.setZero();
            for (std::size_t c = 0; c < fes.getVectorDimension(); ++c)
              displacement(c) = uq(c);

            const auto traced =
              detail::traceMovedPoint(
                  pt, displacement, phi, lsrParams.fieldEvaluation);
            if (traced.exited)
              continue;
            const Real r =
              phi.getValue(traced.point) + traced.correction - s;
            const Real W =
              std::exp(-s * s / (2 * params.deltaW * params.deltaW));
            Real lossContrib;
            if (params.lossSigma > Real(0) && params.useWelschEnergy)
            {
              // Family A: Welsch energy ρ(r) = (σ²/2)(1 − exp(−r²/σ²)).
              const Real sigma2 = params.lossSigma * params.lossSigma;
              lossContrib = sigma2 * (Real(1) - std::exp(-r * r / sigma2));
            }
            else
            {
              // Family B (or quadratic baseline): line search sees ½r²
              // — the IRLS-W tangent has already attenuated, but the
              // energy reported here is the raw quadratic loss.
              lossContrib = r * r;
            }
            energy += Real(0.5)
              * qf.getWeight(q) * pt.getDistortion()
              * params.rhoS * W * params.normalizer * lossContrib;
          }
        }

        if (barrierEnabled)
        {
          for (auto cellIt2 = mesh.getCell(); cellIt2; ++cellIt2)
          {
            const auto& barrierCell = *cellIt2;
            const Real e = computeBarrierSampledCellEnergy(
                barrierCell, u, params.shapeWeight, barrierParams);
            if (!std::isfinite(e))
              return std::numeric_limits<Real>::infinity();
            energy += e;
          }
        }
	        if (params.h1RegularizationWeight != Real(0))
	          energy += computeH1RegularizationEnergy(params);
        if (params.interfaceWeight != Real(0))
          energy += computeInterfaceDistanceEnergy(phi, gradPhi, params);
        if (params.outsideBandTikhonovWeight != Real(0))
          energy += computeOutsideBandTikhonovEnergy(psi, params);
	        return energy;
	      }

      template <class PhiDerived, class GradDerived>
      Real computeInterfaceDistanceEnergy(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& gradPhi,
          const LSRParameters& params) const
      {
        const auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto lsrParams = makeLSRIntegratorParameters(params);

        Real energy = 0;
        for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
        {
          const auto& face = *faceIt;
          const auto attr = face.getAttribute();
          if (!attr || *attr != params.interfaceAttribute)
            continue;

          const auto& fe = fes.getFiniteElement(
              face.getDimension(), face.getIndex());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(
                lsrQuadOrderFor(fe.getOrder(), lsrParams),
                face.getGeometry());
          const auto& quadrature = face.getQuadrature(qf);

          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& pt = quadrature.getPoint(q);
            const auto cellIdx = detail::firstIncidentCell(face);
            if (!cellIdx)
              continue;
            const auto cellIt = mesh.getCell(*cellIdx);
            const auto& cell = *cellIt;
            Math::SpatialPoint rcCell(face.getDimension() + 1);
            cell.getTransformation().inverse(
                rcCell, pt.getPhysicalCoordinates());
            const Geometry::Point original(
                cell, rcCell, pt.getPhysicalCoordinates());

            const auto uq = u.getValue(original);
            Math::SpatialVector<Real> displacement(mesh.getSpaceDimension());
            displacement.setZero();
            for (std::size_t c = 0; c < fes.getVectorDimension(); ++c)
              displacement(c) = uq(c);

            const auto traced =
              detail::traceMovedPoint(
                  original, displacement, phi, lsrParams.fieldEvaluation);
            if (traced.exited)
              continue;

            const Real phi_y =
              phi.getValue(traced.point) + traced.correction;
            const auto grad = gradPhi.getValue(traced.point);
            const Real gradNorm =
              std::max(grad.norm(), params.interfaceGradientFloor);
            const Real distance = phi_y / gradNorm;
            energy += Real(0.5)
              * qf.getWeight(q) * pt.getDistortion()
              * params.interfaceWeight
              * params.interfaceNormalizer
              * distance * distance;
          }
        }
        return energy;
      }

      // E_pull objective: (1/2) integral W(psi(X)) ( phi(X) - psi(X-u))^2 dX
      // + barrier + h1 reg.
      template <class Psi, class PhiDerived, class PsiDispDerived>
      Real computeObjectivePull(
          const Psi& psi,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::RealFunctionBase<PsiDispDerived>& psiDisp,
          const LSRParameters& params,
          const BarrierParameters& barrierParams,
          bool barrierEnabled) const
      {
        const auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto lsrParams = makeLSRIntegratorParameters(params);

        Real energy = 0;
        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        {
          const auto& cell = *cellIt;
          const auto& fe = fes.getFiniteElement(
              cell.getDimension(), cell.getIndex());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(
                lsrQuadOrderFor(fe.getOrder(), lsrParams),
                cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& pt = quadrature.getPoint(q);
            const Variational::IntegrationPoint ip(pt, &qf, q);
            const Real s = psi.getValue(ip);
            const Real phi_X = phi.getValue(ip);

            const auto uq = u.getValue(ip);
            const std::size_t d = cell.getDimension();
            const std::size_t vdim = fes.getVectorDimension();
            Math::SpatialVector<Real> yMinus(d);
            for (std::size_t c = 0; c < vdim; ++c)
              yMinus(c) = pt.getCoordinates()(c) - uq(c);
            const Geometry::Point pulledPoint(
                pt.getPolytope(),
                pt.getReferenceCoordinates(),
                yMinus);
            const Real psi_disp = psiDisp.getValue(pulledPoint);

            const Real r = phi_X - psi_disp;
            const Real W =
              std::exp(-s * s / (2 * params.deltaW * params.deltaW));
            energy += Real(0.5)
              * qf.getWeight(q) * pt.getDistortion()
              * params.rhoS * W * params.normalizer * r * r;
          }
        }

        if (barrierEnabled)
        {
          for (auto cellIt2 = mesh.getCell(); cellIt2; ++cellIt2)
          {
            const auto& barrierCell = *cellIt2;
            const Real e = computeBarrierSampledCellEnergy(
                barrierCell, u, params.shapeWeight, barrierParams);
            if (!std::isfinite(e))
              return std::numeric_limits<Real>::infinity();
            energy += e;
          }
        }
        if (params.h1RegularizationWeight != Real(0))
          energy += computeH1RegularizationEnergy(params);
        return energy;
      }

		      Real computeH1RegularizationEnergy(const LSRParameters& params) const
	      {
	        using Variational::IntegrationPoint;
	        using Variational::Jacobian;

	        const auto& u = m_u.get();
	        const auto& fes = u.getFiniteElementSpace();
	        const auto& mesh = fes.getMesh();
	        auto gradU = Jacobian(u);

	        Real energy = 0;
	        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
	        {
	          const auto& cell = *cellIt;
	          const auto& fe =
	            fes.getFiniteElement(cell.getDimension(), cell.getIndex());
	          const auto& qf =
	            QF::PolytopeQuadratureFormula::get(
	                lsrQuadOrderFor(
	                    fe.getOrder(), makeLSRIntegratorParameters(params)),
	                cell.getGeometry());
	          const auto& quadrature = cell.getQuadrature(qf);
	          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
	          {
	            const auto& pt = quadrature.getPoint(q);
	            const IntegrationPoint ip(pt, &qf, q);
	            energy += Real(0.5)
	              * params.h1RegularizationWeight
	              * qf.getWeight(q)
	              * pt.getDistortion()
	              * gradU.getValue(ip).squaredNorm();
	          }
	        }
	        return energy;
	      }

      // E_damp = 0.5 * w_damp * int_Omega (1 - W(psi))^2 |u|^2 dX.
      template <class Psi>
      Real computeOutsideBandTikhonovEnergy(
          const Psi& psi, const LSRParameters& params) const
      {
        const auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        Real energy = 0;
        const Real deltaW = params.deltaW;
        const Real w_damp = params.outsideBandTikhonovWeight * params.normalizer;
        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        {
          const auto& cell = *cellIt;
          const auto& fe =
            fes.getFiniteElement(cell.getDimension(), cell.getIndex());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(
                lsrQuadOrderFor(
                    fe.getOrder(), makeLSRIntegratorParameters(params)),
                cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& pt = quadrature.getPoint(q);
            const Real s = psi.getValue(pt);
            const Real W = std::exp(-s * s / (Real(2) * deltaW * deltaW));
            const Real complement = Real(1) - W;
            const auto uq = u.getValue(pt);
            Real uSq = 0;
            for (int c = 0; c < uq.size(); ++c)
              uSq += uq(c) * uq(c);
            energy += Real(0.5)
              * w_damp * complement * complement
              * qf.getWeight(q)
              * pt.getDistortion()
              * uSq;
          }
        }
        return energy;
      }


      template <class Psi>
      Real computeWeightedBandMeasure(const Psi& psi, Real deltaW) const
      {
        const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
        Real weightedBandMeasure = 0;

        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        {
          const auto& cell = *cellIt;
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(cell.getDimension(), cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& pt = quadrature.getPoint(q);
            Math::SpatialMatrix<Real> J;
            cell.getTransformation().jacobian(J, pt.getReferenceCoordinates());
            const Real detJ = std::abs(J.determinant());
            const Real wq = qf.getWeight(q) * detJ;
            const Real s = psi.getValue(pt);
            const Real W = std::exp(-s * s / (2 * deltaW * deltaW));
            weightedBandMeasure += W * wq;
          }
        }
        return weightedBandMeasure;
      }

      Real computeDomainMeasure() const
      {
        const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
        Real measure = 0;
        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        {
          const auto& cell = *cellIt;
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(cell.getDimension(), cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& pt = quadrature.getPoint(q);
            Math::SpatialMatrix<Real> J;
            cell.getTransformation().jacobian(J, pt.getReferenceCoordinates());
            measure += qf.getWeight(q) * std::abs(J.determinant());
          }
        }
        return measure;
      }

      Real computeInterfaceMeasure(Geometry::Attribute attr) const
      {
        const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
        Real measure = 0;
        for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
        {
          const auto& face = *faceIt;
          const auto faceAttr = face.getAttribute();
          if (!faceAttr || *faceAttr != attr)
            continue;
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(
                face.getDimension(), face.getGeometry());
          const auto& quadrature = face.getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
            measure += qf.getWeight(q) * quadrature.getPoint(q).getDistortion();
        }
        return measure;
      }

      template <class Psi>
      void completeParameters(LSRParameters& params, const Psi& psi) const
      {
        const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
        const Real domainMeasure = computeDomainMeasure();

        if (params.hRef <= 0)
          params.hRef = std::sqrt(domainMeasure / static_cast<Real>(mesh.getCellCount()));
        if (params.deltaW <= 0)
          params.deltaW = params.hRef;
        if (params.normalizer <= 0)
        {
          const Real weightedBandMeasure = computeWeightedBandMeasure(psi, params.deltaW);
          if (weightedBandMeasure <= 0)
            throw std::runtime_error("LSR: weighted level-set band is empty.");
          params.normalizer = Real(1) / (weightedBandMeasure * params.hRef * params.hRef);
        }
        if (params.interfaceWeight != Real(0)
            && params.interfaceNormalizer <= 0)
        {
          const Real interfaceMeasure =
            computeInterfaceMeasure(params.interfaceAttribute);
          if (interfaceMeasure <= 0)
            throw std::runtime_error("LSR: interface-distance term is enabled but the interface set is empty.");
          params.interfaceNormalizer =
            Real(1) / (interfaceMeasure * params.hRef * params.hRef);
        }
      }

      Real initialGuessScaleFactor(const LSRParameters& params) const
      {
        switch (params.initialGuessScaling)
        {
          case LSRInitialGuessScaling::Unnormalized:
            return Real(1);
          case LSRInitialGuessScaling::EnergyNormalized:
            return params.normalizer;
          case LSRInitialGuessScaling::BandNormalized:
            return params.normalizer * params.hRef * params.hRef;
        }
        return Real(1);
      }

      template <class Psi, class PhiDerived, class GradDerived,
                class ShapeWeightDerived>
      void applyHilbertInitialGuess(
          const Psi& psi,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& gradPhi,
          const Variational::RealFunctionBase<ShapeWeightDerived>& shapeWeight,
          const BarrierParameters& barrierParams,
          const LSRParameters& params)
      {
        using Variational::DirichletBC;
        using Variational::Integral;
        using Variational::Jacobian;
        using Variational::LinearElasticityIntegral;
        using Variational::Problem;
        using Variational::RealFunction;
        using Variational::TestFunction;
        using Variational::TrialFunction;
        using Variational::VectorFunction;
        using Variational::Zero;

        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();

        TrialFunction duH(fes);
        TestFunction vH(fes);
        auto zero = VectorFunction{ Zero(), Zero() };
        const Real scale = initialGuessScaleFactor(params);

        RealFunction rhsScalar(
            [&](const Geometry::Point& p) -> Real
            {
              const Real phiVal = phi.getValue(p);
              const Real sVal = psi.getValue(p);
              const Real W =
                std::exp(-sVal * sVal / (2 * params.deltaW * params.deltaW));
              return scale * params.rhoS * W * (phiVal - sVal);
            });

        u.getData().setZero();

        Problem hilbert(duH, vH);
        switch (params.initialGuessMetric)
        {
          case LSRHilbertMetric::Harmonic:
            if (params.bandHilbertStiffness == Real(0))
            {
              hilbert =
                  Integral(Jacobian(duH), Jacobian(vH))
                + Integral(rhsScalar * gradPhi, vH)
                + DirichletBC(duH, zero);
            }
            else
            {
              // Band-weighted Harmonic metric:
              //   a(u, v) = int (1 + s_H * (1 - W(psi))^2) grad u : grad v.
              const Real sH = params.bandHilbertStiffness;
              const Real deltaW = params.deltaW;
              RealFunction muBand(
                  [&, sH, deltaW](const Geometry::Point& p) -> Real
                  {
                    const Real sVal = psi.getValue(p);
                    const Real W = std::exp(
                        -sVal * sVal / (Real(2) * deltaW * deltaW));
                    const Real c = Real(1) - W;
                    return Real(1) + sH * c * c;
                  });
              hilbert =
                  Integral(muBand * Jacobian(duH), Jacobian(vH))
                + Integral(rhsScalar * gradPhi, vH)
                + DirichletBC(duH, zero);
            }
            break;
          case LSRHilbertMetric::Elasticity:
            hilbert =
                LinearElasticityIntegral(duH, vH)(
                    params.initialGuessElasticityLambda,
                    params.initialGuessElasticityMu)
              + Integral(rhsScalar * gradPhi, vH)
              + DirichletBC(duH, zero);
            break;
          case LSRHilbertMetric::ShapeHessian:
          {
            JacobianAdmissibilityBarrierSampled barrier(
                duH, vH, u, params.quadratureOrder);
            hilbert =
                barrier.TangentPSDProjected(shapeWeight, barrierParams)
              + Integral(rhsScalar * gradPhi, vH)
              + DirichletBC(duH, zero);
            break;
          }
        }
        Solver::SparseLU(hilbert).solve();

        const Math::Vector<Real> lift = duH.getSolution().getData();
        const Math::Vector<Real> zeroData = u.getData();
        const Real jLineSearchRatio =
          std::max(params.jMinRatio,
                   params.lineSearchSafetyMargin * params.jSafeRatio);

        Real alpha = params.alphaInit;
        std::size_t backtracks = 0;
        while (alpha >= params.alphaMin)
        {
          const Math::Vector<Real> trial = alpha * lift;
          const LSRAdmissibilityReport adm =
            evaluateLSRAdmissibilitySampled(
                u, trial, params.jMinRatio, params.quadratureOrder);

          if (adm.minJRatio > jLineSearchRatio
              && adm.inadmissibleCount == 0)
          {
            u.getData() = trial;
            m_report.initialGuessAlpha = alpha;
            m_report.initialGuessBacktracks = backtracks;
            m_report.initialGuessMinJRatio = adm.minJRatio;
            return;
          }

          alpha *= params.alphaReduction;
          ++backtracks;
          u.getData() = zeroData;
        }

        u.getData() = zeroData;
        m_report.initialGuessAlpha = Real(0);
        m_report.initialGuessBacktracks = backtracks;
        m_report.initialGuessMinJRatio = Real(1);
      }

      std::reference_wrapper<Displacement> m_u;
      LSRParameters m_params;
      LSRReport m_report;
      Real m_prevAcceptedAlpha = Real(0);
  };

  template <class Displacement>
  LSR(Displacement&) -> LSR<Displacement>;
}

#endif
