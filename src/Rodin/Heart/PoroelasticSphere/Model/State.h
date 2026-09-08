/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file State.h
 * @brief State, input and intermediate data structures for the 0D
 * solid-incompressible poroelastic sphere.
 *
 * The wall is the variationally consistent 0D reduction of the thick
 * poroelastic spherical shell (uniform Lagrangian porosity, exact pointwise
 * solid incompressibility, Galerkin projection on the tangent fields
 * @f$ v_1, v_2 @f$). The passive, viscous and active laws, the valves and the
 * Windkessel are those of the CCMLC2014 model
 * (Caruel et al., Biomech Model Mechanobiol 2014, doi:10.1007/s10237-013-0544-6),
 * so the CCMLC2014 parameter sets carry over; in the membrane limit
 * @f$ d_0 / R_0 \to 0 @f$, @f$ J = 1 @f$ the wall law reduces to the CCMLC2014
 * one.
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_MODEL_STATE_H
#define RODIN_HEART_POROELASTICSPHERE_MODEL_STATE_H

#include <cstddef>
#include <functional>

#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/PartialPivLU.h"
#include "Rodin/Heart/CCMLC2014/Model/State.h"

namespace Rodin::Heart::PoroelasticSphere::Model
{
  /**
   * @brief Indices of the nonlinear unknown vector for one time step.
   *
   * The wall is quasi-static: the CCMLC2014 pair @f$ (y, v) @f$ collapses to
   * the endocardial displacement @f$ y = r_{in} - R_{in} @f$, and the
   * Lagrangian porosity @f$ \Phi @f$ is added as a differential unknown.
   */
  enum Variable : size_t
  {
    RadialDisplacement = 0, ///< @f$ y_{n+1} @f$ — endocardial radial displacement.
    Porosity,               ///< @f$ \Phi^{n+1} @f$ — Lagrangian porosity (uniform).
    VentricularPressure,    ///< @f$ p_v^{n+1} @f$ — left-ventricular pressure.
    ArterialPressure,       ///< @f$ p_{ar}^{n+1} @f$ — proximal arterial pressure.
    DistalPressure,         ///< @f$ p_d^{n+1} @f$ — distal pressure.
    FiberDeformation, ///< @f$ e_c^{n+1} @f$ — contractile deformation.
    ActiveStiffness, ///< @f$ k_c^{n+1} = \gamma^2 @f$.
    ActiveStress, ///< @f$ \tau_c^{n+1} = \gamma\beta @f$.
    LoadDependentRelaxation, ///< @f$ w^{n+1} @f$.
    NumberOfVariables       ///< Total number of unknowns in the coupled system.
  };

  /// @brief Implicit time-integration scheme (shared with CCMLC2014).
  using TimeScheme = CCMLC2014::Model::TimeScheme;

  /// @brief Rheology of the Windkessel branch flows (shared with CCMLC2014).
  using WindkesselRheology = CCMLC2014::Model::WindkesselRheology;

  /**
   * @brief Dynamic state variables at a given time.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  struct StateT
  {
    Scalar y = 0.0;     ///< Endocardial radial displacement @f$ y = r_{in} - R_{in} @f$.
    Scalar phi = 0.0;   ///< Lagrangian porosity @f$ \Phi @f$ (fluid volume per unit reference volume).
    Scalar pv = 0.0;    ///< Left-ventricular pressure @f$ p_v @f$.
    Scalar par = 0.0;   ///< Arterial/Windkessel proximal pressure @f$ p_{ar} @f$.
    Scalar pd = 0.0;    ///< Distal pressure @f$ p_d @f$.

    Scalar ec = 0.0;    ///< Contractile internal variable @f$ e_c @f$ (fiber deformation).
    Scalar kc = 0.0;    ///< Active stiffness-like scalar @f$ k_c = \gamma^2 @f$.
    Scalar tauc = 0.0;  ///< Active stress-like scalar @f$ \tau_c = \gamma\beta @f$.
    Scalar gamma = 0.0; ///< Active state variable @f$ \gamma @f$.
    Scalar beta = 0.0;  ///< Active state variable @f$ \beta @f$.
    Scalar w = 1.0; ///< Load-dependent relaxation multiplier @f$ w @f$.

    Scalar lambdaBar = 0.0; ///< Wall-averaged incompressibility multiplier @f$ \bar\lambda @f$ (explicit).
    Scalar pf = 0.0;        ///< Interstitial fluid pressure @f$ \tilde p = \Psi_P'(\Phi) - \bar\lambda @f$ (explicit).

    Scalar t = 0.0;     ///< Time associated with this state.
  };

  /**
   * @brief Model parameters and external forcings for one simulation.
   *
   * @tparam Scalar Scalar numeric type.
   * @tparam PassiveEnergyLaw Passive constitutive energy law type.
   */
  template <class Scalar, class PassiveEnergyLaw>
  struct InputT
  {
    Scalar d0 = 1.0;      ///< Reference wall thickness @f$ R_{out} - R_{in} @f$.
    Scalar R0 = 1.0;      ///< Reference cavity (endocardial) radius @f$ R_{in} @f$.
    Scalar phi0 = Scalar(0.1); ///< Reference porosity @f$ \phi_0 @f$.

    Scalar Es = 1.0;      ///< Elastic stiffness parameter.
    Scalar eta = 0.0;     ///< Viscous coefficient.
    Scalar mu = 0.0;      ///< Active-viscous coupling coefficient.
    Scalar alpha = 0.0;   ///< Active length-rate coupling coefficient.
    Scalar alphaR =
      Scalar(0.12); ///< Relaxation time scale @f$ \alpha_r @f$ for @f$ w @f$.
    Scalar k0 = 0.0;      ///< Active stiffness rate coefficient.
    Scalar sigma0 = 0.0;  ///< Active stress rate coefficient.

    Scalar KPhi = 0.0;     ///< Fluid storage modulus: @f$ \Psi_P(\Phi) = \tfrac12 K_\Phi (\Phi - \phi_0)^2 @f$ (0: pure Terzaghi).
    Scalar gammaAr = 0.0;  ///< Arterial perfusion conductance @f$ \gamma_{ar} @f$.
    Scalar gammaVen = 0.0; ///< Venous perfusion conductance @f$ \gamma_{ven} @f$.

    Scalar Cp = 1.0;      ///< Proximal (arterial) compliance.
    Scalar Cd = 1.0;      ///< Distal compliance.
    Scalar Rp = 1.0;      ///< Peripheral/proximal resistance.
    Scalar Rd = 1.0;      ///< Distal resistance.

    Scalar proximalRadius = 0.0; ///< Proximal branch radius.
    Scalar proximalLength = 0.0; ///< Proximal branch length.
    Scalar distalRadius = 0.0; ///< Distal branch radius.
    Scalar distalLength = 0.0; ///< Distal branch length.

    Scalar m = 0.0; ///< Power-law consistency coefficient.
    Scalar n = 0.0; ///< Power-law, Cross, or Carreau-Yasuda exponent.

    Scalar mu_0 = 0.0; ///< Zero-shear viscosity.
    Scalar mu_Inf = 0.0; ///< Infinite-shear viscosity.
    Scalar lambda = 0.0; ///< Rheology time-scale parameter.
    Scalar yasuda = 0.0; ///< Carreau-Yasuda transition parameter.

    Scalar k_0 = 0.0; ///< Quemada low-shear coefficient.
    Scalar k_Inf = 0.0; ///< Quemada high-shear coefficient.
    Scalar phi_quemada = 0.0; ///< Quemada hematocrit volume fraction.
    Scalar mu_plasma = 0.0; ///< Plasma viscosity for Quemada rheology.
    Scalar gamma_c = 0.0; ///< Critical shear rate for Quemada rheology.

    Scalar Kat = 0.0;     ///< Atrioventricular (mitral) conductance.
    Scalar Kp = 0.0;      ///< Valve leakage conductance.
    Scalar Kar = 0.0;     ///< Aortic conductance.

    Scalar cavityCapacity = Scalar(5e-12); ///< Cavity compliance-like capacity term.

    Scalar absRegularization = Scalar(1e-14); ///< Regularization scalar for absolute values.
    size_t wallQuadraturePoints = 8; ///< Gauss-Legendre points across the wall @f$ [R_{in}, R_{out}] @f$.
    TimeScheme timeScheme =
      TimeScheme::BDF2; ///< Implicit time discretization for dynamic states.
    WindkesselRheology windkesselRheology =
      WindkesselRheology::Newtonian; ///< Rheology for reduced Windkessel branches.

    Scalar initFibDef = 0.0;          ///< Initial fiber deformation @f$ e_c @f$.
    Scalar initActiveStiffness = 0.0; ///< Initial active stiffness-like value @f$ k_c @f$.
    Scalar initActiveStress = 0.0;    ///< Initial active stress-like value @f$ \tau_c @f$.

    std::function<Scalar(Scalar)> u =
      [](Scalar) { return Scalar(0); }; ///< Active input drive @f$ u(t) @f$.

    std::function<Scalar(Scalar)> m0 = [](Scalar) {
      return Scalar(1);
    }; ///< Relaxation target @f$ m_0(e_c) @f$.

    std::function<Scalar(Scalar)> dm0 = [](Scalar) {
      return Scalar(0);
    }; ///< Derivative @f$ m_0'(e_c) @f$.

    std::function<Scalar(Scalar)> pAt =
      [](Scalar) { return Scalar(0); }; ///< Atrial pressure boundary condition.

    std::function<Scalar(Scalar)> pSv =
      [](Scalar) { return Scalar(0); }; ///< Venous pressure boundary condition.

    std::function<Scalar(Scalar)> qArterialExternal = [](Scalar) {
      return Scalar(0);
    }; ///< Flow drawn from the proximal Windkessel node by an external coronary model (e.g. a 3D tree), evaluated at @f$ t_{n+1} @f$.

    std::function<Scalar(Scalar)> qPerfusionExternal = [](Scalar) {
      return Scalar(0);
    }; ///< Flow delivered into the interstitium by an external coronary model, evaluated at @f$ t_{n+1} @f$.

    PassiveEnergyLaw passiveEnergy; ///< Passive reduced constitutive law.
  };

  /**
   * @brief Newton report for one global dynamic step.
   *
   * @tparam Scalar Scalar numeric type.
   * @tparam DenseLinearSystem Dense linear system type.
   */
  template <class Scalar, class DenseLinearSystem>
  struct ReportT
  {
    bool converged = false;     ///< Whether the nonlinear solve converged.
    size_t iterations = 0;      ///< Number of nonlinear iterations used.
    Scalar finalResidual = 0.0; ///< Final residual norm.
    Scalar finalStepNorm = 0.0; ///< Final step norm.
    /**
     * @brief Termination reason of the nonlinear solve.
     *
     * Fully qualified: inside Rodin::Heart::PoroelasticSphere the name
     * Solver:: resolves to the nested PoroelasticSphere::Solver namespace,
     * not Rodin::Solver.
     */
    using ConvergenceReason = typename ::Rodin::Solver::NewtonSolver<
      ::Rodin::Solver::PartialPivLU<DenseLinearSystem>>::ConvergedReason;
    ConvergenceReason reason = ConvergenceReason::MaxIterations; ///< Nonlinear convergence reason.
  };

  /**
   * @brief Active-law data cached during residual/Jacobian evaluation.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  struct ActiveDataT
  {
    Scalar activationDrive = 0.0;             ///< Value of @f$ u(t_{n+1}) @f$.
    Scalar activationDrivePositivePart = 0.0; ///< Positive part of activation drive.
    Scalar activationDriveNegativePart =
      0.0; ///< Negative part magnitude of activation drive.
    Scalar relaxationTarget = 1.0; ///< Target @f$ m_0(e_c) @f$ for @f$ w @f$.
    Scalar relaxationDrive = 0.0; ///< Effective decay drive @f$ |u|_+ + w |u|_- + \alpha |\dot e_c| @f$.
    Scalar recruitmentFraction = 0.0;         ///< Recruitment fraction @f$ n_0 @f$.

    Scalar activeStressOneDimensional = 0.0;             ///< Series-element fiber stress @f$ \sigma_{1D} @f$.
    Scalar partialActiveStressWrtDisplacement = 0.0;     ///< @f$ \partial \sigma_{1D} / \partial y @f$.
    Scalar partialActiveStressWrtPorosity = 0.0;         ///< @f$ \partial \sigma_{1D} / \partial \Phi @f$.
    Scalar partialActiveStressWrtFiberDeformation = 0.0; ///< @f$ \partial \sigma_{1D} / \partial e_c @f$.
  };

  /**
   * @brief Intermediate data cached during residual/Jacobian evaluation.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  struct EvalDataT
  {
    StateT<Scalar> sn;   ///< State at time level @f$ n @f$.
    StateT<Scalar> snm1; ///< State at time level @f$ n-1 @f$.

    Scalar tnp1 = 0.0; ///< Target time @f$ t_{n+1} @f$.
    Scalar dt = 0.0;   ///< Time-step size.

    Scalar y = 0.0;   ///< Candidate @f$ y_{n+1} @f$.
    Scalar phi = 0.0; ///< Candidate @f$ \Phi^{n+1} @f$.
    Scalar pv = 0.0;  ///< Candidate @f$ p_v^{n+1} @f$.
    Scalar par = 0.0; ///< Candidate @f$ p_{ar}^{n+1} @f$.
    Scalar pd = 0.0;  ///< Candidate @f$ p_d^{n+1} @f$.
    Scalar ec = 0.0; ///< Candidate @f$ e_c^{n+1} @f$.
    Scalar kc = 0.0; ///< Candidate @f$ k_c^{n+1} @f$.
    Scalar tauc = 0.0; ///< Candidate @f$ \tau_c^{n+1} @f$.
    Scalar w = 1.0; ///< Candidate @f$ w^{n+1} @f$.

    Scalar yDot = 0.0;   ///< Discrete rate @f$ \dot y @f$ (BE/BDF2).
    Scalar phiDot = 0.0; ///< Discrete rate @f$ \dot\Phi @f$ (BE/BDF2).
    Scalar a0 = 0.0;     ///< Leading time-derivative coefficient @f$ \partial \dot q / \partial q @f$.

    Scalar pvMid = 0.0;  ///< Pressure seen by the Windkessel evaluator (fully implicit).
    Scalar parMid = 0.0; ///< Arterial pressure seen by the Windkessel evaluator.
    Scalar pdMid = 0.0;  ///< Distal pressure seen by the Windkessel evaluator.

    Scalar rin = 0.0; ///< Deformed endocardial radius @f$ r_{in} = R_{in} + y @f$.
    Scalar J = 1.0;   ///< Uniform volume ratio @f$ J = 1 - \phi_0 + \Phi @f$.
    Scalar Vw0 = 0.0; ///< Reference wall volume @f$ V_{w0} = \tfrac{4\pi}{3}(R_{out}^3 - R_{in}^3) @f$.

    Scalar lambdaMid = 1.0;  ///< Mid-wall hoop stretch @f$ \lambda_\theta(R_m) @f$.
    Scalar strain1D = 0.0;   ///< Mid-wall Green-Lagrange hoop strain @f$ e = \tfrac12(\lambda_\theta^2 - 1) @f$.
    Scalar diffGreen = 0.0;  ///< @f$ \partial e / \partial y @f$.
    Scalar diffGreenWrtPorosity = 0.0; ///< @f$ \partial e / \partial \Phi @f$.

    Scalar wallPressure = 0.0;      ///< Cavity pressure sustained by the wall, @f$ P_v(r_{in}, \Phi) @f$ from the @f$ v_1 @f$ projection.
    Scalar dWallPressure_dy = 0.0;  ///< Total derivative wrt @f$ y @f$ (rates included).
    Scalar dWallPressure_dphi = 0.0; ///< Total derivative wrt @f$ \Phi @f$ (rates included).
    Scalar dWallPressure_dec = 0.0; ///< Derivative wrt @f$ e_c @f$.

    Scalar lambdaBar = 0.0;      ///< Wall-averaged multiplier @f$ \bar\lambda @f$ from the @f$ v_2 @f$ projection.
    Scalar dLambdaBar_dy = 0.0;  ///< Total derivative wrt @f$ y @f$.
    Scalar dLambdaBar_dphi = 0.0; ///< Total derivative wrt @f$ \Phi @f$.
    Scalar dLambdaBar_dec = 0.0; ///< Derivative wrt @f$ e_c @f$.

    Scalar pf = 0.0;      ///< Interstitial pressure @f$ \tilde p = K_\Phi(\Phi - \phi_0) - \bar\lambda @f$.
    Scalar dPf_dy = 0.0;  ///< Total derivative wrt @f$ y @f$.
    Scalar dPf_dphi = 0.0; ///< Total derivative wrt @f$ \Phi @f$.
    Scalar dPf_dec = 0.0; ///< Derivative wrt @f$ e_c @f$.

    Scalar perfusionInflow = 0.0;  ///< Arterial perfusion inflow @f$ \gamma_{ar}(p_{ar} - \tilde p) @f$.
    Scalar perfusionOutflow = 0.0; ///< Venous perfusion outflow @f$ \gamma_{ven}(\tilde p - p_{sv}) @f$.
    Scalar externalArterialOutflow = 0.0; ///< External flow drawn from the proximal Windkessel node (@ref InputT::qArterialExternal).
    Scalar externalPerfusionInflow = 0.0; ///< External flow delivered into the interstitium (@ref InputT::qPerfusionExternal).

    Scalar pAtCur = 0.0;  ///< Atrial pressure at @f$ t_{n+1} @f$.
    Scalar pSvMid = 0.0;  ///< Venous pressure at @f$ t_{n+1} @f$.

    Scalar cavityFluxCur = 0.0;      ///< Current cavity valve flux.
    Scalar dCavityFluxCur_dPv = 0.0; ///< Derivative of current cavity flux wrt ventricular pressure.
    Scalar dCavityFluxCur_dPar = 0.0; ///< Derivative of current cavity flux wrt arterial pressure.

    Scalar windkesselOutflow = 0.0;      ///< Outflow term toward Windkessel branch.
    Scalar dWindkesselOutflow_dPv = 0.0; ///< Derivative of outflow wrt ventricular pressure.
    Scalar dWindkesselOutflow_dPar = 0.0; ///< Derivative of outflow wrt arterial pressure.

    Scalar windkesselflowP = 0.0; ///< flow term toward Windkessel branch.
    Scalar windkesselflowD = 0.0; ///< flow term toward Windkessel branch.
    Scalar dWindkesselflowP_dPar = 0.0; ///< Derivative of flow wrt arterial pressure.
    Scalar dWindkesselflowP_dPd = 0.0; ///< Derivative of flow wrt distal pressure.
    Scalar dWindkesselflowD_dPar = 0.0; ///< Derivative of flow wrt arterial pressure.
    Scalar dWindkesselflowD_dPd = 0.0; ///< Derivative of flow wrt distal pressure.

    ActiveDataT<Scalar> active; ///< Active-law data.
  };
}

#endif
