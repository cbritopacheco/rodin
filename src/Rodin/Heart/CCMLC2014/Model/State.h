/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file State.h
 * @brief State, input and intermediate data structures for the CCMLC2014 0D model.
 *
 * The variable names documented in this file follow the reduced 0D notation from:
 * M. Caruel et al., Biomechanics and Modeling in Mechanobiology (2014),
 * doi:10.1007/s10237-013-0544-6, HAL: hal-00872746.
 */
#ifndef RODIN_HEART_CCMLC2014_MODEL_STATE_H
#define RODIN_HEART_CCMLC2014_MODEL_STATE_H

#include <cstddef>
#include <functional>

#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/PartialPivLU.h"

namespace Rodin::Heart::CCMLC2014::Model
{
  /**
   * @brief Indices of the nonlinear unknown vector for one time step.
   *
   * Matches the ordering of the coupled system unknowns
   * @f$ (y, p_v, p_{ar}, p_d) @f$ from Caruel et al. (2014),
   * Biomechanics and Modeling in Mechanobiology, §4 (0D model).
   */
  enum Variable : size_t
  {
    RadialDisplacement = 0, ///< @f$ y_{n+1} @f$ — radial displacement.
    RadialVelocity,         ///< @f$ v_{n+1} @f$ — radial velocity.
    VentricularPressure,    ///< @f$ p_v^{n+1} @f$ — left-ventricular pressure.
    ArterialPressure,       ///< @f$ p_{ar}^{n+1} @f$ — proximal arterial pressure.
    DistalPressure,         ///< @f$ p_d^{n+1} @f$ — distal pressure.
    FiberDeformation,       ///< @f$ e_c^{n+1} @f$ — contractile deformation.
    ActiveStiffness,        ///< @f$ k_c^{n+1} = \gamma^2 @f$.
    ActiveStress,           ///< @f$ \tau_c^{n+1} = \gamma\beta @f$.
    LoadDependentRelaxation, ///< @f$ w^{n+1} @f$.
    NumberOfVariables       ///< Total number of unknowns in the coupled system.
  };

  /**
   * @brief Implicit time-integration scheme for the 0D state equations.
   */
  enum class TimeScheme
  {
    BackwardEuler, ///< First-order L-stable backward Euler.
    BDF2           ///< Second-order stiffly-damped BDF2 with BE startup.
  };

  /**
   * @brief Rheology used in the reduced Windkessel branch flow laws.
   */
  enum class WindkesselRheology
  {
    Newtonian,              ///< Linear resistance branch flows.
    CarreauYasuda,           ///< Non-Newtonian Carreau-Yasuda tube flow.
    Cross,
    PowerLaw,
    Quemada
  };

  /**
   * @brief CCMLC2014 dynamic state variables at a given time.
   *
   * The short symbols are intentionally retained to match the publication notation.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  struct StateT
  {
    Scalar y = 0.0;     ///< Radial displacement-like kinematic unknown @f$ y @f$.
    Scalar v = 0.0;     ///< Radial velocity-like quantity @f$ \dot{y} @f$.
    Scalar pv = 0.0;    ///< Left-ventricular pressure @f$ p_v @f$.
    Scalar par = 0.0;   ///< Arterial/Windkessel proximal pressure @f$ p_{ar} @f$.
    Scalar pd = 0.0;    ///< Distal pressure @f$ p_d @f$.

    Scalar ec = 0.0;    ///< Contractile internal variable @f$ e_c @f$ (fiber deformation).
    Scalar kc = 0.0;    ///< Active stiffness-like scalar @f$ k_c = \gamma^2 @f$.
    Scalar tauc = 0.0;  ///< Active stress-like scalar @f$ \tau_c = \gamma\beta @f$.
    Scalar gamma = 0.0; ///< Active state variable @f$ \gamma @f$.
    Scalar beta = 0.0;  ///< Active state variable @f$ \beta @f$.
    Scalar w = 1.0;     ///< Load-dependent relaxation multiplier @f$ w @f$.

    Scalar t = 0.0;     ///< Time associated with this state.
  };

  /**
   * @brief Model parameters and external forcings for one CCMLC2014 simulation.
   *
   * @tparam Scalar Scalar numeric type.
   * @tparam PassiveEnergyLaw Passive constitutive energy law type.
   */
  template <class Scalar, class PassiveEnergyLaw>
  struct InputT
  {
    Scalar rho = 1.0;     ///< Wall density/inertia coefficient.
    Scalar d0 = 1.0;      ///< Reference wall thickness.
    Scalar R0 = 1.0;      ///< Reference cavity radius.

    Scalar Es = 1.0;      ///< Elastic stiffness parameter.
    Scalar eta = 0.0;     ///< Viscous coefficient.
    Scalar mu = 0.0;      ///< Active-viscous coupling coefficient.
    Scalar alpha = 0.0;   ///< Active length-rate coupling coefficient.
    Scalar alphaR = Scalar(0.12); ///< Relaxation time scale @f$ \alpha_r @f$ for @f$ w @f$.
    Scalar k0 = 0.0;      ///< Active stiffness rate coefficient.
    Scalar sigma0 = 0.0;  ///< Active stress rate coefficient.

    Scalar Cp = 1.0;      ///< Proximal (arterial) compliance.
    Scalar Cd = 1.0;      ///< Distal compliance.
    Scalar Rp = 1.0;      ///< Peripheral/proximal resistance.
    Scalar Rd = 1.0;      ///< Distal resistance.

    Scalar proximalRadius = 0.0;
    Scalar proximalLength = 0.0;
    Scalar distalRadius = 0.0;
    Scalar distalLength = 0.0;

    Scalar m = 0.0;
    Scalar n = 0.0;

    Scalar mu_0    = 0.0;
    Scalar mu_Inf  = 0.0;
    Scalar lambda = 0.0;
    Scalar yasuda = 0.0;

    Scalar k_0 = 0.0;
    Scalar k_Inf = 0.0;
    Scalar phi_quemada = 0.0;
    Scalar mu_plasma = 0.0;
    Scalar gamma_c = 0.0;

    Scalar Kat = 0.0;     ///< Atrioventricular (mitral) conductance.
    Scalar Kp = 0.0;      ///< Valve leakage conductance.
    Scalar Kar = 0.0;     ///< Aortic conductance.

    Scalar cavityCapacity = Scalar(5e-12); ///< Cavity compliance-like capacity term.

    Scalar localTolerance = Scalar(1e-12); ///< Absolute tolerance for local active solve.
    size_t localMaxIterations = 50;         ///< Maximum iterations for local active solve.
    Scalar localDamping = Scalar(1.0);      ///< Damping factor for local active Newton updates.
    Scalar absRegularization = Scalar(1e-14); ///< Regularization scalar for absolute values.
    TimeScheme timeScheme = TimeScheme::BDF2; ///< Implicit time discretization for dynamic states.
    WindkesselRheology windkesselRheology =
      WindkesselRheology::Newtonian; ///< Rheology for reduced Windkessel branches.

    Scalar initFibDef = 0.0;          ///< Initial fiber deformation @f$ e_c @f$.
    Scalar initActiveStiffness = 0.0; ///< Initial active stiffness-like value @f$ k_c @f$.
    Scalar initActiveStress = 0.0;    ///< Initial active stress-like value @f$ \tau_c @f$.

    std::function<Scalar(Scalar)> u =
      [](Scalar) { return Scalar(0); }; ///< Active input drive @f$ u(t) @f$.

    std::function<Scalar(Scalar)> m0 =
      [](Scalar) { return Scalar(1); }; ///< Relaxation target @f$ m_0(e_c) @f$.

    std::function<Scalar(Scalar)> dm0 =
      [](Scalar) { return Scalar(0); }; ///< Derivative @f$ m_0'(e_c) @f$.

    std::function<Scalar(Scalar)> pAt =
      [](Scalar) { return Scalar(0); }; ///< Atrial pressure boundary condition.

    std::function<Scalar(Scalar)> pSv =
      [](Scalar) { return Scalar(0); }; ///< Venous pressure boundary condition.

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
    using ConvergenceReason =
      typename Solver::NewtonSolver<::Rodin::Solver::PartialPivLU<DenseLinearSystem>>::ConvergedReason;
    ConvergenceReason reason = ConvergenceReason::MaxIterations; ///< Nonlinear convergence reason.
  };

  /**
   * @brief Internal data for the local active-dynamics solve at one global step.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  struct LocalActiveDataT
  {
    Scalar fiberDeformationPrevious = 0.0; ///< Previous local fiber deformation.
    Scalar fiberDeformationCurrent = 0.0;  ///< Current local fiber deformation iterate.
    Scalar fiberDeformationMidpoint = 0.0; ///< Midpoint local fiber deformation.

    Scalar gammaPrevious = 0.0; ///< Previous global gamma state.
    Scalar betaPrevious = 0.0;  ///< Previous global beta state.
    Scalar wPrevious = 1.0;     ///< Previous load-dependent relaxation state.
    Scalar gammaCurrent = 0.0;  ///< Updated gamma state.
    Scalar betaCurrent = 0.0;   ///< Updated beta state.
    Scalar wCurrent = 1.0;      ///< Updated load-dependent relaxation state.

    Scalar activationDrive = 0.0;             ///< Value of @f$ u(t_{n+1}) @f$.
    Scalar activationDrivePositivePart = 0.0; ///< Positive part of activation drive.
    Scalar activationDriveNegativePart = 0.0; ///< Negative part magnitude of activation drive.
    Scalar relaxationTarget = 1.0;            ///< Target @f$ m_0(e_c^n) @f$ for @f$ w @f$.
    Scalar relaxationDrive = 0.0;             ///< Effective decay drive @f$ |u|_+ + w |u|_- @f$.
    Scalar recruitmentFraction = 0.0;         ///< Recruitment fraction @f$ n_0 @f$.

    Scalar partialResidualWrtDisplacement = 0.0;     ///< Coupling residual derivative wrt displacement.
    Scalar partialResidualWrtFiberDeformation = 0.0; ///< Local residual derivative wrt fiber deformation.
    Scalar fiberDeformationNewtonStep = 0.0;         ///< Newton correction term for fiber deformation.

    Scalar activeStressOneDimensional = 0.0;             ///< Local one-dimensional active stress.
    Scalar partialActiveStressWrtDisplacement = 0.0;     ///< Derivative wrt displacement.
    Scalar partialActiveStressWrtFiberDeformation = 0.0; ///< Derivative wrt fiber deformation.

    Scalar activeStress = 0.0;                   ///< Effective active stress in global balance.
    Scalar dActiveStressWrtDisplacement = 0.0;   ///< Effective active-stress derivative wrt displacement.

    bool converged = false; ///< Whether local active solve converged.
    size_t iterations = 0;  ///< Number of local active iterations.
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
    Scalar v = 0.0;   ///< Candidate @f$ v_{n+1} @f$.
    Scalar pv = 0.0;  ///< Candidate @f$ p_v^{n+1} @f$.
    Scalar par = 0.0; ///< Candidate @f$ p_{ar}^{n+1} @f$.
    Scalar pd = 0.0;  ///< Candidate @f$ p_d^{n+1} @f$.
    Scalar ec = 0.0;  ///< Candidate @f$ e_c^{n+1} @f$.
    Scalar kc = 0.0;  ///< Candidate @f$ k_c^{n+1} @f$.
    Scalar tauc = 0.0; ///< Candidate @f$ \tau_c^{n+1} @f$.
    Scalar w = 1.0;   ///< Candidate @f$ w^{n+1} @f$.

    Scalar yPrev = 0.0;   ///< @f$ y_n @f$.
    Scalar vPrev = 0.0;   ///< @f$ v_n @f$.
    Scalar pvPrev = 0.0;  ///< @f$ p_v^n @f$.
    Scalar parPrev = 0.0; ///< @f$ p_{ar}^n @f$.
    Scalar pdPrev = 0.0;  ///< @f$ p_d^n @f$.

    Scalar yPrevPrev = 0.0; ///< @f$ y_{n-1} @f$.

    Scalar yMid = 0.0;   ///< Midpoint displacement.
    Scalar pvMid = 0.0;  ///< Midpoint ventricular pressure.
    Scalar parMid = 0.0; ///< Midpoint arterial pressure.
    Scalar pdMid = 0.0;  ///< Midpoint distal pressure.
    Scalar vel = 0.0;    ///< Midpoint radial velocity approximation.

    Scalar sqrtC = 0.0;    ///< @f$ \sqrt{C} @f$.
    Scalar C = 0.0;        ///< Reduced right Cauchy-Green scalar @f$ C @f$.
    Scalar strain1D = 0.0; ///< Reduced Green-Lagrange strain.
    Scalar diffGreen = 0.0; ///< Derivative of reduced strain mapping wrt displacement.

    Scalar stressPassive = 0.0;      ///< Passive stress contribution.
    Scalar diffStressPassive = 0.0;  ///< Derivative of passive stress wrt displacement.
    Scalar stressViscous = 0.0;      ///< Viscous stress contribution.
    Scalar diffStressViscous = 0.0;  ///< Derivative of viscous stress wrt displacement.
    Scalar diffStressViscousWrtVelocity = 0.0; ///< Derivative of viscous stress wrt velocity.

    Scalar pAtCur = 0.0;  ///< Atrial pressure at @f$ t_{n+1} @f$.
    Scalar pAtPrev = 0.0; ///< Atrial pressure at @f$ t_n @f$.
    Scalar pSvMid = 0.0;  ///< Venous pressure at midpoint time.

    Scalar cavityFluxCur = 0.0;      ///< Current cavity valve flux.
    Scalar cavityFluxPrev = 0.0;     ///< Previous cavity valve flux.
    Scalar dCavityFluxCur_dPv = 0.0; ///< Derivative of current cavity flux wrt ventricular pressure.
    Scalar dCavityFluxCur_dPar = 0.0; ///< Derivative of current cavity flux wrt arterial pressure.

    Scalar windkesselOutflow = 0.0;      ///< Outflow term toward Windkessel branch.
    Scalar dWindkesselOutflow_dPv = 0.0; ///< Derivative of outflow wrt ventricular pressure.
    Scalar dWindkesselOutflow_dPar = 0.0; ///< Derivative of outflow wrt arterial pressure.

    Scalar windkesselflowP = 0.0;    ///< flow term toward Windkessel branch.
    Scalar windkesselflowD = 0.0;    ///< flow term toward Windkessel branch.
    Scalar dWindkesselflowP_dPar = 0.0;  ///< Derivative of flow wrt ventricular pressure.
    Scalar dWindkesselflowP_dPd = 0.0;  ///< Derivative of flow wrt ventricular pressure.
    Scalar dWindkesselflowD_dPar = 0.0;  ///< Derivative of flow wrt arterial pressure.
    Scalar dWindkesselflowD_dPd = 0.0;  ///< Derivative of flow wrt arterial pressure.


    LocalActiveDataT<Scalar> active; ///< Local active-dynamics data.
  };
}

#endif
