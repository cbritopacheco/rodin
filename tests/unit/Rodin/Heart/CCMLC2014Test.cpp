#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "Rodin/Heart/CCMLC2014.h"
#include "Rodin/Heart/CCMLC2014/Numerics/DynamicSystem.h"
#include "Rodin/Heart/CCMLC2014/Physics/Kinematics.h"
#include "Rodin/Heart/CCMLC2014/Physics/PassiveStress.h"
#include "Rodin/Heart/CCMLC2014/Physics/ValveLaw.h"

using namespace Rodin;

namespace
{
  using Model = Heart::CCMLC2014T<>;
  using PassiveLaw = Heart::CCMLC2014PassiveLaw<Real>;
  namespace CCMLC2014Vars = Heart::CCMLC2014::Model;

  Model::Input makeGenericCardiacInput()
  {
    Model::Input cardiacInput;

    cardiacInput.rho = 1.0e3;
    cardiacInput.R0 = 2.36e-2;
    cardiacInput.d0 = 1.42e-2;

    cardiacInput.Es = 3.0e5;
    cardiacInput.mu = 70.0;
    cardiacInput.eta = 70.0;
    cardiacInput.alpha = 3.0;
    cardiacInput.k0 = 1.0e5;
    cardiacInput.sigma0 = 5.0e5;

    cardiacInput.Rp = 8.0e6;
    cardiacInput.Cp = 5.0e-9;
    cardiacInput.Rd = 1.0e8;
    cardiacInput.Cd = 1.0e-8;

    cardiacInput.Kat = 8.0e-7;
    cardiacInput.Kp  = 5.0e-10;
    cardiacInput.Kar = 1.3e-5;

    cardiacInput.cavityCapacity = 5.0e-12;
    cardiacInput.localTolerance = 1e-12;
    cardiacInput.localMaxIterations = 50;
    cardiacInput.localDamping = 1.0;
    cardiacInput.absRegularization = 1e-14;

    cardiacInput.u = [](Real) { return 25.0; };
    cardiacInput.pAt = [](Real) { return 900.0; };
    cardiacInput.pSv = [](Real) { return 1000.0; };

    using PassiveEnergyType = std::decay_t<decltype(cardiacInput.passiveEnergy)>;
    typename PassiveEnergyType::Parameters passiveParameters;
    passiveParameters.mu1 = 0.0;
    passiveParameters.mu2 = 0.0;
    passiveParameters.C0 = 1.9e3;
    passiveParameters.C1 = 1.1e-1;
    passiveParameters.C2 = 1.9e3;
    passiveParameters.C3 = 1.1e-1;
    cardiacInput.passiveEnergy = PassiveEnergyType(passiveParameters);

    return cardiacInput;
  }

  Model::DenseVector makeCandidateState(Real scale = 1.0)
  {
    Model::DenseVector candidateState(CCMLC2014Vars::NumberOfVariables);
    candidateState.setZero();
    candidateState[CCMLC2014Vars::RadialDisplacement] = scale * 8e-5;
    candidateState[CCMLC2014Vars::RadialVelocity] = scale * 1.5e-3;
    candidateState[CCMLC2014Vars::VentricularPressure] = scale * 1.0e4;
    candidateState[CCMLC2014Vars::ArterialPressure] = scale * 9.0e3;
    candidateState[CCMLC2014Vars::DistalPressure] = scale * 8.0e3;
    candidateState[CCMLC2014Vars::FiberDeformation] = scale * 0.012;
    candidateState[CCMLC2014Vars::ActiveStiffness] = scale * 0.045;
    candidateState[CCMLC2014Vars::ActiveStress] = scale * 0.065;
    candidateState[CCMLC2014Vars::LoadDependentRelaxation] = 1.2;
    return candidateState;
  }

  Model::State makeCurrentState(Real scale = 1.0)
  {
    Model::State currentState;
    currentState.t = 0.1;
    currentState.y = scale * 7e-5;
    currentState.v = scale * 1.0e-3;
    currentState.pv = scale * 9.8e3;
    currentState.par = scale * 8.8e3;
    currentState.pd = scale * 7.9e3;
    currentState.ec = scale * 0.01;
    currentState.kc = scale * 0.04;
    currentState.tauc = scale * 0.06;
    currentState.gamma = std::sqrt(std::max<Real>(currentState.kc, 0.0));
    currentState.beta =
      (currentState.gamma > 0.0) ? currentState.tauc / currentState.gamma : 0.0;
    currentState.w = 1.0;
    return currentState;
  }

  Model::State makePreviousState(const Model::State& currentState, Real scale = 1.0)
  {
    Model::State previousState = currentState;
    previousState.t = currentState.t - 1e-3;
    previousState.y = scale * 6e-5;
    previousState.v = scale * 8e-4;
    previousState.ec = scale * 0.009;
    previousState.kc = scale * 0.05;
    previousState.tauc = scale * 0.07;
    previousState.gamma = std::sqrt(std::max<Real>(previousState.kc, 0.0));
    previousState.beta =
      (previousState.gamma > 0.0) ? previousState.tauc / previousState.gamma : 0.0;
    previousState.w = 0.98;
    return previousState;
  }

  Real examplePeriodicActivation(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);

    if (tau < 0.13)
      return 0.0;
    if (tau < 0.141)
      return 35.0 * ((tau - 0.13) / 0.011);
    if (tau < 0.281)
      return 35.0;
    if (tau < 0.361)
      return 35.0 - 55.0 * ((tau - 0.281) / 0.08);
    if (tau < 0.45)
      return -20.0;
    return 0.0;
  }

  Real exampleLoadDependentRelaxationM0(Real ec)
  {
    const Real lowEc = 0.0;
    const Real highEc = 2.0;
    const Real lowValue = 1.6;
    const Real highValue = 1.0;

    if (ec <= lowEc)
      return lowValue;
    if (ec >= highEc)
      return highValue;

    const Real s = (ec - lowEc) / (highEc - lowEc);
    return (1.0 - s) * lowValue + s * highValue;
  }

  Real exampleLoadDependentRelaxationDM0(Real ec)
  {
    const Real lowEc = 0.0;
    const Real highEc = 2.0;
    const Real lowValue = 1.6;
    const Real highValue = 1.0;

    if (ec <= lowEc || ec >= highEc)
      return 0.0;
    return (highValue - lowValue) / (highEc - lowEc);
  }

  Real exampleAtrialPressure(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);

    const Real minValue = 500.0;
    const Real maxValue = 1000.0;
    const Real secondThreshold = 1250.0;

    const Real t1 = 0.02;
    const Real t2 = 0.15;
    const Real t3 = 0.17;
    const Real t4 = 0.56;
    const Real t5 = 0.62;
    const Real t6 = 0.85;

    Real alpha = 0.0;
    Real value = minValue;

    if (tau < t1)
    {
      alpha = -(tau - t1) / t1;
      value = alpha * minValue + (1.0 - alpha) * maxValue;
    }
    else if (tau < t2)
    {
      value = maxValue;
    }
    else if (tau < t3)
    {
      alpha = -(tau - t3) / (t3 - t2);
      value = alpha * maxValue + (1.0 - alpha) * minValue;
    }
    else if (tau < t4)
    {
      alpha = -(tau - t4) / (t4 - t3);
      value = alpha * minValue + (1.0 - alpha) * secondThreshold;
    }
    else if (tau < t5)
    {
      value = secondThreshold;
    }
    else if (tau < t6)
    {
      alpha = -(tau - t6) / (t6 - t5);
      value = alpha * secondThreshold + (1.0 - alpha) * minValue;
    }

    return value;
  }

  Model::Input makeExampleCardiacInput()
  {
    Model::Input cardiacInput;

    cardiacInput.rho = 1.0e3;
    cardiacInput.R0 = 2.36e-2;
    cardiacInput.d0 = 1.42e-2;

    cardiacInput.Es = 3.0e7;
    cardiacInput.mu = 70.0;
    cardiacInput.eta = 70.0;
    cardiacInput.alpha = 1.5;
    cardiacInput.alphaR = 0.12;
    cardiacInput.k0 = 1.0e5;
    cardiacInput.sigma0 = 1.24e5;

    cardiacInput.Rp = 8.0e6;
    cardiacInput.Cp = 2.5e-9;
    cardiacInput.Rd = 1.0e8;
    cardiacInput.Cd = 1.0e-8;

    cardiacInput.mu_0 = 0.04868;
    cardiacInput.mu_Inf = 0.003605;
    cardiacInput.lambda = 3.39;
    cardiacInput.n = 0.198;
    cardiacInput.yasuda = 1.235;
    cardiacInput.proximalRadius = 0.015;
    cardiacInput.proximalLength = 0.4;
    cardiacInput.distalRadius = 0.0007;
    cardiacInput.distalLength = 0.004;
    cardiacInput.windkesselRheology =
      Heart::CCMLC2014::Model::WindkesselRheology::CarreauYasuda;

    cardiacInput.Kat = 9.0e-6;
    cardiacInput.Kp = 5.0e-10;
    cardiacInput.Kar = 1.3e-5;

    cardiacInput.cavityCapacity = 5.0e-12;
    cardiacInput.localTolerance = 1e-12;
    cardiacInput.localMaxIterations = 50;
    cardiacInput.localDamping = 1.0;
    cardiacInput.absRegularization = 1e-14;

    cardiacInput.initFibDef = 0.0;
    cardiacInput.initActiveStiffness = 0.0;
    cardiacInput.initActiveStress = 0.0;

    cardiacInput.pSv = [](Real) { return 1.0e3; };
    cardiacInput.pAt = exampleAtrialPressure;
    cardiacInput.u = examplePeriodicActivation;
    cardiacInput.m0 = exampleLoadDependentRelaxationM0;
    cardiacInput.dm0 = exampleLoadDependentRelaxationDM0;

    using PassiveEnergyType = std::decay_t<decltype(cardiacInput.passiveEnergy)>;
    typename PassiveEnergyType::Parameters passiveParameters;
    passiveParameters.mu1 = 0.0;
    passiveParameters.mu2 = 0.0;
    passiveParameters.C0 = 1.9e3;
    passiveParameters.C1 = 1.1e-1;
    passiveParameters.C2 = 1.9e3;
    passiveParameters.C3 = 1.1e-1;
    cardiacInput.passiveEnergy = PassiveEnergyType(passiveParameters);

    return cardiacInput;
  }

  Model::State makeExampleInitialState(const Model::Input& cardiacInput)
  {
    Model::State initialState;
    initialState.t = 0.0;
    initialState.y = 0.0;
    initialState.v = 0.0;
    initialState.pv = cardiacInput.pAt(0.0) - 100.0;
    initialState.par = 11000.0;
    initialState.pd = 10000.0;
    initialState.ec = cardiacInput.initFibDef;
    initialState.gamma = std::sqrt(std::max<Real>(cardiacInput.initActiveStiffness, 0.0));
    initialState.beta = (initialState.gamma > 0.0)
      ? (cardiacInput.initActiveStress / initialState.gamma)
      : 0.0;
    initialState.kc = initialState.gamma * initialState.gamma;
    initialState.tauc = initialState.gamma * initialState.beta;
    initialState.w = cardiacInput.m0(initialState.ec);
    return initialState;
  }

  /**
   * @brief Computes Jacobian agreement error against a central finite-difference approximation.
   *
   * Builds the analytical dynamic-system Jacobian at the provided evaluation state and
   * compares it to a central finite-difference Jacobian assembled by perturbing each unknown
   * by @p relativePerturbation * max(1, |x_j|).
   *
   * @param[in] cardiacInput Model input parameters.
   * @param[in] candidateState Candidate nonlinear state vector at which to evaluate.
   * @param[in] currentState Current model state.
   * @param[in] previousState Previous model state.
   * @param[in] dt Time step used for Jacobian/residual assembly.
   * @param[in] relativePerturbation Relative perturbation factor for finite differences.
   * @returns Relative matrix error
   * @f$ \|J_{analytic} - J_{fd}\| / \max(\|J_{fd}\|, 10^{-14}) @f$.
   */
  Real computeDynamicJacobianRelativeError(
      const Model::Input& cardiacInput,
      const Model::DenseVector& candidateState,
      const Model::State& currentState,
      const Model::State& previousState,
      Real dt,
      Real relativePerturbation)
  {
    using InputType = Model::Input;
    Heart::CCMLC2014::Numerics::DynamicSystem<PassiveLaw, InputType>
        dynamicSystem(cardiacInput);
    const Real nextTime = currentState.t + dt;

    Model::EvalData evaluationData;
    dynamicSystem.buildEvalData(
        candidateState, currentState, previousState, nextTime, dt, evaluationData);

    Model::DenseMatrix analyticalJacobian;
    dynamicSystem.evaluateJacobian(evaluationData, analyticalJacobian, dt);

    Model::DenseMatrix finiteDifferenceJacobian(
        CCMLC2014Vars::NumberOfVariables, CCMLC2014Vars::NumberOfVariables);
    finiteDifferenceJacobian.setZero();

    for (Index j = 0; j < CCMLC2014Vars::NumberOfVariables; ++j)
    {
      const Real perturbation =
        relativePerturbation * std::max<Real>(1.0, std::abs(candidateState[j]));
      auto statePlus = candidateState;
      auto stateMinus = candidateState;
      statePlus[j] += perturbation;
      stateMinus[j] -= perturbation;

      Model::EvalData evaluationDataPlus;
      Model::EvalData evaluationDataMinus;
      dynamicSystem.buildEvalData(
          statePlus, currentState, previousState, nextTime, dt, evaluationDataPlus);
      dynamicSystem.buildEvalData(
          stateMinus, currentState, previousState, nextTime, dt, evaluationDataMinus);

      Model::DenseVector residualPlus;
      Model::DenseVector residualMinus;
      dynamicSystem.evaluateResidual(evaluationDataPlus, residualPlus);
      dynamicSystem.evaluateResidual(evaluationDataMinus, residualMinus);
      finiteDifferenceJacobian.col(j) =
        (residualPlus - residualMinus) / (2.0 * perturbation);
    }

    return (analyticalJacobian - finiteDifferenceJacobian).norm()
      / std::max<Real>(finiteDifferenceJacobian.norm(), 1e-14);
  }
}

/// @brief Verifies initialize uses input active defaults when initial active state is zero for CCMLC 2014 test by checking tolerance-based numerical results.
TEST(CCMLC2014Test, InitializeUsesInputActiveDefaultsWhenInitialActiveStateIsZero)
{
  auto cardiacInput = makeGenericCardiacInput();
  cardiacInput.initFibDef = 0.12;
  cardiacInput.initActiveStiffness = 0.25;
  cardiacInput.initActiveStress = 0.5;

  Model model(cardiacInput);

  Model::State initial;
  initial.t = 0.0;
  model.initialize(initial);

  const auto& state = model.getState();
  EXPECT_NEAR(state.ec, 0.12, 1e-14);
  EXPECT_NEAR(state.gamma, 0.5, 1e-14);
  EXPECT_NEAR(state.beta, 1.0, 1e-14);
  EXPECT_NEAR(state.kc, 0.25, 1e-14);
  EXPECT_NEAR(state.tauc, 0.5, 1e-14);
}

/// @brief Verifies initialize uses provided gamma and beta for CCMLC 2014 test by checking tolerance-based numerical results.
TEST(CCMLC2014Test, InitializeUsesProvidedGammaAndBeta)
{
  auto cardiacInput = makeGenericCardiacInput();

  Model model(cardiacInput);

  Model::State initial;
  initial.t = 0.0;
  initial.ec = 0.03;
  initial.gamma = 0.4;
  initial.beta = 0.2;
  initial.kc = 1.0;
  initial.tauc = 1.0;
  model.initialize(initial);

  const auto& state = model.getState();
  EXPECT_NEAR(state.ec, 0.03, 1e-14);
  EXPECT_NEAR(state.gamma, 0.4, 1e-14);
  EXPECT_NEAR(state.beta, 0.2, 1e-14);
  EXPECT_NEAR(state.kc, 0.16, 1e-14);
  EXPECT_NEAR(state.tauc, 0.08, 1e-14);
}

/// @brief Verifies initialize uses relaxation target when initial W is zero for CCMLC 2014 test by checking tolerance-based numerical results.
TEST(CCMLC2014Test, InitializeUsesRelaxationTargetWhenInitialWIsZero)
{
  auto cardiacInput = makeGenericCardiacInput();
  cardiacInput.m0 = [](Real) { return 1.4; };

  Model model(cardiacInput);

  Model::State initial;
  initial.t = 0.0;
  initial.ec = 0.2;
  initial.w = 0.0;
  model.initialize(initial);

  const auto& state = model.getState();
  EXPECT_NEAR(state.ec, 0.2, 1e-14);
  EXPECT_NEAR(state.w, 1.4, 1e-14);
}

/// @brief Verifies step converges and advances time for CCMLC 2014 test by checking tolerance-based numerical results, true predicates.
TEST(CCMLC2014Test, StepConvergesAndAdvancesTime)
{
  auto cardiacInput = makeGenericCardiacInput();
  cardiacInput.initFibDef = 0.0;
  cardiacInput.initActiveStiffness = 0.0;
  cardiacInput.initActiveStress = 0.0;

  Model model(cardiacInput);
  model.setMaxIterations(200)
       .setAbsoluteTolerance(1e-8)
       .setRelativeTolerance(1e-8)
       .setStepTolerance(1e-10)
       .setDampingFactor(1.0);

  Model::State initial;
  initial.t = 0.0;
  initial.y = 0.0;
  initial.v = 0.0;
  initial.pv = cardiacInput.pAt(0.0) - 100.0;
  initial.par = 11000.0;
  initial.pd = 10000.0;
  model.initialize(initial);

  const Real dt = 1e-3;
  const auto report = model.step(dt);
  EXPECT_TRUE(report.converged);
  EXPECT_NEAR(model.getState().t, dt, 1e-14);
}

/// @brief Verifies example setup completes three cycles with physical state for CCMLC 2014 test by checking tolerance-based numerical results, true predicates.
TEST(CCMLC2014Test, ExampleSetupCompletesThreeCyclesWithPhysicalState)
{
  auto cardiacInput = makeExampleCardiacInput();
  Model model(cardiacInput);
  model.setMaxIterations(200)
    .setAbsoluteTolerance(1e-8)
    .setRelativeTolerance(1e-8)
    .setStepTolerance(1e-10)
    .setDampingFactor(1.0);

  model.initialize(makeExampleInitialState(cardiacInput));

  const Real dt = 1e-3;
  const int nsteps = 3 * static_cast<int>(0.85 / dt);
  Real maxFinalResidual = 0.0;
  size_t maxIterations = 0;

  for (int step = 0; step < nsteps; ++step)
  {
    const auto report = model.step(dt);
    ASSERT_TRUE(report.converged)
      << "Example setup failed at step " << step << ", t = " << model.getState().t
      << ", residual = " << report.finalResidual;
    maxFinalResidual = std::max(maxFinalResidual, report.finalResidual);
    maxIterations = std::max(maxIterations, report.iterations);
  }

  const auto& state = model.getState();
  EXPECT_NEAR(state.t, nsteps * dt, 1e-12);
  EXPECT_TRUE(std::isfinite(state.y));
  EXPECT_TRUE(std::isfinite(state.v));
  EXPECT_TRUE(std::isfinite(state.pv));
  EXPECT_TRUE(std::isfinite(state.par));
  EXPECT_TRUE(std::isfinite(state.pd));
  EXPECT_TRUE(std::isfinite(state.ec));
  EXPECT_TRUE(std::isfinite(state.kc));
  EXPECT_TRUE(std::isfinite(state.tauc));
  EXPECT_TRUE(std::isfinite(state.w));
  EXPECT_GT(state.kc, 0.0);
  EXPECT_GT(state.w, 0.0);
  EXPECT_LT(maxFinalResidual, 5e-2);
  EXPECT_LT(maxIterations, 20u);
}

/// @brief Verifies example setup evolves windkessel and relaxation states for CCMLC 2014 test by checking true predicates.
TEST(CCMLC2014Test, ExampleSetupEvolvesWindkesselAndRelaxationStates)
{
  auto cardiacInput = makeExampleCardiacInput();
  Model model(cardiacInput);
  model.setMaxIterations(200)
    .setAbsoluteTolerance(1e-8)
    .setRelativeTolerance(1e-8)
    .setStepTolerance(1e-10)
    .setDampingFactor(1.0);

  const auto initialState = makeExampleInitialState(cardiacInput);
  model.initialize(initialState);

  const Real dt = 1e-3;
  const int nsteps = static_cast<int>(0.85 / dt);

  Real minPar = initialState.par;
  Real maxPar = initialState.par;
  Real minPd = initialState.pd;
  Real maxPd = initialState.pd;
  Real minW = initialState.w;
  Real maxW = initialState.w;
  Real maxActiveStress = initialState.tauc;

  for (int step = 0; step < nsteps; ++step)
  {
    const auto report = model.step(dt);
    ASSERT_TRUE(report.converged) << "Failed at step " << step;

    const auto& state = model.getState();
    minPar = std::min(minPar, state.par);
    maxPar = std::max(maxPar, state.par);
    minPd = std::min(minPd, state.pd);
    maxPd = std::max(maxPd, state.pd);
    minW = std::min(minW, state.w);
    maxW = std::max(maxW, state.w);
    maxActiveStress = std::max(maxActiveStress, state.tauc);
  }

  EXPECT_GT(maxPar - minPar, 100.0);
  EXPECT_GT(maxPd - minPd, 100.0);
  EXPECT_GT(maxW - minW, 1e-2);
  EXPECT_GT(maxActiveStress, 100.0);
}

/// @brief Verifies non newtonian example trajectory differs from newtonian for CCMLC 2014 test by checking true predicates.
TEST(CCMLC2014Test, NonNewtonianExampleTrajectoryDiffersFromNewtonian)
{
  auto nonNewtonianInput = makeExampleCardiacInput();
  auto newtonianInput = nonNewtonianInput;
  newtonianInput.windkesselRheology =
    Heart::CCMLC2014::Model::WindkesselRheology::Newtonian;

  Model nonNewtonianModel(nonNewtonianInput);
  nonNewtonianModel.setMaxIterations(200)
    .setAbsoluteTolerance(1e-8)
    .setRelativeTolerance(1e-8)
    .setStepTolerance(1e-10)
    .setDampingFactor(1.0);
  nonNewtonianModel.initialize(makeExampleInitialState(nonNewtonianInput));

  Model newtonianModel(newtonianInput);
  newtonianModel.setMaxIterations(200)
    .setAbsoluteTolerance(1e-8)
    .setRelativeTolerance(1e-8)
    .setStepTolerance(1e-10)
    .setDampingFactor(1.0);
  newtonianModel.initialize(makeExampleInitialState(newtonianInput));

  const Real dt = 1e-3;
  const int nsteps = static_cast<int>(0.85 / dt);
  for (int step = 0; step < nsteps; ++step)
  {
    const auto nonNewtonianReport = nonNewtonianModel.step(dt);
    const auto newtonianReport = newtonianModel.step(dt);
    ASSERT_TRUE(nonNewtonianReport.converged) << "Non-Newtonian failed at step " << step;
    ASSERT_TRUE(newtonianReport.converged) << "Newtonian failed at step " << step;
  }

  EXPECT_GT(
    std::abs(nonNewtonianModel.getState().par - newtonianModel.getState().par), 10.0);
  EXPECT_GT(
    std::abs(nonNewtonianModel.getState().pd - newtonianModel.getState().pd), 10.0);
}

/// @brief Verifies load dependent relaxation accelerates negative activation decay for CCMLC 2014 test.
TEST(CCMLC2014Test, LoadDependentRelaxationAcceleratesNegativeActivationDecay)
{
  auto baseInput = makeGenericCardiacInput();
  baseInput.u = [](Real) { return -20.0; };
  baseInput.m0 = [](Real) { return 1.0; };
  baseInput.alphaR = 0.12;

  auto relaxedInput = baseInput;
  relaxedInput.m0 = [](Real) { return 1.6; };

  const Real dt = 1e-3;
  Model::DenseVector baseCandidateState = makeCandidateState();
  baseCandidateState[CCMLC2014Vars::LoadDependentRelaxation] = 1.0;

  Model::DenseVector relaxedCandidateState = baseCandidateState;
  relaxedCandidateState[CCMLC2014Vars::LoadDependentRelaxation] =
    (1.0 + (dt / relaxedInput.alphaR) * 1.6) / (1.0 + dt / relaxedInput.alphaR);

  Model::State currentState = makeCurrentState();
  currentState.w = 1.0;
  Model::State previousState = makePreviousState(currentState);
  previousState.w = 1.0;

  using InputType = Model::Input;
  Heart::CCMLC2014::Numerics::DynamicSystem<PassiveLaw, InputType> baseSystem(baseInput);
  Heart::CCMLC2014::Numerics::DynamicSystem<PassiveLaw, InputType> relaxedSystem(
    relaxedInput);

  Model::EvalData baseData;
  Model::EvalData relaxedData;
  baseSystem.buildEvalData(
    baseCandidateState, currentState, previousState, currentState.t + dt, dt, baseData);
  relaxedSystem.buildEvalData(relaxedCandidateState, currentState, previousState,
    currentState.t + dt, dt, relaxedData);

  Model::DenseVector baseResidual;
  Model::DenseVector relaxedResidual;
  baseSystem.evaluateResidual(baseData, baseResidual);
  relaxedSystem.evaluateResidual(relaxedData, relaxedResidual);

  EXPECT_GT(relaxedData.active.wCurrent, baseData.active.wCurrent);
  EXPECT_GT(relaxedData.active.relaxationDrive, baseData.active.relaxationDrive);
  EXPECT_GT(relaxedResidual[CCMLC2014Vars::ActiveStiffness],
    baseResidual[CCMLC2014Vars::ActiveStiffness]);
}

/// @brief Verifies dynamic jacobian matches finite difference for CCMLC 2014 test.
TEST(CCMLC2014Test, DynamicJacobianMatchesFiniteDifference)
{
  auto cardiacInput = makeGenericCardiacInput();

  Model::DenseVector candidateState = makeCandidateState();
  Model::State currentState = makeCurrentState();
  Model::State previousState = makePreviousState(currentState);

  const Real dt = 1e-3;
  const Real relativeError = computeDynamicJacobianRelativeError(
      cardiacInput, candidateState, currentState, previousState, dt, 1e-7);
  EXPECT_LT(relativeError, 1e-3);
}

/// @brief Recruitment segments and instantaneous relaxation satisfy independent residual and tangent formulas.
TEST(CCMLC2014Test, RecruitmentSegmentsAndInstantaneousRelaxation)
{
  auto input = makeGenericCardiacInput();
  input.alpha = 0;
  input.alphaR = 0;
  // Keep active residual magnitudes comparable to their strain derivatives.
  input.k0 = 1;
  input.sigma0 = 2;
  input.m0 = [](Real ec) { return 1.2 + 0.25 * ec; };
  input.dm0 = [](Real) { return 0.25; };
  Heart::CCMLC2014::Numerics::DynamicSystem<PassiveLaw, Model::Input> system(input);
  constexpr Real dt = 1e-3;
  constexpr Real activation = 25;
  constexpr Real perturbation = 1e-7;
  constexpr std::array<Real, 7> strains{-0.45, -0.05, 0.515, 0.865, 1.15, 1.85, 2.5};
  constexpr std::array<Real, 7> recruitment{0, 0.19, 0.56, 0.87, 1, 0.5, 0};
  constexpr std::array<Real, 7> slopes{
    0, 0.38 / 0.7, 0.36 / 0.43, 0.26 / 0.27, 0, -1 / 1.1, 0};
  for (size_t i = 0; i < strains.size(); ++i)
  {
    SCOPED_TRACE(strains[i]);
    auto candidate = makeCandidateState();
    auto current = makeCurrentState();
    auto previous = makePreviousState(current);
    candidate[CCMLC2014Vars::FiberDeformation] = strains[i];
    current.ec = previous.ec = strains[i];
    Model::EvalData data;
    system.buildEvalData(candidate, current, previous, current.t + dt, dt, data);
    Model::DenseVector residual;
    Model::DenseMatrix jacobian;
    system.evaluateResidual(data, residual);
    system.evaluateJacobian(data, jacobian, dt);
    const Real kc = candidate[CCMLC2014Vars::ActiveStiffness];
    const Real tauc = candidate[CCMLC2014Vars::ActiveStress];
    const Real kcDot = (1.5 * kc - 2 * current.kc + 0.5 * previous.kc) / dt;
    const Real taucDot = (1.5 * tauc - 2 * current.tauc + 0.5 * previous.tauc) / dt;
    EXPECT_NEAR(residual[CCMLC2014Vars::ActiveStiffness],
      kcDot + activation * kc - recruitment[i] * input.k0 * activation, 1e-8);
    EXPECT_NEAR(residual[CCMLC2014Vars::ActiveStress],
      taucDot + activation * tauc - recruitment[i] * input.sigma0 * activation, 1e-8);
    EXPECT_NEAR(residual[CCMLC2014Vars::LoadDependentRelaxation],
      candidate[CCMLC2014Vars::LoadDependentRelaxation] - (1.2 + 0.25 * strains[i]),
      1e-14);
    for (const Index column :
      {CCMLC2014Vars::FiberDeformation, CCMLC2014Vars::LoadDependentRelaxation})
    {
      auto plus = candidate;
      auto minus = candidate;
      plus[column] += perturbation;
      minus[column] -= perturbation;
      Model::EvalData dataPlus, dataMinus;
      Model::DenseVector residualPlus, residualMinus;
      system.buildEvalData(plus, current, previous, current.t + dt, dt, dataPlus);
      system.buildEvalData(minus, current, previous, current.t + dt, dt, dataMinus);
      system.evaluateResidual(dataPlus, residualPlus);
      system.evaluateResidual(dataMinus, residualMinus);
      for (const Index row : {CCMLC2014Vars::ActiveStiffness, CCMLC2014Vars::ActiveStress,
             CCMLC2014Vars::LoadDependentRelaxation})
      {
        const Real derivative =
          (residualPlus[row] - residualMinus[row]) / (2 * perturbation);
        EXPECT_NEAR(jacobian(row, column), derivative,
          1e-6 * std::max<Real>(1, std::abs(derivative)));
      }
    }
    EXPECT_NEAR(jacobian(CCMLC2014Vars::ActiveStiffness, CCMLC2014Vars::FiberDeformation),
      -slopes[i] * input.k0 * activation, 1e-8);
    EXPECT_NEAR(jacobian(CCMLC2014Vars::ActiveStress, CCMLC2014Vars::FiberDeformation),
      -slopes[i] * input.sigma0 * activation - kc * 1.5 / dt, 1e-8);
    EXPECT_EQ(
      jacobian(CCMLC2014Vars::LoadDependentRelaxation, CCMLC2014Vars::FiberDeformation),
      -0.25);
    EXPECT_EQ(jacobian(CCMLC2014Vars::LoadDependentRelaxation,
                CCMLC2014Vars::LoadDependentRelaxation),
      1);
  }
}

/// @brief Verifies windkessel residual and jacobian use branch flows and BDF 2 for CCMLC 2014 test by checking tolerance-based numerical results.
TEST(CCMLC2014Test, WindkesselResidualAndJacobianUseBranchFlowsAndBDF2)
{
  auto cardiacInput = makeGenericCardiacInput();

  Model::DenseVector candidateState = makeCandidateState();
  candidateState[CCMLC2014Vars::VentricularPressure] = 1.3e4;
  candidateState[CCMLC2014Vars::ArterialPressure] = 9.0e3;
  candidateState[CCMLC2014Vars::DistalPressure] = 8.0e3;

  Model::State currentState = makeCurrentState();
  currentState.par = 8.8e3;
  currentState.pd = 7.9e3;

  Model::State previousState = makePreviousState(currentState);
  previousState.par = 8.7e3;
  previousState.pd = 7.8e3;

  using InputType = Model::Input;
  Heart::CCMLC2014::Numerics::DynamicSystem<PassiveLaw, InputType> dynamicSystem(
    cardiacInput);

  const Real dt = 1e-3;
  const Real nextTime = currentState.t + dt;

  Model::EvalData evaluationData;
  dynamicSystem.buildEvalData(
    candidateState, currentState, previousState, nextTime, dt, evaluationData);

  Model::DenseVector residual;
  dynamicSystem.evaluateResidual(evaluationData, residual);

  Model::DenseMatrix jacobian;
  dynamicSystem.evaluateJacobian(evaluationData, jacobian, dt);

  const Real a0 = 1.5 / dt;
  const Real parDot = a0 * candidateState[CCMLC2014Vars::ArterialPressure] -
    2.0 / dt * currentState.par + 0.5 / dt * previousState.par;
  const Real pdDot = a0 * candidateState[CCMLC2014Vars::DistalPressure] -
    2.0 / dt * currentState.pd + 0.5 / dt * previousState.pd;

  const Real expectedOutflow = cardiacInput.Kar *
    (candidateState[CCMLC2014Vars::VentricularPressure] -
      candidateState[CCMLC2014Vars::ArterialPressure]);
  const Real expectedProximalFlow = (candidateState[CCMLC2014Vars::ArterialPressure] -
                                      candidateState[CCMLC2014Vars::DistalPressure]) /
      cardiacInput.Rp -
    expectedOutflow;
  const Real expectedDistalFlow = (candidateState[CCMLC2014Vars::DistalPressure] -
                                    candidateState[CCMLC2014Vars::ArterialPressure]) /
      cardiacInput.Rp -
    (cardiacInput.pSv(nextTime) - candidateState[CCMLC2014Vars::DistalPressure]) /
      cardiacInput.Rd;

  EXPECT_NEAR(residual[CCMLC2014Vars::ArterialPressure],
    cardiacInput.Cp * parDot + expectedProximalFlow, 1e-13);
  EXPECT_NEAR(residual[CCMLC2014Vars::DistalPressure],
    cardiacInput.Cd * pdDot + expectedDistalFlow, 1e-13);

  EXPECT_NEAR(
    jacobian(CCMLC2014Vars::ArterialPressure, CCMLC2014Vars::VentricularPressure),
    -cardiacInput.Kar, 1e-14);
  EXPECT_NEAR(jacobian(CCMLC2014Vars::ArterialPressure, CCMLC2014Vars::ArterialPressure),
    cardiacInput.Cp * a0 + 1.0 / cardiacInput.Rp + cardiacInput.Kar, 1e-14);
  EXPECT_NEAR(jacobian(CCMLC2014Vars::ArterialPressure, CCMLC2014Vars::DistalPressure),
    -1.0 / cardiacInput.Rp, 1e-14);
  EXPECT_NEAR(jacobian(CCMLC2014Vars::DistalPressure, CCMLC2014Vars::ArterialPressure),
    -1.0 / cardiacInput.Rp, 1e-14);
  EXPECT_NEAR(jacobian(CCMLC2014Vars::DistalPressure, CCMLC2014Vars::DistalPressure),
    cardiacInput.Cd * a0 + 1.0 / cardiacInput.Rp + 1.0 / cardiacInput.Rd, 1e-14);
}

/// @brief Verifies non newtonian windkessel path produces finite consistent jacobian for CCMLC 2014 test by checking true predicates.
TEST(CCMLC2014Test, NonNewtonianWindkesselPathProducesFiniteConsistentJacobian)
{
  auto cardiacInput = makeGenericCardiacInput();
  cardiacInput.windkesselRheology =
    Heart::CCMLC2014::Model::WindkesselRheology::CarreauYasuda;
  cardiacInput.proximalRadius = 0.015;
  cardiacInput.proximalLength = 0.4;
  cardiacInput.distalRadius = 0.0007;
  cardiacInput.distalLength = 0.004;
  cardiacInput.mu_0 = 0.04868;
  cardiacInput.mu_Inf = 0.003605;
  cardiacInput.lambda = 3.39;
  cardiacInput.n = 0.198;
  cardiacInput.yasuda = 1.235;

  Model::DenseVector candidateState = makeCandidateState();
  candidateState[CCMLC2014Vars::VentricularPressure] = 1.3e4;
  candidateState[CCMLC2014Vars::ArterialPressure] = 9.0e3;
  candidateState[CCMLC2014Vars::DistalPressure] = 8.0e3;

  Model::State currentState = makeCurrentState();
  Model::State previousState = makePreviousState(currentState);

  const Real relativeError = computeDynamicJacobianRelativeError(
    cardiacInput, candidateState, currentState, previousState, 1e-3, 1e-7);
  EXPECT_TRUE(std::isfinite(relativeError));
  EXPECT_LT(relativeError, 3e-3);
}

/// @brief Verifies dynamic jacobian matches finite difference across perturbation scales for CCMLC 2014 test.
TEST(CCMLC2014Test, DynamicJacobianMatchesFiniteDifferenceAcrossPerturbationScales)
{
  auto cardiacInput = makeGenericCardiacInput();

  Model::DenseVector candidateState = makeCandidateState();
  Model::State currentState = makeCurrentState();
  Model::State previousState = makePreviousState(currentState);

  const Real dt = 1e-3;
  const std::array<Real, 3> perturbations{{1e-6, 1e-7, 1e-8}};
  for (const Real relativePerturbation : perturbations)
  {
    const Real relativeError = computeDynamicJacobianRelativeError(
        cardiacInput,
        candidateState,
        currentState,
        previousState,
        dt,
        relativePerturbation);
    EXPECT_LT(relativeError, 2e-3);
  }
}

/// @brief Verifies dynamic jacobian matches finite difference across data scales for CCMLC 2014 test.
TEST(CCMLC2014Test, DynamicJacobianMatchesFiniteDifferenceAcrossDataScales)
{
  auto cardiacInput = makeGenericCardiacInput();

  const std::array<Real, 3> pressureScales{{0.2, 1.0, 5.0}};
  for (const Real pressureScale : pressureScales)
  {
    Model::DenseVector candidateState = makeCandidateState(pressureScale);
    Model::State currentState = makeCurrentState(pressureScale);
    Model::State previousState = makePreviousState(currentState, pressureScale);

    const Real relativeError = computeDynamicJacobianRelativeError(
        cardiacInput,
        candidateState,
        currentState,
        previousState,
        1e-3,
        1e-7);
    EXPECT_LT(relativeError, 3e-3);
  }
}

namespace
{
  template <class Rheology>
  void checkBranchFlowDerivative(const Model::Input& input, Real dp)
  {
    constexpr Real length = 0.1;
    constexpr Real radius = 0.004;
    constexpr Real resistance = 1e8;
    constexpr Real relativeStep = 1e-5;
    constexpr Real relativeTolerance = 1e-5;
    const auto [flow, derivative] =
      Rheology::flowLaw(dp, length, radius, resistance, input);
    const auto [reverseFlow, reverseDerivative] =
      Rheology::flowLaw(-dp, length, radius, resistance, input);
    ASSERT_GT(flow, 0);
    ASSERT_GT(derivative, 0);
    EXPECT_DOUBLE_EQ(reverseFlow, -flow);
    EXPECT_DOUBLE_EQ(reverseDerivative, derivative);
    const Real step = relativeStep * dp;
    const Real finiteDifference =
      (Rheology::flowLaw(dp + step, length, radius, resistance, input).first -
        Rheology::flowLaw(dp - step, length, radius, resistance, input).first) /
      (2 * step);
    EXPECT_NEAR(derivative / finiteDifference, 1, relativeTolerance);
    for (const auto& [invalidLength, invalidRadius] :
      {std::pair<Real, Real>{0, radius}, {length, 0}, {-length, radius}})
    {
      const auto fallback =
        Rheology::flowLaw(dp, invalidLength, invalidRadius, resistance, input);
      EXPECT_DOUBLE_EQ(fallback.first, dp / resistance);
      EXPECT_DOUBLE_EQ(fallback.second, 1 / resistance);
    }
  }
}

/// @brief Checks power-law flow against Poiseuille's law and a nonlinear derivative.
TEST(CCMLC2014Test, PowerLawBranchFlowAndNewtonianLimit)
{
  using Law = Heart::CCMLC2014::Physics::Rheology::PowerLaw;
  Model::Input input;
  input.m = 0.004;
  input.n = 1;
  constexpr Real length = 0.1;
  constexpr Real radius = 0.004;
  constexpr Real pressure = 1000;
  constexpr Real fallbackResistance = 1e8;
  const Real conductance =
    std::numbers::pi_v<Real> * std::pow(radius, 4) / (8 * input.m * length);
  const auto result = Law::flowLaw(pressure, length, radius, fallbackResistance, input);
  EXPECT_NEAR(result.first / (conductance * pressure), 1, 1e-12);
  EXPECT_NEAR(result.second / conductance, 1, 1e-12);
  input.n = 0.75;
  checkBranchFlowDerivative<Law>(input, pressure);
  const auto zero = Law::flowLaw(0, length, radius, fallbackResistance, input);
  EXPECT_EQ(zero.first, 0);
  EXPECT_EQ(zero.second, 1 / fallbackResistance);
}

/// @brief Checks Cross flow symmetry, derivative and its constant-viscosity limit.
TEST(CCMLC2014Test, CrossBranchFlowAndNewtonianLimit)
{
  using Law = Heart::CCMLC2014::Physics::Rheology::Cross;
  Model::Input input;
  input.mu_0 = 0.02;
  input.mu_Inf = 0.004;
  input.lambda = 1;
  input.n = 0.6;
  checkBranchFlowDerivative<Law>(input, 1000);
  input.mu_0 = input.mu_Inf;
  constexpr Real length = 0.1;
  constexpr Real radius = 0.004;
  const Real conductance =
    std::numbers::pi_v<Real> * std::pow(radius, 4) / (8 * input.mu_0 * length);
  const auto result = Law::flowLaw(1000, length, radius, 1e8, input);
  EXPECT_NEAR(result.first / (1000 * conductance), 1, 1e-7);
  EXPECT_NEAR(result.second / conductance, 1, 1e-7);
  const auto zero = Law::flowLaw(0, length, radius, 1e8, input);
  EXPECT_EQ(zero.first, 0);
  EXPECT_NEAR(zero.second / conductance, 1, 1e-12);
}

/// @brief Checks both Quemada pressure regimes and the near-zero linear fallback.
TEST(CCMLC2014Test, QuemadaBranchFlowDerivativeAcrossPressureRegimes)
{
  using Law = Heart::CCMLC2014::Physics::Rheology::Quemada;
  Model::Input input;
  input.k_0 = 3.5;
  input.k_Inf = 2;
  input.phi_quemada = 0.45;
  input.gamma_c = 10;
  input.mu_plasma = 0.0012;
  input.mu_Inf = 0.004;
  checkBranchFlowDerivative<Law>(input, 0.1);
  checkBranchFlowDerivative<Law>(input, 1000);
  const auto zero = Law::flowLaw(0, 0.1, 0.004, 1e8, input);
  EXPECT_EQ(zero.first, 0);
  EXPECT_EQ(zero.second, 1e-8);
}

namespace
{
  /// @brief Minimal evaluation record consumed by the standalone physics evaluators.
  struct PhysicsEvaluation
  {
      Real y, yPrev, yMid, pv, pvPrev, pvMid, par, parPrev, parMid, pd, pdPrev, pdMid;
      Real dt, vel, sqrtC, C, strain1D, diffGreen;
      Real stressPassive, diffStressPassive, stressViscous, diffStressViscous;
      Real pAtCur, pAtPrev, cavityFluxCur, cavityFluxPrev;
      Real dCavityFluxCur_dPv, dCavityFluxCur_dPar;
  };
}

/// @brief Midpoint kinematics satisfy displacement, pressure and strain identities.
TEST(CCMLC2014Test, StandaloneMidpointKinematics)
{
  auto input = makeGenericCardiacInput();
  input.R0 = 2;
  PhysicsEvaluation data{};
  data.y = 0.4;
  data.yPrev = 0.1;
  data.pv = 6;
  data.pvPrev = 2;
  data.par = 8;
  data.parPrev = 4;
  data.pd = 10;
  data.pdPrev = 2;
  data.dt = 0.3;
  Heart::CCMLC2014::Physics::MidpointKinematics model(input);
  model.evaluate(data);
  EXPECT_DOUBLE_EQ(data.yMid, 0.25);
  EXPECT_DOUBLE_EQ(data.pvMid, 4);
  EXPECT_DOUBLE_EQ(data.parMid, 6);
  EXPECT_DOUBLE_EQ(data.pdMid, 6);
  EXPECT_NEAR(data.vel, 1, 1e-14);
  EXPECT_DOUBLE_EQ(data.sqrtC, 1.125);
  EXPECT_DOUBLE_EQ(data.C, 1.265625);
  EXPECT_DOUBLE_EQ(data.strain1D, 0.1328125);
  EXPECT_DOUBLE_EQ(data.diffGreen, 0.5625);
}

/// @brief Passive stress's stored tangent differentiates with respect to the new displacement.
TEST(CCMLC2014Test, StandalonePassiveStressChainRule)
{
  auto input = makeGenericCardiacInput();
  input.R0 = 2;
  PhysicsEvaluation data{};
  data.y = 0.4;
  data.yPrev = 0.1;
  data.dt = 0.3;
  Heart::CCMLC2014::Physics::MidpointKinematics kinematics(input);
  Heart::CCMLC2014::Physics::PassiveStressEvaluator<PassiveLaw, Model::Input> stress(
    input);
  auto evaluate = [&](PhysicsEvaluation& state) {
    kinematics.evaluate(state);
    stress.evaluate(state);
  };
  evaluate(data);
  constexpr Real Step = 1e-6;
  auto plus = data, minus = data;
  plus.y += Step;
  minus.y -= Step;
  evaluate(plus);
  evaluate(minus);
  const Real derivative = (plus.stressPassive - minus.stressPassive) / (2 * Step);
  EXPECT_GT(std::abs(data.stressPassive), 1);
  EXPECT_GT(std::abs(derivative), 1);
  EXPECT_NEAR(data.diffStressPassive / derivative, 1, 1e-7);
}

/// @brief Viscous stress and its midpoint derivative agree with an independent scalar law.
TEST(CCMLC2014Test, StandaloneViscousStressMidpointDerivative)
{
  auto input = makeGenericCardiacInput();
  input.R0 = 2;
  input.eta = 0.7;
  PhysicsEvaluation data{};
  data.y = 0.4;
  data.yPrev = 0.1;
  data.dt = 0.3;
  Heart::CCMLC2014::Physics::MidpointKinematics kinematics(input);
  Heart::CCMLC2014::Physics::ViscousStressEvaluator stress(input);
  kinematics.evaluate(data);
  stress.evaluate(data);
  const Real stretch = 1.125;
  const Real exact = 2 * input.eta / input.R0 * (stretch + 2 * std::pow(stretch, -11));
  EXPECT_NEAR(data.stressViscous, exact, 1e-12);
  constexpr Real Step = 1e-6;
  auto value = [&](Real midpoint) {
    const Real s = 1 + midpoint / input.R0;
    const Real velocity = 2 * (midpoint - data.yPrev) / data.dt;
    return 2 * input.eta / input.R0 * (s + 2 * std::pow(s, -11)) * velocity;
  };
  const Real derivative =
    (value(data.yMid + Step) - value(data.yMid - Step)) / (2 * Step);
  EXPECT_NEAR(data.diffStressViscous, derivative, 1e-8);
}

/// @brief Current and previous valve branches are independent and retain exact pressure slopes.
TEST(CCMLC2014Test, StandaloneValveBranchTransitions)
{
  auto input = makeGenericCardiacInput();
  input.Kat = 3;
  input.Kp = 0.2;
  input.Kar = 1.7;
  Heart::CCMLC2014::Physics::CavityFluxEvaluator model(input);
  constexpr std::array<Real, 5> Pressure{1, 2, 3, 5, 6};
  constexpr std::array<Real, 5> Flux{-3, 0, 0.2, 0.6, 2.3};
  constexpr std::array<Real, 5> VentricularSlope{3, 3, 0.2, 0.2, 1.7};
  constexpr std::array<Real, 5> ArterialSlope{0, 0, 0, 0, -1.5};
  for (size_t current = 0; current < Pressure.size(); ++current)
  {
    for (size_t previous = 0; previous < Pressure.size(); ++previous)
    {
      PhysicsEvaluation data{};
      data.pAtCur = data.pAtPrev = 2;
      data.par = data.parPrev = 5;
      data.pv = Pressure[current];
      data.pvPrev = Pressure[previous];
      model.evaluate(data);
      EXPECT_NEAR(data.cavityFluxCur, Flux[current], 1e-14);
      EXPECT_NEAR(data.cavityFluxPrev, Flux[previous], 1e-14);
      EXPECT_DOUBLE_EQ(data.dCavityFluxCur_dPv, VentricularSlope[current]);
      EXPECT_DOUBLE_EQ(data.dCavityFluxCur_dPar, ArterialSlope[current]);
      if (current == 1 || current == 3)
        continue;
      constexpr Real Step = 1e-5;
      auto plus = data, minus = data;
      plus.pv += Step;
      minus.pv -= Step;
      model.evaluate(plus);
      model.evaluate(minus);
      EXPECT_NEAR((plus.cavityFluxCur - minus.cavityFluxCur) / (2 * Step),
        data.dCavityFluxCur_dPv, 1e-9);
      plus = minus = data;
      plus.par += Step;
      minus.par -= Step;
      model.evaluate(plus);
      model.evaluate(minus);
      EXPECT_NEAR((plus.cavityFluxCur - minus.cavityFluxCur) / (2 * Step),
        data.dCavityFluxCur_dPar, 1e-9);
    }
  }
}
