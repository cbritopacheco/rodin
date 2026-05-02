#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "Rodin/Heart/CCMLC2014.h"
#include "Rodin/Heart/CCMLC2014/Numerics/DynamicSystem.h"

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
  Heart::CCMLC2014::Numerics::DynamicSystem<PassiveLaw, InputType>
      baseSystem(baseInput);
  Heart::CCMLC2014::Numerics::DynamicSystem<PassiveLaw, InputType>
      relaxedSystem(relaxedInput);

  Model::EvalData baseData;
  Model::EvalData relaxedData;
  baseSystem.buildEvalData(
      baseCandidateState, currentState, previousState, currentState.t + dt, dt, baseData);
  relaxedSystem.buildEvalData(
      relaxedCandidateState, currentState, previousState, currentState.t + dt, dt, relaxedData);

  Model::DenseVector baseResidual;
  Model::DenseVector relaxedResidual;
  baseSystem.evaluateResidual(baseData, baseResidual);
  relaxedSystem.evaluateResidual(relaxedData, relaxedResidual);

  EXPECT_GT(relaxedData.active.wCurrent, baseData.active.wCurrent);
  EXPECT_GT(relaxedData.active.relaxationDrive, baseData.active.relaxationDrive);
  EXPECT_GT(
      relaxedResidual[CCMLC2014Vars::ActiveStiffness],
      baseResidual[CCMLC2014Vars::ActiveStiffness]);
}

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
