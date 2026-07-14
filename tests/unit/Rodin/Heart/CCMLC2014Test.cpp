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

  Real examplePeriodicActivation(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);

    if (tau < 0.13)  return 0.0;
    if (tau < 0.141) return 35.0 * ((tau - 0.13) / 0.011);
    if (tau < 0.281) return 35.0;
    if (tau < 0.361) return 35.0 - 55.0 * ((tau - 0.281) / 0.08);
    if (tau < 0.45)  return -20.0;
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
    cardiacInput.Kp  = 5.0e-10;
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
    initialState.gamma =
      std::sqrt(std::max<Real>(cardiacInput.initActiveStiffness, 0.0));
    initialState.beta =
      (initialState.gamma > 0.0)
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
      << "Example setup failed at step " << step
      << ", t = " << model.getState().t
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
    ASSERT_TRUE(nonNewtonianReport.converged)
      << "Non-Newtonian failed at step " << step;
    ASSERT_TRUE(newtonianReport.converged)
      << "Newtonian failed at step " << step;
  }

  EXPECT_GT(
      std::abs(nonNewtonianModel.getState().par - newtonianModel.getState().par),
      10.0);
  EXPECT_GT(
      std::abs(nonNewtonianModel.getState().pd - newtonianModel.getState().pd),
      10.0);
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
  Heart::CCMLC2014::Numerics::DynamicSystem<PassiveLaw, InputType>
      dynamicSystem(cardiacInput);

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
  const Real parDot =
      a0 * candidateState[CCMLC2014Vars::ArterialPressure]
    - 2.0 / dt * currentState.par
    + 0.5 / dt * previousState.par;
  const Real pdDot =
      a0 * candidateState[CCMLC2014Vars::DistalPressure]
    - 2.0 / dt * currentState.pd
    + 0.5 / dt * previousState.pd;

  const Real expectedOutflow =
    cardiacInput.Kar
    * (candidateState[CCMLC2014Vars::VentricularPressure]
       - candidateState[CCMLC2014Vars::ArterialPressure]);
  const Real expectedProximalFlow =
    (candidateState[CCMLC2014Vars::ArterialPressure]
     - candidateState[CCMLC2014Vars::DistalPressure]) / cardiacInput.Rp
    - expectedOutflow;
  const Real expectedDistalFlow =
    (candidateState[CCMLC2014Vars::DistalPressure]
     - candidateState[CCMLC2014Vars::ArterialPressure]) / cardiacInput.Rp
    - (cardiacInput.pSv(nextTime)
       - candidateState[CCMLC2014Vars::DistalPressure]) / cardiacInput.Rd;

  EXPECT_NEAR(
      residual[CCMLC2014Vars::ArterialPressure],
      cardiacInput.Cp * parDot + expectedProximalFlow,
      1e-13);
  EXPECT_NEAR(
      residual[CCMLC2014Vars::DistalPressure],
      cardiacInput.Cd * pdDot + expectedDistalFlow,
      1e-13);

  EXPECT_NEAR(
      jacobian(
        CCMLC2014Vars::ArterialPressure,
        CCMLC2014Vars::VentricularPressure),
      -cardiacInput.Kar,
      1e-14);
  EXPECT_NEAR(
      jacobian(
        CCMLC2014Vars::ArterialPressure,
        CCMLC2014Vars::ArterialPressure),
      cardiacInput.Cp * a0 + 1.0 / cardiacInput.Rp + cardiacInput.Kar,
      1e-14);
  EXPECT_NEAR(
      jacobian(
        CCMLC2014Vars::ArterialPressure,
        CCMLC2014Vars::DistalPressure),
      -1.0 / cardiacInput.Rp,
      1e-14);
  EXPECT_NEAR(
      jacobian(
        CCMLC2014Vars::DistalPressure,
        CCMLC2014Vars::ArterialPressure),
      -1.0 / cardiacInput.Rp,
      1e-14);
  EXPECT_NEAR(
      jacobian(
        CCMLC2014Vars::DistalPressure,
        CCMLC2014Vars::DistalPressure),
      cardiacInput.Cd * a0 + 1.0 / cardiacInput.Rp + 1.0 / cardiacInput.Rd,
      1e-14);
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
      cardiacInput,
      candidateState,
      currentState,
      previousState,
      1e-3,
      1e-7);
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
