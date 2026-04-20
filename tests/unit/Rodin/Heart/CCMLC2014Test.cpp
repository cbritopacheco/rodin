#include <gtest/gtest.h>

#include <algorithm>
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

TEST(CCMLC2014Test, DynamicJacobianMatchesFiniteDifference)
{
  auto cardiacInput = makeGenericCardiacInput();

  using InputType = Model::Input;
  Heart::CCMLC2014::Numerics::DynamicSystem<PassiveLaw, InputType>
      dynamicSystem(cardiacInput);

  Model::DenseVector candidateState(CCMLC2014Vars::NumberOfVariables);
  candidateState[CCMLC2014Vars::RadialDisplacement] = 8e-5;
  candidateState[CCMLC2014Vars::VentricularPressure] = 1.0e4;
  candidateState[CCMLC2014Vars::ArterialPressure] = 9.0e3;
  candidateState[CCMLC2014Vars::DistalPressure] = 8.0e3;

  Model::State currentState;
  currentState.t = 0.1;
  currentState.y = 7e-5;
  currentState.pv = 9.8e3;
  currentState.par = 8.8e3;
  currentState.pd = 7.9e3;
  currentState.ec = 0.01;
  currentState.gamma = 0.2;
  currentState.beta = 0.3;
  currentState.kc = currentState.gamma * currentState.gamma;
  currentState.tauc = currentState.gamma * currentState.beta;

  Model::State previousState = currentState;
  previousState.t = 0.099;
  previousState.y = 6e-5;

  const Real dt = 1e-3;
  const Real nextTime = currentState.t + dt;

  Model::EvalData evaluationData;
  dynamicSystem.buildEvalData(
      candidateState, currentState, previousState, nextTime, dt, evaluationData);

  Model::DenseMatrix analyticalJacobian;
  dynamicSystem.evaluateJacobian(evaluationData, analyticalJacobian, dt);

  Model::DenseMatrix finiteDifferenceJacobian(
      CCMLC2014Vars::NumberOfVariables, CCMLC2014Vars::NumberOfVariables);
  finiteDifferenceJacobian.setZero();

  const Real relativePerturbation = 1e-7;
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

  const Real relativeError =
    (analyticalJacobian - finiteDifferenceJacobian).norm()
    / std::max<Real>(finiteDifferenceJacobian.norm(), 1e-14);
  EXPECT_LT(relativeError, 1e-3);
}
