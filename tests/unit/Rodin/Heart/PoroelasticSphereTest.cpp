#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "Rodin/Heart/CCMLC2014.h"
#include "Rodin/Heart/PoroelasticSphere.h"
#include "Rodin/Heart/PoroelasticSphere/Numerics/DynamicSystem.h"

using namespace Rodin;

namespace
{
  using Model = Heart::PoroelasticSphereT<>;
  using PassiveLaw = Heart::PoroelasticSpherePassiveLaw<Real>;
  namespace Vars = Heart::PoroelasticSphere::Model;

  Model::Input makeGenericCardiacInput()
  {
    Model::Input cardiacInput;

    cardiacInput.R0 = 2.36e-2;
    cardiacInput.d0 = 1.42e-2;
    cardiacInput.phi0 = 0.1;

    cardiacInput.Es = 3.0e5;
    cardiacInput.mu = 70.0;
    cardiacInput.eta = 70.0;
    cardiacInput.alpha = 3.0;
    cardiacInput.k0 = 1.0e5;
    cardiacInput.sigma0 = 5.0e5;

    cardiacInput.KPhi = 2.0e5;
    cardiacInput.gammaAr = 7.0e-10;
    cardiacInput.gammaVen = 7.0e-10;

    cardiacInput.Rp = 8.0e6;
    cardiacInput.Cp = 5.0e-9;
    cardiacInput.Rd = 1.0e8;
    cardiacInput.Cd = 1.0e-8;

    cardiacInput.Kat = 8.0e-7;
    cardiacInput.Kp  = 5.0e-10;
    cardiacInput.Kar = 1.3e-5;

    cardiacInput.cavityCapacity = 5.0e-12;
    cardiacInput.absRegularization = 1e-14;
    cardiacInput.wallQuadraturePoints = 8;

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
    Model::DenseVector candidateState(Vars::NumberOfVariables);
    candidateState.setZero();
    candidateState[Vars::RadialDisplacement] = scale * 3e-3;
    candidateState[Vars::Porosity] = 0.12;
    candidateState[Vars::VentricularPressure] = scale * 1.0e4;
    candidateState[Vars::ArterialPressure] = scale * 9.0e3;
    candidateState[Vars::DistalPressure] = scale * 8.0e3;
    candidateState[Vars::FiberDeformation] = scale * 0.012;
    candidateState[Vars::ActiveStiffness] = scale * 0.045;
    candidateState[Vars::ActiveStress] = scale * 0.065;
    candidateState[Vars::LoadDependentRelaxation] = 1.2;
    return candidateState;
  }

  Model::State makeCurrentState(Real scale = 1.0)
  {
    Model::State currentState;
    currentState.t = 0.1;
    currentState.y = scale * 2.8e-3;
    currentState.phi = 0.118;
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
    previousState.y = scale * 2.5e-3;
    previousState.phi = 0.117;
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

    cardiacInput.R0 = 2.4e-2;
    cardiacInput.d0 = 1.45e-2;
    cardiacInput.phi0 = 0.1;

    cardiacInput.Es = 3.0e6;
    cardiacInput.mu = 70.0;
    cardiacInput.eta = 70.0;
    cardiacInput.alpha = 1.5;
    cardiacInput.alphaR = 0.12;
    cardiacInput.k0 = 1.0e5;
    cardiacInput.sigma0 = 1.25e5;

    cardiacInput.KPhi = 2.0e5;
    cardiacInput.gammaAr = 7.0e-10;
    cardiacInput.gammaVen = 7.0e-10;

    cardiacInput.Rp = 5.0e7;
    cardiacInput.Cp = 6e-9;
    cardiacInput.Rd = 1.0e8;
    cardiacInput.Cd = 1.0e-9;

    cardiacInput.mu_0 = 5.35;
    cardiacInput.mu_Inf = 0.0033;
    cardiacInput.lambda = 14.445;
    cardiacInput.n = 0.8;
    cardiacInput.m = 0.003;
    cardiacInput.yasuda = 0.62;
    cardiacInput.proximalRadius = 0.0125;
    cardiacInput.proximalLength = 0.35;
    cardiacInput.distalRadius = 0.002;
    cardiacInput.distalLength = 0.55;
    cardiacInput.windkesselRheology = Vars::WindkesselRheology::Cross;

    cardiacInput.Kat = 6.0e-7;
    cardiacInput.Kp = 5.0e-11;
    cardiacInput.Kar = 1.0e-7;

    cardiacInput.cavityCapacity = 5.0e-12;
    cardiacInput.absRegularization = 1e-14;
    cardiacInput.wallQuadraturePoints = 8;

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
    initialState.phi = cardiacInput.phi0;
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
   * @brief Evaluates the model intermediate data at a candidate state with
   * zero rates and a given porosity (steady wall, backward Euler).
   */
  Model::EvalData evaluateSteadyWall(
      const Model::Input& cardiacInput, Real y, Real phi, Real ec)
  {
    Heart::PoroelasticSphere::Numerics::DynamicSystem<PassiveLaw, Model::Input>
        dynamicSystem(cardiacInput);

    Model::State steady;
    steady.t = 0.0;
    steady.y = y;
    steady.phi = phi;
    steady.ec = ec;

    Model::DenseVector x(Vars::NumberOfVariables);
    x.setZero();
    x[Vars::RadialDisplacement] = y;
    x[Vars::Porosity] = phi;
    x[Vars::FiberDeformation] = ec;

    Model::EvalData data;
    dynamicSystem.buildEvalData(x, steady, steady, 1e-3, 1e-3, data);
    return data;
  }

  /**
   * @brief Computes Jacobian agreement error against a central finite-difference approximation.
   *
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
    Heart::PoroelasticSphere::Numerics::DynamicSystem<PassiveLaw, InputType>
        dynamicSystem(cardiacInput);
    const Real nextTime = currentState.t + dt;

    Model::EvalData evaluationData;
    dynamicSystem.buildEvalData(
        candidateState, currentState, previousState, nextTime, dt, evaluationData);

    Model::DenseMatrix analyticalJacobian;
    dynamicSystem.evaluateJacobian(evaluationData, analyticalJacobian, dt);

    Model::DenseMatrix finiteDifferenceJacobian(
        Vars::NumberOfVariables, Vars::NumberOfVariables);
    finiteDifferenceJacobian.setZero();

    for (Index j = 0; j < Vars::NumberOfVariables; ++j)
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

/// @brief Verifies the passive Piola tangent matches finite differences of the components.
TEST(PoroelasticSphereTest, PassivePiolaTangentMatchesFiniteDifference)
{
  const auto cardiacInput = makeGenericCardiacInput();
  PassiveLaw law;
  const Real lr = 0.83;
  const Real lt = 1.17;
  const Real h = 1e-6;

  const auto P = law(cardiacInput.passiveEnergy, lr, lt);
  const auto Prp = law(cardiacInput.passiveEnergy, lr + h, lt);
  const auto Prm = law(cardiacInput.passiveEnergy, lr - h, lt);
  const auto Ptp = law(cardiacInput.passiveEnergy, lr, lt + h);
  const auto Ptm = law(cardiacInput.passiveEnergy, lr, lt - h);

  const Real scale = std::max<Real>(1.0, std::abs(P.dPtt_dlt));
  EXPECT_NEAR(P.dPRR_dlr, (Prp.PRR - Prm.PRR) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.dPRR_dlt, (Ptp.PRR - Ptm.PRR) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.dPtt_dlr, (Prp.Ptt - Prm.Ptt) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.dPtt_dlt, (Ptp.Ptt - Ptm.Ptt) / (2 * h), 1e-6 * scale);
}

/// @brief Verifies the thick incompressible neo-Hookean sphere matches the closed-form inflation law.
TEST(PoroelasticSphereTest, IncompressibleNeoHookeanMatchesClosedForm)
{
  // W = (mu / 2) J_1 with J = 1 is the incompressible neo-Hookean law; the
  // exact inflation pressure is the Green-Ogden formula
  // P_v = (mu / 2) [ (4 / l_out + 1 / l_out^4) - (4 / l_in + 1 / l_in^4) ].
  auto cardiacInput = makeGenericCardiacInput();
  cardiacInput.R0 = 1.0;
  cardiacInput.d0 = 0.5;
  cardiacInput.phi0 = 0.1;
  cardiacInput.Es = 0.0;
  cardiacInput.eta = 0.0;
  cardiacInput.wallQuadraturePoints = 12;
  const Real muNH = 2.0e3;
  using PassiveEnergyType = std::decay_t<decltype(cardiacInput.passiveEnergy)>;
  typename PassiveEnergyType::Parameters p;
  p.mu1 = 0.5 * muNH;
  cardiacInput.passiveEnergy = PassiveEnergyType(p);

  const Real rin = 1.3;
  const auto data = evaluateSteadyWall(cardiacInput, rin - 1.0, cardiacInput.phi0, 0.0);

  const Real Rout = 1.5;
  const Real lin = rin;
  const Real lout = std::cbrt(rin * rin * rin + Rout * Rout * Rout - 1.0) / Rout;
  const Real expected = 0.5 * muNH *
    ((4.0 / lout + 1.0 / std::pow(lout, 4)) - (4.0 / lin + 1.0 / std::pow(lin, 4)));

  EXPECT_NEAR(data.J, 1.0, 1e-15);
  EXPECT_NEAR(data.wallPressure, expected, 1e-9 * std::abs(expected));
}

/// @brief Verifies the membrane limit at J = 1 recovers the CCMLC2014 wall law (passive, viscous, active).
TEST(PoroelasticSphereTest, MembraneLimitRecoversCCMLC2014Law)
{
  auto cardiacInput = makeGenericCardiacInput();
  cardiacInput.d0 = 1e-4 * cardiacInput.R0;
  cardiacInput.wallQuadraturePoints = 4;

  const Real dt = 1e-3;
  const Real y = 4e-3;
  const Real yPrev = 3.9e-3;
  const Real ec = 0.02;

  Heart::PoroelasticSphere::Numerics::DynamicSystem<PassiveLaw, Model::Input>
    dynamicSystem(cardiacInput);

  // Backward Euler start (snm1 == sn), so yDot = (y - yPrev) / dt.
  Model::State sn;
  sn.t = 0.0;
  sn.y = yPrev;
  sn.phi = cardiacInput.phi0;
  const Model::State snm1 = sn;

  Model::DenseVector x(Vars::NumberOfVariables);
  x.setZero();
  x[Vars::RadialDisplacement] = y;
  x[Vars::Porosity] = cardiacInput.phi0;
  x[Vars::FiberDeformation] = ec;

  Model::EvalData data;
  dynamicSystem.buildEvalData(x, sn, snm1, dt, dt, data);

  // CCMLC2014 membrane law: (d0 / R0) sqrt(C) Sigma = p_v C.
  const Real R0 = cardiacInput.R0;
  const Real sqrtC = 1.0 + y / R0;
  const Real C = sqrtC * sqrtC;
  const Real v = (y - yPrev) / dt;

  Real sigmaPassive = 0.0;
  Real dSigmaPassive = 0.0;
  Heart::CCMLC2014PassiveLaw<Real> membraneLaw;
  membraneLaw(cardiacInput.passiveEnergy, C, 2.0 * sqrtC / R0, sigmaPassive, dSigmaPassive);

  const Real sigmaViscous =
    cardiacInput.eta * v / R0 * (2.0 * sqrtC + 4.0 * std::pow(sqrtC, -11.0));

  const Real e = 0.5 * (C - 1.0);
  const Real h = 1.0 + 2.0 * ec;
  const Real sigmaActive = cardiacInput.Es / (h * h) * (e - ec);

  const Real expected =
    cardiacInput.d0 / R0 * (sigmaPassive + sigmaViscous + sigmaActive) / sqrtC;

  EXPECT_NEAR(data.J, 1.0, 1e-15);
  EXPECT_NEAR(data.wallPressure, expected, 1e-3 * std::abs(expected));
  EXPECT_NEAR(data.strain1D, e, 1e-3 * std::abs(e));
}

/// @brief Verifies the averaged multiplier equals the volume average of the pointwise field recovered from radial equilibrium.
TEST(PoroelasticSphereTest, AveragedMultiplierMatchesPointwiseRecovery)
{
  auto cardiacInput = makeGenericCardiacInput();
  cardiacInput.R0 = 1.0;
  cardiacInput.d0 = 0.5;
  cardiacInput.phi0 = 0.1;
  cardiacInput.eta = 0.0;
  cardiacInput.wallQuadraturePoints = 16;

  const Real phi = 0.2;
  const Real J = 1.0 - cardiacInput.phi0 + phi;
  const Real rin = 1.3;
  const Real ec = 0.05;
  const auto data = evaluateSteadyWall(cardiacInput, rin - 1.0, phi, ec);
  const Real Pv = data.wallPressure;
  const Real sigma = data.active.activeStressOneDimensional;

  // sigma_rr(r) = -P_v + int_{r_in}^r 2 (sigma_tt - sigma_rr) / s ds,
  // lambda(r) = sigma_rr - sigma_M,rr, averaged over the reference wall.
  const Real Rin = cardiacInput.R0;
  const Real Rout = cardiacInput.R0 + cardiacInput.d0;
  PassiveLaw law;
  auto stress = [&](Real R, Real& sigmaMrr, Real& diff)
  {
    const Real r = std::cbrt(rin * rin * rin + J * (R * R * R - Rin * Rin * Rin));
    const Real lt = r / R;
    const Real lr = J / (lt * lt);
    const auto P = law(cardiacInput.passiveEnergy, lr, lt);
    const Real Ptt = P.Ptt + 0.5 * sigma * lt;
    sigmaMrr = lr * P.PRR / J;
    diff = (lt * Ptt - lr * P.PRR) / J;
    return r;
  };

  const size_t n = 2000;
  const Real dR = (Rout - Rin) / static_cast<Real>(n);
  Real sigmaRR = -Pv;
  Real lambdaSum = 0.0;
  Real rPrev, mPrev, dPrev;
  rPrev = stress(Rin, mPrev, dPrev);
  lambdaSum += 0.5 * (sigmaRR - mPrev) * Rin * Rin;
  for (size_t i = 1; i <= n; ++i)
  {
    const Real R = Rin + dR * static_cast<Real>(i);
    Real m, d;
    const Real r = stress(R, m, d);
    sigmaRR += 0.5 * (2.0 * dPrev / rPrev + 2.0 * d / r) * (r - rPrev);
    const Real weight = (i == n) ? 0.5 : 1.0;
    lambdaSum += weight * (sigmaRR - m) * R * R;
    rPrev = r;
    dPrev = d;
  }
  const Real lambdaAverage = 3.0 * lambdaSum * dR / (Rout * Rout * Rout - Rin * Rin * Rin);

  EXPECT_NEAR(sigmaRR, 0.0, 1e-6 * std::abs(Pv));
  EXPECT_NEAR(data.lambdaBar, lambdaAverage, 1e-5 * std::abs(lambdaAverage));
  EXPECT_LT(data.lambdaBar, 0.0);
}

/// @brief Verifies dynamic jacobian matches finite difference.
TEST(PoroelasticSphereTest, DynamicJacobianMatchesFiniteDifference)
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
        cardiacInput, candidateState, currentState, previousState, dt, relativePerturbation);
    EXPECT_LT(relativeError, 2e-3);
  }
}

/// @brief Verifies dynamic jacobian matches finite difference across data scales and the non-Newtonian path.
TEST(PoroelasticSphereTest, DynamicJacobianMatchesFiniteDifferenceAcrossDataScales)
{
  auto cardiacInput = makeGenericCardiacInput();
  cardiacInput.windkesselRheology = Vars::WindkesselRheology::CarreauYasuda;
  cardiacInput.proximalRadius = 0.015;
  cardiacInput.proximalLength = 0.4;
  cardiacInput.distalRadius = 0.0007;
  cardiacInput.distalLength = 0.004;
  cardiacInput.mu_0 = 0.04868;
  cardiacInput.mu_Inf = 0.003605;
  cardiacInput.lambda = 3.39;
  cardiacInput.n = 0.198;
  cardiacInput.yasuda = 1.235;

  const std::array<Real, 3> pressureScales{{0.2, 1.0, 5.0}};
  for (const Real pressureScale : pressureScales)
  {
    Model::DenseVector candidateState = makeCandidateState(pressureScale);
    Model::State currentState = makeCurrentState(pressureScale);
    Model::State previousState = makePreviousState(currentState, pressureScale);

    const Real relativeError = computeDynamicJacobianRelativeError(
        cardiacInput, candidateState, currentState, previousState, 1e-3, 1e-7);
    EXPECT_TRUE(std::isfinite(relativeError));
    EXPECT_LT(relativeError, 3e-3);
  }
}

/// @brief Verifies step converges and advances time.
TEST(PoroelasticSphereTest, StepConvergesAndAdvancesTime)
{
  auto cardiacInput = makeGenericCardiacInput();

  Model model(cardiacInput);
  model.setMaxIterations(200)
       .setAbsoluteTolerance(1e-8)
       .setRelativeTolerance(1e-8)
       .setStepTolerance(1e-10)
       .setDampingFactor(1.0);

  Model::State initial;
  initial.t = 0.0;
  initial.y = 0.0;
  initial.pv = cardiacInput.pAt(0.0) - 100.0;
  initial.par = 11000.0;
  initial.pd = 10000.0;
  model.initialize(initial);
  EXPECT_NEAR(model.getState().phi, cardiacInput.phi0, 1e-14);

  const Real dt = 1e-3;
  const auto report = model.step(dt);
  EXPECT_TRUE(report.converged);
  EXPECT_NEAR(model.getState().t, dt, 1e-14);
  EXPECT_TRUE(std::isfinite(model.getState().lambdaBar));
  EXPECT_TRUE(std::isfinite(model.getState().pf));
}

/// @brief Verifies example setup completes three cycles with physical state and systolic perfusion impediment.
TEST(PoroelasticSphereTest, ExampleSetupCompletesThreeCyclesWithPhysicalState)
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
  Real minPhi = 1.0, maxPhi = 0.0;
  Real minPv = 1e9, maxPv = -1e9;
  Real minLambdaBar = 1e9;
  Real minInflow = 1e9, maxInflow = -1e9;

  for (int step = 0; step < nsteps; ++step)
  {
    const auto report = model.step(dt);
    ASSERT_TRUE(report.converged)
      << "Example setup failed at step " << step << ", t = " << model.getState().t
      << ", residual = " << report.finalResidual;
    maxFinalResidual = std::max(maxFinalResidual, report.finalResidual);
    maxIterations = std::max(maxIterations, report.iterations);

    const auto& s = model.getState();
    minPhi = std::min(minPhi, s.phi);
    maxPhi = std::max(maxPhi, s.phi);
    minPv = std::min(minPv, s.pv);
    maxPv = std::max(maxPv, s.pv);
    minLambdaBar = std::min(minLambdaBar, s.lambdaBar);
    const Real inflow = cardiacInput.gammaAr * (s.par - s.pf);
    minInflow = std::min(minInflow, inflow);
    maxInflow = std::max(maxInflow, inflow);
  }

  const auto& state = model.getState();
  EXPECT_NEAR(state.t, nsteps * dt, 1e-12);
  EXPECT_TRUE(std::isfinite(state.y));
  EXPECT_TRUE(std::isfinite(state.phi));
  EXPECT_TRUE(std::isfinite(state.pv));
  EXPECT_TRUE(std::isfinite(state.par));
  EXPECT_TRUE(std::isfinite(state.pd));
  EXPECT_TRUE(std::isfinite(state.pf));
  EXPECT_GT(state.kc, 0.0);
  EXPECT_GT(state.w, 0.0);
  EXPECT_GT(minPhi, 0.0);
  EXPECT_LT(maxPhi, 1.0);
  EXPECT_GT(maxPv - minPv, 5.0e3);
  EXPECT_LT(minLambdaBar, 0.0);
  EXPECT_GT(maxInflow - minInflow, 0.1 * maxInflow);
  EXPECT_LT(maxFinalResidual, 5e-2);
  EXPECT_LT(maxIterations, 20u);
}
