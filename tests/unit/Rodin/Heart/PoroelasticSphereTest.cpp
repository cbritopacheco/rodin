#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

#include "Rodin/Heart/CCMLC2014.h"
#include "Rodin/Heart/PoroelasticSphere.h"
#include "Rodin/Heart/PoroelasticSphere/Numerics/DynamicSystem.h"

using namespace Rodin;

namespace
{
  using Model = Heart::PoroelasticSphereT<>;
  using PassiveLaw = Heart::UsykSpherePassiveLaw<Real>;
  using HolzapfelLaw = Heart::PoroelasticSpherePassiveLaw<Real>;
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

    cardiacInput.storageCompliance = 5.0e-6;
    cardiacInput.storagePressureScale = 1.0e4;
    cardiacInput.gammaAr = 7.0e-10;
    cardiacInput.permeability = 1.0e-11;
    cardiacInput.venousConductance = 2.67e-12;
    cardiacInput.qPerfusionExternal = [](Real) { return 4.0e-7; };
    cardiacInput.qArterialExternal = [](Real) { return 4.0e-7; };

    cardiacInput.Rp = 8.0e6;
    cardiacInput.Cp = 5.0e-9;
    cardiacInput.Rd = 1.0e8;
    cardiacInput.Cd = 1.0e-8;

    cardiacInput.Kat = 8.0e-7;
    cardiacInput.Kp = 5.0e-10;
    cardiacInput.Kar = 1.3e-5;

    cardiacInput.cavityCapacity = 5.0e-12;
    cardiacInput.absRegularization = 1e-14;

    cardiacInput.u = [](Real) { return 25.0; };
    cardiacInput.pAt = [](Real) { return 900.0; };
    cardiacInput.pSv = [](Real) { return 1000.0; };

    using PassiveEnergyType = std::decay_t<decltype(cardiacInput.passiveEnergy)>;
    typename PassiveEnergyType::Parameters passiveParameters;
    passiveParameters.c0 = 100.0;
    passiveParameters.brr = 1.24;
    passiveParameters.btt = 8.4;
    passiveParameters.brt = 0.6;
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

    cardiacInput.storageCompliance = 5.0e-6;
    cardiacInput.storagePressureScale = 1.0e4;
    cardiacInput.gammaAr = 2.0e-10;
    cardiacInput.permeability = 1.0e-11;
    cardiacInput.venousConductance = 2.67e-12;

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
    passiveParameters.c0 = 100.0;
    passiveParameters.brr = 1.24;
    passiveParameters.btt = 8.4;
    passiveParameters.brt = 0.6;
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
    const Model::Input& cardiacInput, Real y, Real phi, Real ec, Real pv = 0.0)
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
    x[Vars::VentricularPressure] = pv;

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
  Real computeDynamicJacobianRelativeError(const Model::Input& cardiacInput,
    const Model::DenseVector& candidateState, const Model::State& currentState,
    const Model::State& previousState, Real dt, Real relativePerturbation)
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

    return (analyticalJacobian - finiteDifferenceJacobian).norm() /
      std::max<Real>(finiteDifferenceJacobian.norm(), 1e-14);
  }
}

/// @brief Verifies the Usyk Piola components are the restricted-energy derivatives and their tangent.
TEST(PoroelasticSphereTest, UsykPiolaComponentsAndTangentMatchFiniteDifference)
{
  const auto cardiacInput = makeGenericCardiacInput();
  const auto& law = cardiacInput.passiveEnergy;
  const auto& p = law.getParameters();
  const Real lr = 0.83;
  const Real lt = 1.17;
  const Real h = 1e-6;

  auto energy = [&](Real a, Real b) {
    const Real er = a * a - 1.0, et = b * b - 1.0;
    const Real Q = 0.25 * (p.brr * er * er + 2.0 * p.btt * et * et + 2.0 * p.brt * er * et);
    return 0.5 * p.c0 * (std::exp(Q) - 1.0);
  };

  PassiveLaw stress;
  const auto P = stress(law, lr, lt);
  const auto Prp = stress(law, lr + h, lt);
  const auto Prm = stress(law, lr - h, lt);
  const auto Ptp = stress(law, lr, lt + h);
  const auto Ptm = stress(law, lr, lt - h);

  const Real scale = std::max<Real>(1.0, std::abs(P.dPtt_dlt));
  EXPECT_NEAR(P.PRR, (energy(lr + h, lt) - energy(lr - h, lt)) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.Ptt, 0.5 * (energy(lr, lt + h) - energy(lr, lt - h)) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.dPRR_dlr, (Prp.PRR - Prm.PRR) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.dPRR_dlt, (Ptp.PRR - Ptm.PRR) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.dPtt_dlr, (Prp.Ptt - Prm.Ptt) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.dPtt_dlt, (Ptp.Ptt - Ptm.Ptt) / (2 * h), 1e-6 * scale);
}

/// @brief Verifies the reduced-invariant passive tangent matches finite differences.
TEST(PoroelasticSphereTest, ReducedInvariantPiolaTangentMatchesFiniteDifference)
{
  Heart::HolzapfelReducedLaw<Real>::Parameters hp;
  hp.C0 = 1.9e3;
  hp.C1 = 1.1e-1;
  hp.C2 = 1.9e3;
  hp.C3 = 1.1e-1;
  const Heart::HolzapfelReducedLaw<Real> energy(hp);
  HolzapfelLaw law;
  const Real lr = 0.83, lt = 1.17, h = 1e-6;
  const auto P = law(energy, lr, lt);
  const auto Prp = law(energy, lr + h, lt);
  const auto Prm = law(energy, lr - h, lt);
  const auto Ptp = law(energy, lr, lt + h);
  const auto Ptm = law(energy, lr, lt - h);
  const Real scale = std::max<Real>(1.0, std::abs(P.dPtt_dlt));
  EXPECT_NEAR(P.dPRR_dlr, (Prp.PRR - Prm.PRR) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.dPRR_dlt, (Ptp.PRR - Ptm.PRR) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.dPtt_dlr, (Prp.Ptt - Prm.Ptt) / (2 * h), 1e-6 * scale);
  EXPECT_NEAR(P.dPtt_dlt, (Ptp.Ptt - Ptm.Ptt) / (2 * h), 1e-6 * scale);
}

/// @brief Verifies shell kinematics, the radial-stress closure and the pressure law in closed form.
TEST(PoroelasticSphereTest, ShellPressureLawMatchesClosedForm)
{
  auto cardiacInput = makeGenericCardiacInput();
  cardiacInput.eta = 0.0;
  const Real y = 2.5e-3, phi = 0.115, ec = 0.01;

  const auto probe = evaluateSteadyWall(cardiacInput, y, phi, ec, 0.0);
  const Real R0 = cardiacInput.R0, h0 = cardiacInput.d0;
  const Real lt = 1.0 + y / R0;
  const Real J = 1.0 - cardiacInput.phi0 + phi;
  const Real lr = J / (lt * lt);
  EXPECT_NEAR(probe.lt, lt, 1e-15);
  EXPECT_NEAR(probe.lr, lr, 1e-15);
  EXPECT_NEAR(probe.thickness, h0 * lr, 1e-15);

  // Static shell law with the multiplier eliminated:
  // p_v (1 + h0 lr / (R0 lt)) = 2 h0 / (R0 lt^2) [P_Mt - P_MR lr / lt].
  PassiveLaw stress;
  const auto P = stress(cardiacInput.passiveEnergy, lr, lt);
  const Real e = 0.5 * (lt * lt - 1.0);
  const Real hh = 1.0 + 2.0 * ec;
  const Real sigma = cardiacInput.Es / (hh * hh) * (e - ec);
  const Real PMt = P.Ptt + 0.5 * sigma * lt;
  const Real PMR = P.PRR;
  const Real pv = 2.0 * h0 / (R0 * lt * lt) * (PMt - PMR * lr / lt) /
    (1.0 + h0 * lr / (R0 * lt));

  const auto data = evaluateSteadyWall(cardiacInput, y, phi, ec, pv);
  EXPECT_NEAR(data.wallPressure, pv, 1e-10 * std::abs(pv));
  EXPECT_NEAR(data.lambdaBar, -0.5 * pv - PMR / (lt * lt), 1e-10 * std::abs(pv));

  const Real p0 = cardiacInput.storagePressureScale;
  const Real dPhiMax = cardiacInput.storageCompliance * p0;
  const Real pP = -p0 * std::log(1.0 - (phi - cardiacInput.phi0) / dPhiMax);
  EXPECT_NEAR(data.pP, pP, 1e-10 * std::abs(pP));
  EXPECT_NEAR(data.pf, pP + 0.5 * pv + PMR / (lt * lt), 1e-10 * std::abs(pv));
}

/// @brief Verifies the exponential storage law: p_P(phi0) = 0, p_P' = 1/C0 there, and its inverse.
TEST(PoroelasticSphereTest, ExponentialStorageLaw)
{
  const auto cardiacInput = makeGenericCardiacInput();
  Heart::PoroelasticSphere::Physics::PerfusionEvaluator<Model::Input> perfusion(
    cardiacInput);
  Real pP = 1.0, dpP = 0.0;
  perfusion.storage(cardiacInput.phi0, pP, dpP);
  EXPECT_NEAR(pP, 0.0, 1e-12);
  EXPECT_NEAR(dpP * cardiacInput.storageCompliance, 1.0, 1e-12);

  const Real p0 = cardiacInput.storagePressureScale;
  const Real dPhiMax = cardiacInput.storageCompliance * p0;
  for (const Real p : {-5.0e3, 2.0e3, 1.5e4})
  {
    const Real phi = cardiacInput.phi0 + dPhiMax * (1.0 - std::exp(-p / p0));
    perfusion.storage(phi, pP, dpP);
    EXPECT_NEAR(pP, p, 1e-9 * p0);
    EXPECT_NEAR(dpP, std::exp(p / p0) / cardiacInput.storageCompliance, 1e-9 * dpP);
  }
}

/// @brief Verifies the venous drainage scales inversely with the pore viscosity.
TEST(PoroelasticSphereTest, VenousDrainageScalesWithPoreViscosity)
{
  auto cardiacInput = makeGenericCardiacInput();
  const auto a = evaluateSteadyWall(cardiacInput, 2e-3, 0.11, 0.0, 4.0e3);
  cardiacInput.poreViscosity.eta0 *= 2.0;
  cardiacInput.poreViscosity.etaInf *= 2.0;
  const auto b = evaluateSteadyWall(cardiacInput, 2e-3, 0.11, 0.0, 4.0e3);
  EXPECT_NEAR(b.pf, a.pf, 1e-12 * std::abs(a.pf));
  EXPECT_NEAR(b.poreViscosity, 2.0 * a.poreViscosity, 1e-12 * a.poreViscosity);
  EXPECT_NEAR(b.qVenous, 0.5 * a.qVenous, 1e-10 * std::abs(a.qVenous));
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
    const Real relativeError = computeDynamicJacobianRelativeError(cardiacInput,
      candidateState, currentState, previousState, dt, relativePerturbation);
    EXPECT_LT(relativeError, 2e-3);
  }
}

/// @brief Verifies the exact Jacobian rows of the shell, blood-volume and proximal balances entrywise.
TEST(PoroelasticSphereTest, ShellAndPerfusionRowsMatchFiniteDifference)
{
  auto cardiacInput = makeGenericCardiacInput();
  Heart::PoroelasticSphere::Numerics::DynamicSystem<PassiveLaw, Model::Input> sys(
    cardiacInput);
  const auto x = makeCandidateState();
  const auto sn = makeCurrentState();
  const auto snm1 = makePreviousState(sn);
  const Real dt = 1e-3, t = sn.t + dt;

  Model::EvalData data;
  sys.buildEvalData(x, sn, snm1, t, dt, data);
  Model::DenseMatrix Jac;
  sys.evaluateJacobian(data, Jac, dt);

  const std::array<Index, 3> rows{
    {Index(Vars::RadialDisplacement), Index(Vars::Porosity), Index(Vars::ArterialPressure)}};
  for (const Index row : rows)
  {
    Real scale = 1e-12;
    for (Index j = 0; j < Vars::NumberOfVariables; ++j)
      scale = std::max<Real>(scale, std::abs(Jac(row, j)));
    for (Index j = 0; j < Vars::NumberOfVariables; ++j)
    {
      const Real hstep = 1e-7 * std::max<Real>(1.0, std::abs(x[j]));
      auto xp = x, xm = x;
      xp[j] += hstep;
      xm[j] -= hstep;
      Model::EvalData dp, dm;
      sys.buildEvalData(xp, sn, snm1, t, dt, dp);
      sys.buildEvalData(xm, sn, snm1, t, dt, dm);
      Model::DenseVector rp, rm;
      sys.evaluateResidual(dp, rp);
      sys.evaluateResidual(dm, rm);
      const Real fd = (rp[row] - rm[row]) / (2.0 * hstep);
      EXPECT_NEAR(Jac(row, j), fd, 1e-5 * scale) << "row " << row << " col " << j;
    }
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
  EXPECT_GT(model.getState().poreViscosity, 0.0);
}

/// @brief Verifies the example setup completes three cycles, conserves blood volume and raises the pore pressure in systole.
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

  const Real R0 = cardiacInput.R0, h0 = cardiacInput.d0;
  const Real Rin = R0 - 0.5 * h0, Rout = R0 + 0.5 * h0;
  const Real Vw0 =
    4.0 / 3.0 * std::numbers::pi_v<Real> * (Rout * Rout * Rout - Rin * Rin * Rin);

  const Real dt = 1e-3;
  const int nsteps = 3 * static_cast<int>(0.85 / dt);
  size_t maxIterations = 0;
  Real minPhi = 1.0, maxPhi = 0.0;
  Real minPv = 1e9, maxPv = -1e9;
  Real minPf = 1e9, maxPf = -1e9;
  Real firstStepBalance = 0.0;

  for (int step = 0; step < nsteps; ++step)
  {
    const Real phiOld = model.getState().phi;
    const auto report = model.step(dt);
    ASSERT_TRUE(report.converged)
      << "Example setup failed at step " << step << ", t = " << model.getState().t
      << ", residual = " << report.finalResidual;
    maxIterations = std::max(maxIterations, report.iterations);

    const auto& s = model.getState();
    minPhi = std::min(minPhi, s.phi);
    maxPhi = std::max(maxPhi, s.phi);
    minPv = std::min(minPv, s.pv);
    maxPv = std::max(maxPv, s.pv);
    minPf = std::min(minPf, s.pf);
    maxPf = std::max(maxPf, s.pf);
    if (step == 0) // Backward Euler start.
      firstStepBalance = std::abs(Vw0 * (s.phi - phiOld) / dt - (s.qArterial - s.qVenous)) /
        std::max<Real>(std::abs(s.qArterial), 1e-12);
  }

  const auto& state = model.getState();
  EXPECT_NEAR(state.t, nsteps * dt, 1e-12);
  EXPECT_TRUE(std::isfinite(state.pf));
  EXPECT_GT(state.kc, 0.0);
  EXPECT_GT(minPhi, 0.0);
  EXPECT_LT(maxPhi,
    cardiacInput.phi0 + cardiacInput.storageCompliance * cardiacInput.storagePressureScale);
  EXPECT_GT(maxPv - minPv, 5.0e3);
  EXPECT_GT(maxPf - minPf, 0.3 * (maxPv - minPv));
  EXPECT_LT(firstStepBalance, 1e-6);
  EXPECT_LT(maxIterations, 20u);
}
