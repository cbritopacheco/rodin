/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <array>
#include <type_traits>

#include <Rodin/Solid.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>

namespace Rodin::Tests::Unit
{
  // ========================================================================
  // KinematicState tests
  // ========================================================================

  /// @brief Verifies identity for solid kinematic state by checking tolerance-based numerical results.
  TEST(Rodin_Solid_KinematicState, Identity)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);

    EXPECT_NEAR(state.getJacobian(), 1.0, 1e-14);
    EXPECT_NEAR(state.getLogJacobian(), 0.0, 1e-14);

    const auto& F = state.getDeformationGradient();
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(F(i, j), (i == j) ? 1.0 : 0.0, 1e-14);
    }

    const auto& C = state.getRightCauchyGreenTensor();
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(C(i, j), (i == j) ? 1.0 : 0.0, 1e-14);
    }
  }

  /// @brief Verifies simple shear 2 D for solid kinematic state by checking tolerance-based numerical results.
  TEST(Rodin_Solid_KinematicState, SimpleShear2D)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.0;
    H(0, 1) = 0.1;
    H(1, 0) = 0.0;
    H(1, 1) = 0.0;
    state.setDisplacementGradient(H);

    // F = [[1, 0.1], [0, 1]], det(F) = 1
    EXPECT_NEAR(state.getJacobian(), 1.0, 1e-14);

    const auto& F = state.getDeformationGradient();
    EXPECT_NEAR(F(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(F(0, 1), 0.1, 1e-14);
    EXPECT_NEAR(F(1, 0), 0.0, 1e-14);
    EXPECT_NEAR(F(1, 1), 1.0, 1e-14);

    // C = F^T F = [[1, 0.1], [0.1, 1.01]]
    const auto& C = state.getRightCauchyGreenTensor();
    EXPECT_NEAR(C(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(C(0, 1), 0.1, 1e-14);
    EXPECT_NEAR(C(1, 0), 0.1, 1e-14);
    EXPECT_NEAR(C(1, 1), 1.01, 1e-14);
  }

  /// @brief Verifies uniform expansion 3 D for solid kinematic state by checking tolerance-based numerical results.
  TEST(Rodin_Solid_KinematicState, UniformExpansion3D)
  {
    Solid::KinematicState state(3);
    const Real alpha = 0.5;
    Math::SpatialMatrix<Real> I(3, 3);
    I.setIdentity();
    Math::SpatialMatrix<Real> H = alpha * I;
    state.setDisplacementGradient(H);

    // F = (1 + alpha) I, J = (1 + alpha)^3
    const Real expected_J = std::pow(1.0 + alpha, 3);
    EXPECT_NEAR(state.getJacobian(), expected_J, 1e-12);
    EXPECT_NEAR(state.getLogJacobian(), std::log(expected_J), 1e-12);
  }

  /// @brief Verifies inverse consistency for solid kinematic state by checking tolerance-based numerical results.
  TEST(Rodin_Solid_KinematicState, InverseConsistency)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(1, 0) = -0.02;
    H(1, 1) = 0.2;
    state.setDisplacementGradient(H);

    const auto& F = state.getDeformationGradient();
    const auto& Finv = state.getDeformationGradientInverse();

    // F * F^{-1} should be identity
    Math::SpatialMatrix<Real> product = F * Finv;
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(product(i, j), (i == j) ? 1.0 : 0.0, 1e-12);
    }
  }

  // ========================================================================
  // Invariants tests
  // ========================================================================

  /// @brief Verifies identity for solid isotropic invariants by checking tolerance-based numerical results.
  TEST(Rodin_Solid_IsotropicInvariants, Identity)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);

    Solid::IsotropicInvariants inv;
    inv.setState(state);

    // For identity: C = I, I1 = 2, I2 = 1, I3 = 1
    EXPECT_NEAR(inv.getFirstInvariant(), 2.0, 1e-14);
    EXPECT_NEAR(inv.getSecondInvariant(), 1.0, 1e-14);
    EXPECT_NEAR(inv.getThirdInvariant(), 1.0, 1e-14);
  }

  /// @brief Verifies identity 3 D for solid isotropic invariants by checking tolerance-based numerical results.
  TEST(Rodin_Solid_IsotropicInvariants, Identity3D)
  {
    Solid::KinematicState state(3);
    Math::SpatialMatrix<Real> H(3, 3);
    H.setZero();
    state.setDisplacementGradient(H);

    Solid::IsotropicInvariants inv;
    inv.setState(state);

    // For identity 3D: C = I, I1 = 3, I2 = 3, I3 = 1
    EXPECT_NEAR(inv.getFirstInvariant(), 3.0, 1e-14);
    EXPECT_NEAR(inv.getSecondInvariant(), 3.0, 1e-14);
    EXPECT_NEAR(inv.getThirdInvariant(), 1.0, 1e-14);
  }

  /// @brief Verifies aligned fiber for solid fiber invariants by checking tolerance-based numerical results.
  TEST(Rodin_Solid_FiberInvariants, AlignedFiber)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.2;
    H(0, 1) = 0.0;
    H(1, 0) = 0.0;
    H(1, 1) = 0.0;
    state.setDisplacementGradient(H);

    Math::SpatialVector<Real> a0(2);
    a0[0] = 1.0;
    a0[1] = 0.0;
    Solid::FiberInvariants fib(a0);
    fib.setState(state);

    // F = [[1.2, 0],[0, 1]], C = [[1.44, 0],[0, 1]]
    // I4 = a0 . C . a0 = 1.44
    EXPECT_NEAR(fib.getFourthInvariant(), 1.44, 1e-12);

    // I5 = a0 . C^2 . a0 = 1.44^2 = 2.0736
    EXPECT_NEAR(fib.getFifthInvariant(), 1.44 * 1.44, 1e-12);
  }

  /// @brief Verifies normalized direction and strain for solid fiber kinematics by checking tolerance-based numerical results.
  TEST(Rodin_Solid_FiberKinematics, NormalizedDirectionAndStrain)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.2;
    H(0, 1) = 0.0;
    H(1, 0) = 0.0;
    H(1, 1) = 0.0;
    state.setDisplacementGradient(H);

    Solid::ConstitutivePoint cp(state);
    Math::SpatialVector<Real> direction(2);
    direction[0] = 2.0;
    direction[1] = 0.0;
    cp.set<Solid::Tags::FiberDirection>(direction);

    Solid::FiberKinematics fiber(cp);

    EXPECT_NEAR(fiber.direction()[0], 1.0, 1e-14);
    EXPECT_NEAR(fiber.direction()[1], 0.0, 1e-14);
    EXPECT_NEAR(fiber.I4(), 1.44, 1e-12);
    EXPECT_NEAR(fiber.strain(), 0.22, 1e-12);
  }

  // ========================================================================
  // NeoHookean tests
  // ========================================================================

  /// @brief Verifies zero deformation for solid neo hookean by checking tolerance-based numerical results.
  TEST(Rodin_Solid_NeoHookean, ZeroDeformation)
  {
    const Real lambda = 1.0, mu = 0.5;
    Solid::NeoHookean law(lambda, mu);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    // At zero deformation, P should be zero
    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(P(i, j), 0.0, 1e-14);
    }

    // Energy should also be zero
    EXPECT_NEAR(prepared.getStrainEnergyDensity(), 0.0, 1e-14);
  }

  /// @brief Verifies stress symmetry under pure stretch for solid neo hookean by checking tolerance-based numerical results.
  TEST(Rodin_Solid_NeoHookean, StressSymmetryUnderPureStretch)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.1;
    H(0, 1) = 0.0;
    H(1, 0) = 0.0;
    H(1, 1) = 0.1;
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);

    // Under uniform stretch, P should be diagonal
    EXPECT_NEAR(P(0, 1), 0.0, 1e-13);
    EXPECT_NEAR(P(1, 0), 0.0, 1e-13);
    // And P(0,0) == P(1,1) by symmetry
    EXPECT_NEAR(P(0, 0), P(1, 1), 1e-13);
  }

  /// @brief Verifies tangent finite difference for solid neo hookean by checking tolerance-based numerical results.
  TEST(Rodin_Solid_NeoHookean, TangentFiniteDifference)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(1, 0) = -0.02;
    H(1, 1) = 0.15;
    state.setDisplacementGradient(H);

    auto prepared = law;
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    // Test tangent via finite differences
    Math::SpatialMatrix<Real> dF(2, 2);
    dF(0, 0) = 0.3;
    dF(0, 1) = -0.2;
    dF(1, 0) = 0.1;
    dF(1, 1) = 0.4;

    Math::SpatialMatrix<Real> dP_analytical;
    prepared.getMaterialTangent(dP_analytical, dF);

    // Finite difference approximation
    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(2);
    state_plus.setDisplacementGradient(H_plus);
    auto prepared_plus = law;
    Solid::ConstitutivePoint cp_plus(state_plus);
    prepared_plus.setConstitutivePoint(cp_plus);
    Math::SpatialMatrix<Real> P_plus;
    prepared_plus.getFirstPiolaKirchhoffStress(P_plus);

    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);

    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-5);
    }
  }

  // ========================================================================
  // SaintVenantKirchhoff tests
  // ========================================================================

  /// @brief Verifies zero deformation for solid saint venant kirchhoff by checking tolerance-based numerical results.
  TEST(Rodin_Solid_SaintVenantKirchhoff, ZeroDeformation)
  {
    const Real lambda = 1.0, mu = 0.5;
    Solid::SaintVenantKirchhoff law(lambda, mu);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    // At zero deformation, P should be zero (E = 0, S = 0)
    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(P(i, j), 0.0, 1e-14);
    }

    EXPECT_NEAR(prepared.getStrainEnergyDensity(), 0.0, 1e-14);
  }

  /// @brief Verifies tangent finite difference for solid saint venant kirchhoff by checking tolerance-based numerical results.
  TEST(Rodin_Solid_SaintVenantKirchhoff, TangentFiniteDifference)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::SaintVenantKirchhoff law(lambda, mu);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(1, 0) = -0.02;
    H(1, 1) = 0.15;
    state.setDisplacementGradient(H);

    auto prepared = law;
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> dF(2, 2);
    dF(0, 0) = 0.3;
    dF(0, 1) = -0.2;
    dF(1, 0) = 0.1;
    dF(1, 1) = 0.4;

    Math::SpatialMatrix<Real> dP_analytical;
    prepared.getMaterialTangent(dP_analytical, dF);

    // Finite difference
    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(2);
    state_plus.setDisplacementGradient(H_plus);
    auto prepared_plus = law;
    Solid::ConstitutivePoint cp_plus(state_plus);
    prepared_plus.setConstitutivePoint(cp_plus);
    Math::SpatialMatrix<Real> P_plus;
    prepared_plus.getFirstPiolaKirchhoffStress(P_plus);

    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);

    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-5);
    }
  }

  /// @brief For small deformations, SVK should behave like linear elasticity.
  TEST(Rodin_Solid_SaintVenantKirchhoff, ReducesToLinearElasticity)
  {
    // For small deformations, SVK should behave like linear elasticity
    const Real lambda = 1.0, mu = 0.5;
    Solid::SaintVenantKirchhoff law(lambda, mu);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    const Real eps = 1e-6;
    H(0, 0) = eps;
    H(0, 1) = 0.5 * eps;
    H(1, 0) = 0.5 * eps;
    H(1, 1) = 2.0 * eps;
    state.setDisplacementGradient(H);

    auto prepared = law;
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);

    // For infinitesimal strain, P ~ sigma = lambda tr(epsilon) I + 2 mu epsilon
    // where epsilon = 0.5 (H + H^T)
    Math::SpatialMatrix<Real> epsilon = 0.5 * (H + H.transpose());
    Math::SpatialMatrix<Real> I2(2, 2);
    I2.setIdentity();
    Math::SpatialMatrix<Real> sigma = lambda * epsilon.trace() * I2 + 2.0 * mu * epsilon;

    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(P(i, j), sigma(i, j), 1e-10);
    }
  }

  // ========================================================================
  // MooneyRivlin tests
  // ========================================================================

  /// @brief Verifies zero deformation for solid mooney rivlin by checking tolerance-based numerical results.
  TEST(Rodin_Solid_MooneyRivlin, ZeroDeformation)
  {
    const Real c1 = 0.5, c2 = 0.1, kappa = 10.0;
    Solid::MooneyRivlin law(c1, c2, kappa);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    // At zero deformation, P should be zero
    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(P(i, j), 0.0, 1e-12);
    }
  }

  /// @brief Verifies tangent finite difference for solid mooney rivlin by checking tolerance-based numerical results.
  TEST(Rodin_Solid_MooneyRivlin, TangentFiniteDifference)
  {
    const Real c1 = 0.5, c2 = 0.2, kappa = 5.0;
    Solid::MooneyRivlin law(c1, c2, kappa);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(1, 0) = -0.02;
    H(1, 1) = 0.15;
    state.setDisplacementGradient(H);

    auto prepared = law;
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> dF(2, 2);
    dF(0, 0) = 0.3;
    dF(0, 1) = -0.2;
    dF(1, 0) = 0.1;
    dF(1, 1) = 0.4;

    Math::SpatialMatrix<Real> dP_analytical;
    prepared.getMaterialTangent(dP_analytical, dF);

    // Finite difference
    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(2);
    state_plus.setDisplacementGradient(H_plus);
    auto prepared_plus = law;
    Solid::ConstitutivePoint cp_plus(state_plus);
    prepared_plus.setConstitutivePoint(cp_plus);
    Math::SpatialMatrix<Real> P_plus;
    prepared_plus.getFirstPiolaKirchhoffStress(P_plus);

    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);

    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-4);
    }
  }

  // ========================================================================
  // Holzapfel-Ogden and active contraction tests
  // ========================================================================

  /// @brief Verifies zero deformation stress free for solid holzapfel ogden by checking tolerance-based numerical results.
  TEST(Rodin_Solid_HolzapfelOgden, ZeroDeformationStressFree)
  {
    Solid::HolzapfelOgden law(0.5, 0.2, 0.1, 1.5, 0.3, 2.0, 10.0);

    Solid::KinematicState state(3);
    Math::SpatialMatrix<Real> H(3, 3);
    H.setZero();
    state.setDisplacementGradient(H);

    Solid::ConstitutivePoint cp(state);
    Math::SpatialVector<Real> fiber(3);
    fiber[0] = 1.0;
    fiber[1] = 0.0;
    fiber[2] = 0.0;
    cp.set<Solid::Tags::FiberDirection>(fiber);

    auto prepared = law;
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(P(i, j), 0.0, 1e-12);
    }
  }

  /// @brief Verifies tangent finite difference 3 D for solid holzapfel ogden by checking tolerance-based numerical results.
  TEST(Rodin_Solid_HolzapfelOgden, TangentFiniteDifference3D)
  {
    Solid::HolzapfelOgden law(0.5, 0.2, 0.1, 1.5, 0.3, 2.0, 10.0);

    Solid::KinematicState state(3);
    Math::SpatialMatrix<Real> H(3, 3);
    H(0, 0) = 0.08;
    H(0, 1) = 0.03;
    H(0, 2) = 0.01;
    H(1, 0) = -0.02;
    H(1, 1) = 0.06;
    H(1, 2) = 0.04;
    H(2, 0) = 0.01;
    H(2, 1) = -0.03;
    H(2, 2) = 0.05;
    state.setDisplacementGradient(H);

    Solid::ConstitutivePoint cp(state);
    Math::SpatialVector<Real> fiber(3);
    fiber[0] = 0.8;
    fiber[1] = 0.4;
    fiber[2] = 0.2;
    cp.set<Solid::Tags::FiberDirection>(fiber);

    auto prepared = law;
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> dF(3, 3);
    dF(0, 0) = 0.2;
    dF(0, 1) = -0.1;
    dF(0, 2) = 0.05;
    dF(1, 0) = 0.03;
    dF(1, 1) = 0.15;
    dF(1, 2) = -0.07;
    dF(2, 0) = -0.04;
    dF(2, 1) = 0.08;
    dF(2, 2) = 0.12;

    Math::SpatialMatrix<Real> dP;
    prepared.getMaterialTangent(dP, dF);

    const Real eps = 1e-5;
    Solid::KinematicState statePlus(3);
    statePlus.setDisplacementGradient(H + eps * dF);
    Solid::ConstitutivePoint cpPlus(statePlus);
    cpPlus.set<Solid::Tags::FiberDirection>(fiber);
    auto preparedPlus = law;
    preparedPlus.setConstitutivePoint(cpPlus);

    Math::SpatialMatrix<Real> P, PPlus;
    prepared.getFirstPiolaKirchhoffStress(P);
    preparedPlus.getFirstPiolaKirchhoffStress(PPlus);
    Solid::KinematicState stateMinus(3);
    stateMinus.setDisplacementGradient(H + (-eps) * dF);
    Solid::ConstitutivePoint cpMinus(stateMinus);
    cpMinus.set<Solid::Tags::FiberDirection>(fiber);
    auto preparedMinus = law;
    preparedMinus.setConstitutivePoint(cpMinus);
    Math::SpatialMatrix<Real> PMinus;
    preparedMinus.getFirstPiolaKirchhoffStress(PMinus);
    // A centered difference avoids the forward oracle's truncation error.
    Math::SpatialMatrix<Real> dPfd = (0.5 / eps) * PPlus + (-0.5 / eps) * PMinus;

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(dP(i, j), dPfd(i, j), 1e-8);
    }
  }

  /// @brief Verifies adds fiber stress for solid active contraction by checking tolerance-based numerical results.
  TEST(Rodin_Solid_ActiveContraction, AddsFiberStress)
  {
    Solid::NeoHookean passive(0.0, 0.0);
    Solid::ActiveFiberLaw::Parameters activeInput;
    activeInput.stiffness = 10.0;
    activeInput.initial.extension = 0.0;
    Solid::ActiveFiberLaw active(activeInput);
    Solid::ActiveContraction law(passive, active);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    H(0, 0) = 0.1;
    state.setDisplacementGradient(H);

    Solid::ConstitutivePoint cp(state);
    Math::SpatialVector<Real> fiber(2);
    fiber[0] = 1.0;
    fiber[1] = 0.0;
    cp.set<Solid::Tags::FiberDirection>(fiber);
    cp.set<Solid::Tags::ActiveExtension>(0.0);

    auto prepared = law;
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);

    const Real strain = 0.5 * (1.1 * 1.1 - 1.0);
    const Real stress = activeInput.stiffness * strain;
    EXPECT_NEAR(P(0, 0), 1.1 * stress, 1e-12);
    EXPECT_NEAR(P(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(P(1, 0), 0.0, 1e-12);
    EXPECT_NEAR(P(1, 1), 0.0, 1e-12);
  }

  /// @brief Pure active contribution: zero passive parameters so the FD check.
  TEST(Rodin_Solid_ActiveContraction, TangentFiniteDifference2D)
  {
    // Pure active contribution: zero passive parameters so the FD check
    // isolates the active fiber tangent.
    Solid::NeoHookean passive(0.0, 0.0);
    Solid::ActiveFiberLaw::Parameters activeInput;
    activeInput.stiffness = 50.0;
    activeInput.initial.extension = 0.0;
    Solid::ActiveFiberLaw active(activeInput);
    Solid::ActiveContraction law(passive, active);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.07;
    H(0, 1) = 0.02;
    H(1, 0) = -0.03;
    H(1, 1) = 0.05;
    state.setDisplacementGradient(H);

    Solid::ConstitutivePoint cp(state);
    Math::SpatialVector<Real> fiber(2);
    fiber[0] = 0.8;
    fiber[1] = 0.6;
    cp.set<Solid::Tags::FiberDirection>(fiber);
    cp.set<Solid::Tags::ActiveExtension>(-0.06);

    auto prepared = law;
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> dF(2, 2);
    dF(0, 0) = 0.21;
    dF(0, 1) = -0.13;
    dF(1, 0) = 0.07;
    dF(1, 1) = 0.18;

    Math::SpatialMatrix<Real> dP;
    prepared.getMaterialTangent(dP, dF);

    const Real eps = 1e-7;
    Solid::KinematicState statePlus(2);
    statePlus.setDisplacementGradient(H + eps * dF);
    Solid::ConstitutivePoint cpPlus(statePlus);
    cpPlus.set<Solid::Tags::FiberDirection>(fiber);
    cpPlus.set<Solid::Tags::ActiveExtension>(-0.06);
    auto preparedPlus = law;
    preparedPlus.setConstitutivePoint(cpPlus);

    Math::SpatialMatrix<Real> P, PPlus;
    prepared.getFirstPiolaKirchhoffStress(P);
    preparedPlus.getFirstPiolaKirchhoffStress(PPlus);
    Math::SpatialMatrix<Real> dPfd = (1.0 / eps) * PPlus + (-1.0 / eps) * P;

    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dP(i, j), dPfd(i, j), 1e-6);
    }
  }

  /// @brief Mixed passive + active in 3D, fiber not aligned with an axis.
  TEST(Rodin_Solid_ActiveContraction, TangentFiniteDifference3DWithPassive)
  {
    // Mixed passive + active in 3D, fiber not aligned with an axis.
    Solid::NeoHookean passive(20.0, 8.0);
    Solid::ActiveFiberLaw::Parameters activeInput;
    activeInput.stiffness = 30.0;
    activeInput.initial.extension = 0.0;
    Solid::ActiveFiberLaw active(activeInput);
    Solid::ActiveContraction law(passive, active);

    Solid::KinematicState state(3);
    Math::SpatialMatrix<Real> H(3, 3);
    H(0, 0) = 0.04;
    H(0, 1) = 0.02;
    H(0, 2) = -0.01;
    H(1, 0) = -0.03;
    H(1, 1) = 0.06;
    H(1, 2) = 0.02;
    H(2, 0) = 0.01;
    H(2, 1) = -0.02;
    H(2, 2) = 0.05;
    state.setDisplacementGradient(H);

    Solid::ConstitutivePoint cp(state);
    Math::SpatialVector<Real> fiber(3);
    fiber[0] = 0.7;
    fiber[1] = 0.5;
    fiber[2] = 0.3;
    cp.set<Solid::Tags::FiberDirection>(fiber);
    cp.set<Solid::Tags::ActiveExtension>(-0.04);

    auto prepared = law;
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> dF(3, 3);
    dF(0, 0) = 0.18;
    dF(0, 1) = -0.09;
    dF(0, 2) = 0.05;
    dF(1, 0) = 0.04;
    dF(1, 1) = 0.13;
    dF(1, 2) = -0.07;
    dF(2, 0) = -0.05;
    dF(2, 1) = 0.06;
    dF(2, 2) = 0.10;

    Math::SpatialMatrix<Real> dP;
    prepared.getMaterialTangent(dP, dF);

    const Real eps = 1e-7;
    Solid::KinematicState statePlus(3);
    statePlus.setDisplacementGradient(H + eps * dF);
    Solid::ConstitutivePoint cpPlus(statePlus);
    cpPlus.set<Solid::Tags::FiberDirection>(fiber);
    cpPlus.set<Solid::Tags::ActiveExtension>(-0.04);
    auto preparedPlus = law;
    preparedPlus.setConstitutivePoint(cpPlus);

    Math::SpatialMatrix<Real> P, PPlus;
    prepared.getFirstPiolaKirchhoffStress(P);
    preparedPlus.getFirstPiolaKirchhoffStress(PPlus);
    Math::SpatialMatrix<Real> dPfd = (1.0 / eps) * PPlus + (-1.0 / eps) * P;

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(dP(i, j), dPfd(i, j), 1e-6);
    }
  }

  /// @brief Verifies the displacement tangent at fixed active fields.
  TEST(Rodin_Solid_ActiveContraction, DynamicTangentFiniteDifference2D)
  {
    // The displacement block holds all three independent active fields fixed.
    Solid::NeoHookean passive(0.0, 0.0);
    Solid::ActiveFiberLaw::Parameters activeInput;
    activeInput.stiffness = 200.0;
    activeInput.damping = 0.5;
    activeInput.destructionRate = 0.4;
    activeInput.crossBridgeStiffness = 100.0;
    activeInput.contractility = 80.0;
    Solid::ActiveFiberLaw active(activeInput);
    Solid::ActiveContraction law(passive, active);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.06;
    H(0, 1) = 0.02;
    H(1, 0) = -0.03;
    H(1, 1) = 0.04;
    state.setDisplacementGradient(H);

    auto buildPoint = [&](const Solid::KinematicState& s) {
      Solid::ConstitutivePoint cp(s);
      Math::SpatialVector<Real> fiber(2);
      fiber[0] = 0.8;
      fiber[1] = 0.6;
      cp.set<Solid::Tags::FiberDirection>(fiber);
      cp.set<Solid::Tags::TimeStep>(0.01);
      cp.set<Solid::Tags::PreviousActiveExtension>(-0.04);
      cp.set<Solid::Tags::PreviousActiveGamma>(2.0);
      cp.set<Solid::Tags::PreviousActiveBeta>(1.5);
      cp.set<Solid::Tags::ElectricalActivation>(0.7);
      cp.set<Solid::Tags::ActiveExtension>(-0.025);
      cp.set<Solid::Tags::ActiveGamma>(1.8);
      cp.set<Solid::Tags::ActiveBeta>(0.95);
      return cp;
    };

    Solid::ConstitutivePoint cp = buildPoint(state);
    auto prepared = law;
    prepared.setConstitutivePoint(cp);
    EXPECT_TRUE(prepared.isDynamic());
    EXPECT_DOUBLE_EQ(prepared.getActiveExtension(), -0.025);
    EXPECT_DOUBLE_EQ(prepared.getActiveState().gamma, 1.8);
    EXPECT_DOUBLE_EQ(prepared.getActiveState().beta, 0.95);

    Math::SpatialMatrix<Real> dF(2, 2);
    dF(0, 0) = 0.21;
    dF(0, 1) = -0.13;
    dF(1, 0) = 0.07;
    dF(1, 1) = 0.18;

    Math::SpatialMatrix<Real> dP;
    prepared.getMaterialTangent(dP, dF);

    const Real eps = 1e-7;
    Solid::KinematicState statePlus(2);
    statePlus.setDisplacementGradient(H + eps * dF);
    Solid::ConstitutivePoint cpPlus = buildPoint(statePlus);
    auto preparedPlus = law;
    preparedPlus.setConstitutivePoint(cpPlus);

    Math::SpatialMatrix<Real> P, PPlus;
    prepared.getFirstPiolaKirchhoffStress(P);
    preparedPlus.getFirstPiolaKirchhoffStress(PPlus);
    Math::SpatialMatrix<Real> dPfd = (1.0 / eps) * PPlus + (-1.0 / eps) * P;

    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dP(i, j), dPfd(i, j), 5e-5);
    }
  }

  /// @brief Verifies dynamic tangent finite difference 3 D with passive for solid active contraction by checking tolerance-based numerical results, true predicates.
  TEST(Rodin_Solid_ActiveContraction, DynamicTangentFiniteDifference3DWithPassive)
  {
    Solid::NeoHookean passive(20.0, 8.0);
    Solid::ActiveFiberLaw::Parameters activeInput;
    activeInput.stiffness = 60.0;
    activeInput.damping = 0.3;
    activeInput.destructionRate = 0.5;
    activeInput.crossBridgeStiffness = 40.0;
    activeInput.contractility = 30.0;
    Solid::ActiveFiberLaw active(activeInput);
    Solid::ActiveContraction law(passive, active);

    Solid::KinematicState state(3);
    Math::SpatialMatrix<Real> H(3, 3);
    H(0, 0) = 0.04;
    H(0, 1) = 0.02;
    H(0, 2) = -0.01;
    H(1, 0) = -0.03;
    H(1, 1) = 0.05;
    H(1, 2) = 0.02;
    H(2, 0) = 0.01;
    H(2, 1) = -0.02;
    H(2, 2) = 0.03;
    state.setDisplacementGradient(H);

    auto buildPoint = [&](const Solid::KinematicState& s) {
      Solid::ConstitutivePoint cp(s);
      Math::SpatialVector<Real> fiber(3);
      fiber[0] = 0.7;
      fiber[1] = 0.5;
      fiber[2] = 0.3;
      cp.set<Solid::Tags::FiberDirection>(fiber);
      cp.set<Solid::Tags::TimeStep>(0.005);
      cp.set<Solid::Tags::PreviousActiveExtension>(-0.03);
      cp.set<Solid::Tags::PreviousActiveGamma>(1.5);
      cp.set<Solid::Tags::PreviousActiveBeta>(0.8);
      cp.set<Solid::Tags::ElectricalActivation>(1.2);
      cp.set<Solid::Tags::ActiveExtension>(-0.025);
      cp.set<Solid::Tags::ActiveGamma>(1.8);
      cp.set<Solid::Tags::ActiveBeta>(0.95);
      return cp;
    };

    Solid::ConstitutivePoint cp = buildPoint(state);
    auto prepared = law;
    prepared.setConstitutivePoint(cp);
    EXPECT_TRUE(prepared.isDynamic());

    Math::SpatialMatrix<Real> dF(3, 3);
    dF(0, 0) = 0.18;
    dF(0, 1) = -0.09;
    dF(0, 2) = 0.05;
    dF(1, 0) = 0.04;
    dF(1, 1) = 0.13;
    dF(1, 2) = -0.07;
    dF(2, 0) = -0.05;
    dF(2, 1) = 0.06;
    dF(2, 2) = 0.10;

    Math::SpatialMatrix<Real> dP;
    prepared.getMaterialTangent(dP, dF);

    const Real eps = 1e-7;
    Solid::KinematicState statePlus(3);
    statePlus.setDisplacementGradient(H + eps * dF);
    Solid::ConstitutivePoint cpPlus = buildPoint(statePlus);
    auto preparedPlus = law;
    preparedPlus.setConstitutivePoint(cpPlus);

    Math::SpatialMatrix<Real> P, PPlus;
    prepared.getFirstPiolaKirchhoffStress(P);
    preparedPlus.getFirstPiolaKirchhoffStress(PPlus);
    Math::SpatialMatrix<Real> dPfd = (1.0 / eps) * PPlus + (-1.0 / eps) * P;

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(dP(i, j), dPfd(i, j), 5e-5);
    }
  }

  /// @brief Independent active fields remain unchanged away from equilibrium.
  TEST(Rodin_Solid_ActiveContraction, DynamicFieldsAreNotEliminated)
  {
    // Evaluation must preserve the independent fields, even away from equilibrium.
    Solid::NeoHookean passive(0.0, 0.0);
    Solid::ActiveFiberLaw::Parameters activeInput;
    activeInput.stiffness = 200.0;
    activeInput.damping = 0.5;
    activeInput.destructionRate = 0.4;
    activeInput.crossBridgeStiffness = 100.0;
    activeInput.contractility = 80.0;
    Solid::ActiveFiberLaw active(activeInput);
    Solid::ActiveContraction law(passive, active);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.10;
    H(0, 1) = 0.03;
    H(1, 0) = -0.02;
    H(1, 1) = 0.06;
    state.setDisplacementGradient(H);

    Solid::ConstitutivePoint cp(state);
    Math::SpatialVector<Real> fiber(2);
    fiber[0] = 1.0;
    fiber[1] = 0.0;
    cp.set<Solid::Tags::FiberDirection>(fiber);
    cp.set<Solid::Tags::TimeStep>(0.01);
    cp.set<Solid::Tags::PreviousActiveExtension>(-0.05);
    cp.set<Solid::Tags::PreviousActiveGamma>(2.0);
    cp.set<Solid::Tags::PreviousActiveBeta>(1.5);
    cp.set<Solid::Tags::ElectricalActivation>(0.9);
    cp.set<Solid::Tags::ActiveExtension>(-0.025);
    cp.set<Solid::Tags::ActiveGamma>(1.8);
    cp.set<Solid::Tags::ActiveBeta>(0.95);

    auto prepared = law;
    prepared.setConstitutivePoint(cp);

    EXPECT_TRUE(prepared.isDynamic());
    EXPECT_GT(std::abs(prepared.getActiveResponse().residual), 1e-3);
    EXPECT_EQ(
      prepared.getActiveResponse().residual, prepared.getActiveResidual().value[0]);
    EXPECT_DOUBLE_EQ(prepared.getActiveExtension(), -0.025);
  }

  /// @brief Verifies static stress derivative finite difference for solid active fiber law by checking tolerance-based numerical results.
  TEST(Rodin_Solid_ActiveFiberLaw, StaticStressDerivativeFiniteDifference)
  {
    Solid::ActiveFiberLaw::Parameters params;
    params.stiffness = 50.0;
    Solid::ActiveFiberLaw law(params);

    const Real e = 0.07;
    const Real c = -0.06;
    const Real eps = 1e-5;

    // Central differences: O(eps^2) truncation error.
    const Real fdDe = (law.stress(e + eps, c) - law.stress(e - eps, c)) / (2 * eps);
    EXPECT_NEAR(law.dStressDe(c), fdDe, 1e-7);

    const Real fdDc = (law.stress(e, c + eps) - law.stress(e, c - eps)) / (2 * eps);
    EXPECT_NEAR(law.dStressDc(e, c), fdDc, 1e-5);
  }

  /// @brief Verifies initial values for solid active fiber law by checking tolerance-based numerical results.
  TEST(Rodin_Solid_ActiveFiberLaw, InitialValues)
  {
    Solid::ActiveFiberLaw::Parameters parameters;
    parameters.initial.stiffness = 4.0;
    parameters.initial.stress = 6.0;

    Solid::ActiveFiberLaw law(parameters);
    const auto state = law.initialState();

    EXPECT_NEAR(state.gamma, 2.0, 1e-14);
    EXPECT_NEAR(state.beta, 3.0, 1e-14);
  }

  // ========================================================================
  // Hooke tests
  // ========================================================================

  /// @brief Verifies isotropic stress for solid hooke by checking tolerance-based numerical results.
  TEST(Rodin_Solid_Hooke, IsotropicStress)
  {
    Solid::Hooke hooke(1.0, 0.5);
    Math::SpatialMatrix<Real> epsilon(2, 2);
    epsilon(0, 0) = 0.1;
    epsilon(0, 1) = 0.05;
    epsilon(1, 0) = 0.05;
    epsilon(1, 1) = 0.2;

    Math::SpatialMatrix<Real> sigma;
    hooke.getStress(sigma, epsilon);

    // sigma = lambda tr(epsilon) I + 2 mu epsilon
    // tr(epsilon) = 0.3
    // sigma = 0.3 * I + 2 * 0.5 * epsilon = 0.3 I + epsilon
    EXPECT_NEAR(sigma(0, 0), 0.3 + 0.1, 1e-14);
    EXPECT_NEAR(sigma(1, 1), 0.3 + 0.2, 1e-14);
    EXPECT_NEAR(sigma(0, 1), 0.05, 1e-14);
    EXPECT_NEAR(sigma(1, 0), 0.05, 1e-14);
  }

  /// @brief Verifies young poisson conversion for solid hooke by checking tolerance-based numerical results.
  TEST(Rodin_Solid_Hooke, YoungPoissonConversion)
  {
    const Real E = 200.0, nu = 0.3;
    auto hooke = Solid::Hooke::YoungPoisson(E, nu);

    const Real expected_lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const Real expected_mu = E / (2.0 * (1.0 + nu));

    EXPECT_NEAR(hooke.getLameFirstParameter(), expected_lambda, 1e-10);
    EXPECT_NEAR(hooke.getShearModulus(), expected_mu, 1e-10);
  }

  // ========================================================================
  // PostProcessing tests
  // ========================================================================

  /// @brief Verifies zero deformation for solid green lagrange strain by checking tolerance-based numerical results.
  TEST(Rodin_Solid_GreenLagrangeStrain, ZeroDeformation)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);

    Solid::GreenLagrangeStrain gl;
    Math::SpatialMatrix<Real> E;
    gl.getGreenLagrangeStrain(E, state);

    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(E(i, j), 0.0, 1e-14);
    }
  }

  /// @brief Verifies pure extension for solid green lagrange strain by checking tolerance-based numerical results.
  TEST(Rodin_Solid_GreenLagrangeStrain, PureExtension)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.1;
    H(0, 1) = 0.0;
    H(1, 0) = 0.0;
    H(1, 1) = 0.0;
    state.setDisplacementGradient(H);

    Solid::GreenLagrangeStrain gl;
    Math::SpatialMatrix<Real> E;
    gl.getGreenLagrangeStrain(E, state);

    // F = [[1.1, 0],[0, 1]], C = [[1.21, 0],[0, 1]]
    // E = 0.5(C - I) = [[0.105, 0],[0, 0]]
    EXPECT_NEAR(E(0, 0), 0.5 * (1.21 - 1.0), 1e-12);
    EXPECT_NEAR(E(1, 1), 0.0, 1e-14);
    EXPECT_NEAR(E(0, 1), 0.0, 1e-14);
    EXPECT_NEAR(E(1, 0), 0.0, 1e-14);
  }

  /// @brief Verifies symmetry check for solid cauchy stress by checking tolerance-based numerical results.
  TEST(Rodin_Solid_CauchyStress, SymmetryCheck)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(1, 0) = -0.02;
    H(1, 1) = 0.15;
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    Solid::CauchyStress<Solid::NeoHookean> cauchy(prepared);
    Math::SpatialMatrix<Real> sigma;
    cauchy.getCauchyStress(sigma);

    // Cauchy stress must be symmetric
    EXPECT_NEAR(sigma(0, 1), sigma(1, 0), 1e-12);
  }

  // ========================================================================
  // ConstitutivePoint tag tests
  // ========================================================================

  /// @brief Verifies typed tag set get for solid constitutive point by checking tolerance-based numerical results, true predicates, false predicates.
  TEST(Rodin_Solid_ConstitutivePoint, TypedTagSetGet)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);

    // Set and get fiber direction
    Math::SpatialVector<Real> fiber(2);
    fiber[0] = 1.0;
    fiber[1] = 0.0;
    cp.set<Solid::Tags::FiberDirection>(fiber);

    EXPECT_TRUE(cp.has<Solid::Tags::FiberDirection>());
    const auto& retrieved = cp.get<Solid::Tags::FiberDirection>();
    EXPECT_NEAR(retrieved[0], 1.0, 1e-14);
    EXPECT_NEAR(retrieved[1], 0.0, 1e-14);

    // Set and get activation
    cp.set<Solid::Tags::Activation>(0.42);
    EXPECT_TRUE(cp.has<Solid::Tags::Activation>());
    EXPECT_NEAR(cp.get<Solid::Tags::Activation>(), 0.42, 1e-14);

    // Tags not set should return false
    EXPECT_FALSE(cp.has<Solid::Tags::SheetDirection>());
    EXPECT_FALSE(cp.has<Solid::Tags::SheetNormalDirection>());
  }

  /// @brief Verifies reference wrapper semantics for solid constitutive point by checking tolerance-based numerical results, false predicates.
  TEST(Rodin_Solid_ConstitutivePoint, ReferenceWrapperSemantics)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);

    Solid::ConstitutivePoint cp(state);

    // Verify kinematic state reference
    EXPECT_NEAR(cp.getKinematicState().getJacobian(), 1.0, 1e-14);

    // Without geometric point, getPoint() should be empty
    EXPECT_FALSE(cp.getPoint().has_value());
  }

  /// @brief Verifies custom tag for solid constitutive point by checking tolerance-based numerical results, true predicates.
  TEST(Rodin_Solid_ConstitutivePoint, CustomTag)
  {
    struct MyCustomTag
    {
      using Type = Real;
    };

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);

    cp.set<MyCustomTag>(3.14);
    EXPECT_TRUE(cp.has<MyCustomTag>());
    EXPECT_NEAR(cp.get<MyCustomTag>(), 3.14, 1e-14);
  }

  // ========================================================================
  // InternalVirtualWorkResidual tests
  // ========================================================================

  /// @brief Verifies zero displacement zero force for solid internal virtual work residual by checking tolerance-based numerical results.
  TEST(Rodin_Solid_InternalVirtualWork_Residual, ZeroDisplacementZeroForce)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TestFunction v(Vh);

    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);

    GridFunction gf(Vh);
    gf.getData().setZero();

    Solid::InternalVirtualWorkResidual force(law, v, gf);

    auto cellIt = mesh.getCell(0);
    force.setPolytope(*cellIt);

    const size_t ndof = 3 * vdim;
    for (size_t i = 0; i < ndof; ++i)
      EXPECT_NEAR(force.integrate(i), 0.0, 1e-14);
  }

  /// @brief Verifies zero displacement zero force SVK for solid internal virtual work residual by checking tolerance-based numerical results.
  TEST(Rodin_Solid_InternalVirtualWork_Residual, ZeroDisplacementZeroForce_SVK)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TestFunction v(Vh);

    GridFunction gf(Vh);
    gf.getData().setZero();

    Solid::SaintVenantKirchhoff law(2.0, 1.0);
    Solid::InternalVirtualWorkResidual force(law, v, gf);

    auto cellIt = mesh.getCell(0);
    force.setPolytope(*cellIt);

    const size_t ndof = 3 * vdim;
    for (size_t i = 0; i < ndof; ++i)
      EXPECT_NEAR(force.integrate(i), 0.0, 1e-14);
  }

  /// @brief Verifies zero displacement zero force mooney rivlin for solid internal virtual work residual by checking tolerance-based numerical results.
  TEST(Rodin_Solid_InternalVirtualWork_Residual, ZeroDisplacementZeroForce_MooneyRivlin)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TestFunction v(Vh);

    GridFunction gf(Vh);
    gf.getData().setZero();

    Solid::MooneyRivlin law(0.5, 0.2, 5.0);
    Solid::InternalVirtualWorkResidual force(law, v, gf);

    auto cellIt = mesh.getCell(0);
    force.setPolytope(*cellIt);

    const size_t ndof = 3 * vdim;
    for (size_t i = 0; i < ndof; ++i)
      EXPECT_NEAR(force.integrate(i), 0.0, 1e-12);
  }

  /// @brief Verifies non zero displacement non zero force for solid internal virtual work residual.
  TEST(Rodin_Solid_InternalVirtualWork_Residual, NonZeroDisplacementNonZeroForce)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TestFunction v(Vh);

    // Set a non-trivial displacement (uniform extension)
    GridFunction gf(Vh);
    gf.getData().setZero();
    for (Index i = 0; i < static_cast<Index>(Vh.getSize()); ++i)
      gf.getData()(i) = 0.01 * static_cast<Real>(i % 3);

    Solid::NeoHookean law(2.0, 1.0);
    Solid::InternalVirtualWorkResidual force(law, v, gf);

    auto cellIt = mesh.getCell(0);
    force.setPolytope(*cellIt);

    // With non-zero displacement, at least some DOFs should have non-zero force
    const size_t ndof = 3 * vdim;
    Real norm = 0;
    for (size_t i = 0; i < ndof; ++i)
      norm += force.integrate(i) * force.integrate(i);
    EXPECT_GT(norm, 1e-14);
  }

  /// @brief Verifies with input function for solid internal virtual work residual by checking true predicates.
  TEST(Rodin_Solid_InternalVirtualWork_Residual, WithInputFunction)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TestFunction v(Vh);

    GridFunction gf(Vh);
    gf.getData().setZero();

    Solid::NeoHookean law(2.0, 1.0);
    Solid::InternalVirtualWorkResidual force(law, v, gf);
    bool inputCalled = false;
    force.setInput([&](Solid::ConstitutivePoint& cp) {
      inputCalled = true;
      Math::SpatialVector<Real> fiber(2);
      fiber[0] = 1.0;
      fiber[1] = 0.0;
      cp.set<Solid::Tags::FiberDirection>(fiber);
    });

    auto cellIt = mesh.getCell(0);
    force.setPolytope(*cellIt);

    // The input should have been called during setPolytope
    EXPECT_TRUE(inputCalled);
  }

  // ========================================================================
  // InternalVirtualWorkTangent tests
  // ========================================================================

  /// @brief Verifies zero displacement symmetry for solid internal virtual work tangent by checking tolerance-based numerical results.
  TEST(Rodin_Solid_InternalVirtualWork_Tangent, ZeroDisplacementSymmetry)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TrialFunction u(Vh);
    TestFunction v(Vh);

    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);

    GridFunction gf(Vh);
    gf.getData().setZero();

    Solid::InternalVirtualWorkTangent tangent(law, u, v, gf);

    auto cellIt = mesh.getCell(0);
    tangent.setPolytope(*cellIt);

    const size_t ndof = 3 * vdim;
    for (size_t i = 0; i < ndof; ++i)
    {
      for (size_t j = 0; j < ndof; ++j)
        EXPECT_NEAR(tangent.integrate(i, j), tangent.integrate(j, i), 1e-12);
    }
  }

  /// @brief Verifies zero displacement symmetry SVK for solid internal virtual work tangent by checking tolerance-based numerical results.
  TEST(Rodin_Solid_InternalVirtualWork_Tangent, ZeroDisplacementSymmetry_SVK)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TrialFunction u(Vh);
    TestFunction v(Vh);

    GridFunction gf(Vh);
    gf.getData().setZero();

    Solid::SaintVenantKirchhoff law(2.0, 1.0);
    Solid::InternalVirtualWorkTangent tangent(law, u, v, gf);

    auto cellIt = mesh.getCell(0);
    tangent.setPolytope(*cellIt);

    const size_t ndof = 3 * vdim;
    for (size_t i = 0; i < ndof; ++i)
    {
      for (size_t j = 0; j < ndof; ++j)
        EXPECT_NEAR(tangent.integrate(i, j), tangent.integrate(j, i), 1e-12);
    }
  }

  /// @brief Verifies zero displacement symmetry mooney rivlin for solid internal virtual work tangent by checking tolerance-based numerical results.
  TEST(Rodin_Solid_InternalVirtualWork_Tangent, ZeroDisplacementSymmetry_MooneyRivlin)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TrialFunction u(Vh);
    TestFunction v(Vh);

    GridFunction gf(Vh);
    gf.getData().setZero();

    Solid::MooneyRivlin law(0.5, 0.2, 5.0);
    Solid::InternalVirtualWorkTangent tangent(law, u, v, gf);

    auto cellIt = mesh.getCell(0);
    tangent.setPolytope(*cellIt);

    const size_t ndof = 3 * vdim;
    for (size_t i = 0; i < ndof; ++i)
    {
      for (size_t j = 0; j < ndof; ++j)
        EXPECT_NEAR(tangent.integrate(i, j), tangent.integrate(j, i), 1e-12);
    }
  }

  /// @brief Verifies non zero displacement symmetry for solid internal virtual work tangent by checking tolerance-based numerical results.
  TEST(Rodin_Solid_InternalVirtualWork_Tangent, NonZeroDisplacementSymmetry)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TrialFunction u(Vh);
    TestFunction v(Vh);

    // Non-trivial linearization point
    GridFunction gf(Vh);
    gf.getData().setZero();
    for (Index i = 0; i < static_cast<Index>(Vh.getSize()); ++i)
      gf.getData()(i) = 0.01 * static_cast<Real>(i % 3);

    Solid::NeoHookean law(2.0, 1.0);
    Solid::InternalVirtualWorkTangent tangent(law, u, v, gf);

    auto cellIt = mesh.getCell(0);
    tangent.setPolytope(*cellIt);

    const size_t ndof = 3 * vdim;
    for (size_t i = 0; i < ndof; ++i)
    {
      for (size_t j = 0; j < ndof; ++j)
        EXPECT_NEAR(tangent.integrate(i, j), tangent.integrate(j, i), 1e-10);
    }
  }

  /// @brief Verifies non zero displacement non trivial for solid internal virtual work tangent.
  TEST(Rodin_Solid_InternalVirtualWork_Tangent, NonZeroDisplacementNonTrivial)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TrialFunction u(Vh);
    TestFunction v(Vh);

    GridFunction gf(Vh);
    gf.getData().setZero();
    for (Index i = 0; i < static_cast<Index>(Vh.getSize()); ++i)
      gf.getData()(i) = 0.01 * static_cast<Real>(i % 3);

    Solid::NeoHookean law(2.0, 1.0);
    Solid::InternalVirtualWorkTangent tangent(law, u, v, gf);

    auto cellIt = mesh.getCell(0);
    tangent.setPolytope(*cellIt);

    // With non-zero displacement, tangent matrix should have non-trivial entries
    const size_t ndof = 3 * vdim;
    Real norm = 0;
    for (size_t i = 0; i < ndof; ++i)
    {
      for (size_t j = 0; j < ndof; ++j)
        norm += tangent.integrate(i, j) * tangent.integrate(i, j);
    }
    EXPECT_GT(norm, 1e-14);
  }

  /// @brief Verifies with input function for solid internal virtual work tangent by checking true predicates.
  TEST(Rodin_Solid_InternalVirtualWork_Tangent, WithInputFunction)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});

    const size_t vdim = 2;
    P1 Vh(mesh, vdim);
    TrialFunction u(Vh);
    TestFunction v(Vh);

    GridFunction gf(Vh);
    gf.getData().setZero();

    Solid::NeoHookean law(2.0, 1.0);
    Solid::InternalVirtualWorkTangent tangent(law, u, v, gf);
    bool inputCalled = false;
    tangent.setInput([&](Solid::ConstitutivePoint& cp) {
      inputCalled = true;
      cp.set<Solid::Tags::Activation>(1.0);
    });

    auto cellIt = mesh.getCell(0);
    tangent.setPolytope(*cellIt);

    EXPECT_TRUE(inputCalled);
  }

  // ========================================================================
  // 3D tests
  // ========================================================================

  /// @brief Verifies tangent finite difference 3 D for solid neo hookean by checking tolerance-based numerical results.
  TEST(Rodin_Solid_NeoHookean, TangentFiniteDifference3D)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);

    Solid::KinematicState state(3);
    Math::SpatialMatrix<Real> H(3, 3);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(0, 2) = -0.02;
    H(1, 0) = -0.02;
    H(1, 1) = 0.15;
    H(1, 2) = 0.03;
    H(2, 0) = 0.01;
    H(2, 1) = -0.01;
    H(2, 2) = 0.08;
    state.setDisplacementGradient(H);

    auto prepared = law;
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> dF(3, 3);
    dF(0, 0) = 0.3;
    dF(0, 1) = -0.2;
    dF(0, 2) = 0.1;
    dF(1, 0) = 0.1;
    dF(1, 1) = 0.4;
    dF(1, 2) = -0.15;
    dF(2, 0) = -0.05;
    dF(2, 1) = 0.2;
    dF(2, 2) = 0.25;

    Math::SpatialMatrix<Real> dP_analytical;
    prepared.getMaterialTangent(dP_analytical, dF);

    // Finite difference
    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(3);
    state_plus.setDisplacementGradient(H_plus);
    auto prepared_plus = law;
    Solid::ConstitutivePoint cp_plus(state_plus);
    prepared_plus.setConstitutivePoint(cp_plus);
    Math::SpatialMatrix<Real> P_plus;
    prepared_plus.getFirstPiolaKirchhoffStress(P_plus);

    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);

    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-5);
    }
  }

  /// @brief Verifies tangent finite difference 3 D for solid saint venant kirchhoff by checking tolerance-based numerical results.
  TEST(Rodin_Solid_SaintVenantKirchhoff, TangentFiniteDifference3D)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::SaintVenantKirchhoff law(lambda, mu);

    Solid::KinematicState state(3);
    Math::SpatialMatrix<Real> H(3, 3);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(0, 2) = -0.02;
    H(1, 0) = -0.02;
    H(1, 1) = 0.15;
    H(1, 2) = 0.03;
    H(2, 0) = 0.01;
    H(2, 1) = -0.01;
    H(2, 2) = 0.08;
    state.setDisplacementGradient(H);

    auto prepared = law;
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> dF(3, 3);
    dF(0, 0) = 0.3;
    dF(0, 1) = -0.2;
    dF(0, 2) = 0.1;
    dF(1, 0) = 0.1;
    dF(1, 1) = 0.4;
    dF(1, 2) = -0.15;
    dF(2, 0) = -0.05;
    dF(2, 1) = 0.2;
    dF(2, 2) = 0.25;

    Math::SpatialMatrix<Real> dP_analytical;
    prepared.getMaterialTangent(dP_analytical, dF);

    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(3);
    state_plus.setDisplacementGradient(H_plus);
    auto prepared_plus = law;
    Solid::ConstitutivePoint cp_plus(state_plus);
    prepared_plus.setConstitutivePoint(cp_plus);
    Math::SpatialMatrix<Real> P_plus;
    prepared_plus.getFirstPiolaKirchhoffStress(P_plus);

    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);

    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-5);
    }
  }

  /// @brief Verifies tangent finite difference 3 D for solid mooney rivlin by checking tolerance-based numerical results.
  TEST(Rodin_Solid_MooneyRivlin, TangentFiniteDifference3D)
  {
    const Real c1 = 0.5, c2 = 0.2, kappa = 5.0;
    Solid::MooneyRivlin law(c1, c2, kappa);

    Solid::KinematicState state(3);
    Math::SpatialMatrix<Real> H(3, 3);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(0, 2) = -0.02;
    H(1, 0) = -0.02;
    H(1, 1) = 0.15;
    H(1, 2) = 0.03;
    H(2, 0) = 0.01;
    H(2, 1) = -0.01;
    H(2, 2) = 0.08;
    state.setDisplacementGradient(H);

    auto prepared = law;
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    Math::SpatialMatrix<Real> dF(3, 3);
    dF(0, 0) = 0.3;
    dF(0, 1) = -0.2;
    dF(0, 2) = 0.1;
    dF(1, 0) = 0.1;
    dF(1, 1) = 0.4;
    dF(1, 2) = -0.15;
    dF(2, 0) = -0.05;
    dF(2, 1) = 0.2;
    dF(2, 2) = 0.25;

    Math::SpatialMatrix<Real> dP_analytical;
    prepared.getMaterialTangent(dP_analytical, dF);

    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(3);
    state_plus.setDisplacementGradient(H_plus);
    auto prepared_plus = law;
    Solid::ConstitutivePoint cp_plus(state_plus);
    prepared_plus.setConstitutivePoint(cp_plus);
    Math::SpatialMatrix<Real> P_plus;
    prepared_plus.getFirstPiolaKirchhoffStress(P_plus);

    Math::SpatialMatrix<Real> P;
    prepared.getFirstPiolaKirchhoffStress(P);

    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-4);
    }
  }

  // ========================================================================
  // Additional law tests
  // ========================================================================

  /// @brief Verifies energy positive under deformation for solid neo hookean.
  TEST(Rodin_Solid_NeoHookean, EnergyPositiveUnderDeformation)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.2;
    H(0, 1) = 0.05;
    H(1, 0) = -0.02;
    H(1, 1) = 0.3;
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    // Energy should be positive for non-zero deformation
    EXPECT_GT(prepared.getStrainEnergyDensity(), 0.0);
  }

  /// @brief Verifies energy positive under deformation for solid saint venant kirchhoff.
  TEST(Rodin_Solid_SaintVenantKirchhoff, EnergyPositiveUnderDeformation)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::SaintVenantKirchhoff law(lambda, mu);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.2;
    H(0, 1) = 0.0;
    H(1, 0) = 0.0;
    H(1, 1) = 0.3;
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    EXPECT_GT(prepared.getStrainEnergyDensity(), 0.0);
  }

  /// @brief Verifies energy positive under deformation for solid mooney rivlin.
  TEST(Rodin_Solid_MooneyRivlin, EnergyPositiveUnderDeformation)
  {
    const Real c1 = 0.5, c2 = 0.2, kappa = 5.0;
    Solid::MooneyRivlin law(c1, c2, kappa);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.2;
    H(0, 1) = 0.0;
    H(1, 0) = 0.0;
    H(1, 1) = 0.3;
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    EXPECT_GT(prepared.getStrainEnergyDensity(), 0.0);
  }

  /// @brief For isotropic material: K = lambda + 2*mu/d (2D), K = lambda + 2*mu/3 (3D).
  TEST(Rodin_Solid_Hooke, BulkModulusConsistency)
  {
    // For isotropic material: K = lambda + 2*mu/d (2D), K = lambda + 2*mu/3 (3D)
    const Real E = 200.0, nu = 0.3;
    auto hooke = Solid::Hooke::YoungPoisson(E, nu);

    const Real lambda = hooke.getLameFirstParameter();
    const Real mu = hooke.getShearModulus();

    // 3D bulk modulus: K = lambda + 2/3 mu
    const Real K3D = lambda + 2.0 / 3.0 * mu;
    const Real K3D_direct = E / (3.0 * (1.0 - 2.0 * nu));
    EXPECT_NEAR(K3D, K3D_direct, 1e-10);
  }

  /// @brief Verifies symmetry with SVK for solid cauchy stress by checking tolerance-based numerical results.
  TEST(Rodin_Solid_CauchyStress, SymmetryWithSVK)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::SaintVenantKirchhoff law(lambda, mu);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(1, 0) = -0.02;
    H(1, 1) = 0.15;
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    Solid::CauchyStress<Solid::SaintVenantKirchhoff> cauchy(prepared);
    Math::SpatialMatrix<Real> sigma;
    cauchy.getCauchyStress(sigma);

    // Cauchy stress must be symmetric
    EXPECT_NEAR(sigma(0, 1), sigma(1, 0), 1e-12);
  }

  /// @brief Verifies symmetry with mooney rivlin for solid cauchy stress by checking tolerance-based numerical results.
  TEST(Rodin_Solid_CauchyStress, SymmetryWithMooneyRivlin)
  {
    const Real c1 = 0.5, c2 = 0.2, kappa = 5.0;
    Solid::MooneyRivlin law(c1, c2, kappa);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(1, 0) = -0.02;
    H(1, 1) = 0.15;
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    Solid::CauchyStress<Solid::MooneyRivlin> cauchy(prepared);
    Math::SpatialMatrix<Real> sigma;
    cauchy.getCauchyStress(sigma);

    EXPECT_NEAR(sigma(0, 1), sigma(1, 0), 1e-12);
  }

  /// @brief Verifies consistency with law for solid first piola kirchhoff stress by checking tolerance-based numerical results.
  TEST(Rodin_Solid_FirstPiolaKirchhoffStress, ConsistencyWithLaw)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);
    auto prepared = law;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0, 0) = 0.1;
    H(0, 1) = 0.05;
    H(1, 0) = -0.02;
    H(1, 1) = 0.15;
    state.setDisplacementGradient(H);
    Solid::ConstitutivePoint cp(state);
    prepared.setConstitutivePoint(cp);

    // Direct from law
    Math::SpatialMatrix<Real> P_direct;
    prepared.getFirstPiolaKirchhoffStress(P_direct);

    // Via Field wrapper
    Solid::FirstPiolaKirchhoffStress<Solid::NeoHookean> pk1(prepared);
    Math::SpatialMatrix<Real> P_wrapper;
    pk1.getFirstPiolaKirchhoffStress(P_wrapper);

    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(P_direct(i, j), P_wrapper(i, j), 1e-14);
    }
  }

  /// @brief Verifies parameter accessors for solid mooney rivlin by checking tolerance-based numerical results.
  TEST(Rodin_Solid_MooneyRivlin, ParameterAccessors)
  {
    const Real c1 = 0.5, c2 = 0.2, kappa = 5.0;
    Solid::MooneyRivlin law(c1, c2, kappa);

    EXPECT_NEAR(law.getMaterialConstantC1(), c1, 1e-14);
    EXPECT_NEAR(law.getMaterialConstantC2(), c2, 1e-14);
    EXPECT_NEAR(law.getBulkModulus(), kappa, 1e-14);
  }

  /// @brief Verifies parameter accessors for solid neo hookean by checking tolerance-based numerical results.
  TEST(Rodin_Solid_NeoHookean, ParameterAccessors)
  {
    Solid::NeoHookean law(3.0, 1.5);
    EXPECT_NEAR(law.getLameFirstParameter(), 3.0, 1e-14);
    EXPECT_NEAR(law.getShearModulus(), 1.5, 1e-14);
  }

  /// @brief Verifies parameter accessors for solid saint venant kirchhoff by checking tolerance-based numerical results.
  TEST(Rodin_Solid_SaintVenantKirchhoff, ParameterAccessors)
  {
    Solid::SaintVenantKirchhoff law(4.0, 2.0);
    EXPECT_NEAR(law.getLameFirstParameter(), 4.0, 1e-14);
    EXPECT_NEAR(law.getShearModulus(), 2.0, 1e-14);
  }

  // ========================================================================
  // ConstitutivePoint tag storage tests
  //
  // The tag storage is a reserved FlatMap keyed by type_index. These cover
  // the round-trip, the open tag set (tags declared outside the library),
  // and the ReservedTags boundary, since ReservedTags is a capacity hint
  // and must not behave as a limit.
  // ========================================================================

  namespace
  {
    /// @brief A tag declared outside the library, of a declared value type.
    struct UserScalarTag
    {
      using Type = Real;
    };

    /// @brief A tag declared outside the library, of a type no library tag uses.
    struct UserMatrixTag
    {
      using Type = Math::SpatialMatrix<Real>;
    };

    /// @brief A family of distinct tag types, used to exceed ReservedTags.
    template <std::size_t I>
    struct OverflowTag
    {
      using Type = Real;
    };

    Solid::KinematicState makeIdentityState()
    {
      Solid::KinematicState state(3);
      Math::SpatialMatrix<Real> H(3, 3);
      H.setZero();
      state.setDisplacementGradient(H);
      return state;
    }
  }

  /// @brief Verifies a declared tag round-trips through the constitutive point.
  TEST(Rodin_Solid_ConstitutivePoint, DeclaredTagRoundTrips)
  {
    auto state = makeIdentityState();
    Solid::ConstitutivePoint cp(state);

    EXPECT_FALSE(cp.has<Solid::Tags::TimeStep>());

    cp.set<Solid::Tags::TimeStep>(0.25);

    EXPECT_TRUE(cp.has<Solid::Tags::TimeStep>());
    EXPECT_NEAR(cp.get<Solid::Tags::TimeStep>(), 0.25, 1e-14);
  }

  /// @brief Verifies setting a tag twice overwrites rather than duplicates.
  TEST(Rodin_Solid_ConstitutivePoint, SetOverwrites)
  {
    auto state = makeIdentityState();
    Solid::ConstitutivePoint cp(state);

    cp.set<Solid::Tags::TimeStep>(0.25);
    cp.set<Solid::Tags::TimeStep>(0.5);

    EXPECT_NEAR(cp.get<Solid::Tags::TimeStep>(), 0.5, 1e-14);
  }

  /// @brief Verifies distinct tags of the same value type do not alias.
  TEST(Rodin_Solid_ConstitutivePoint, DistinctTagsOfSameTypeDoNotAlias)
  {
    auto state = makeIdentityState();
    Solid::ConstitutivePoint cp(state);

    cp.set<Solid::Tags::PreviousActiveGamma>(1.5);
    cp.set<Solid::Tags::PreviousActiveBeta>(2.5);

    EXPECT_NEAR(cp.get<Solid::Tags::PreviousActiveGamma>(), 1.5, 1e-14);
    EXPECT_NEAR(cp.get<Solid::Tags::PreviousActiveBeta>(), 2.5, 1e-14);
  }

  /**
   * @brief Verifies an unset tag reports absent, and a set tag does not make
   *        an unrelated tag appear present.
   */
  TEST(Rodin_Solid_ConstitutivePoint, UnsetTagIsAbsent)
  {
    auto state = makeIdentityState();
    Solid::ConstitutivePoint cp(state);

    cp.set<Solid::Tags::TimeStep>(0.25);

    EXPECT_TRUE(cp.has<Solid::Tags::TimeStep>());
    EXPECT_FALSE(cp.has<Solid::Tags::ElectricalActivation>());
    EXPECT_FALSE(cp.has<Solid::Tags::FiberDirection>());
  }

  /**
   * @brief Verifies the tag set is open: tags declared outside the library,
   *        including one whose value type no library tag uses, round-trip.
   */
  TEST(Rodin_Solid_ConstitutivePoint, UserDefinedTagRoundTrips)
  {
    auto state = makeIdentityState();
    Solid::ConstitutivePoint cp(state);

    Math::SpatialMatrix<Real> m(3, 3);
    m.setIdentity();

    cp.set<UserScalarTag>(0.42);
    cp.set<UserMatrixTag>(m);

    ASSERT_TRUE(cp.has<UserScalarTag>());
    ASSERT_TRUE(cp.has<UserMatrixTag>());
    EXPECT_NEAR(cp.get<UserScalarTag>(), 0.42, 1e-14);
    EXPECT_NEAR(cp.get<UserMatrixTag>()(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(cp.get<UserMatrixTag>()(0, 1), 0.0, 1e-14);
  }

  /**
   * @brief Boundary: ReservedTags is a capacity hint, not a limit. Setting
   *        more tags than are reserved must still round-trip; the storage
   *        reallocates.
   */
  TEST(Rodin_Solid_ConstitutivePoint, ExceedingReservedTagsStillRoundTrips)
  {
    auto state = makeIdentityState();
    Solid::ConstitutivePoint cp(state);

    // One distinct tag type per index, more than ReservedTags of them.
    [&]<std::size_t... I>(std::index_sequence<I...>) {
      (cp.set<OverflowTag<I>>(static_cast<Real>(I)), ...);
    }(std::make_index_sequence<Solid::ConstitutivePoint::ReservedTags + 8>{});

    // gtest macros expand to statements, so check inside a function body.
    auto check = [&]<std::size_t I>() {
      EXPECT_TRUE(cp.has<OverflowTag<I>>());
      EXPECT_NEAR(cp.get<OverflowTag<I>>(), static_cast<Real>(I), 1e-14);
    };

    [&]<std::size_t... I>(std::index_sequence<I...>) {
      (check.template operator()<I>(), ...);
    }(std::make_index_sequence<Solid::ConstitutivePoint::ReservedTags + 8>{});
  }

  /// @brief Verifies a copied constitutive point carries its tags.
  TEST(Rodin_Solid_ConstitutivePoint, CopyCarriesTags)
  {
    auto state = makeIdentityState();
    Solid::ConstitutivePoint cp(state);

    Math::SpatialVector<Real> fiber(3);
    fiber[0] = 1.0;
    fiber[1] = 0.0;
    fiber[2] = 0.0;
    cp.set<Solid::Tags::FiberDirection>(fiber);
    cp.set<Solid::Tags::TimeStep>(0.25);

    Solid::ConstitutivePoint copy(cp);

    ASSERT_TRUE(copy.has<Solid::Tags::FiberDirection>());
    EXPECT_NEAR(copy.get<Solid::Tags::FiberDirection>()[0], 1.0, 1e-14);
    EXPECT_NEAR(copy.get<Solid::Tags::TimeStep>(), 0.25, 1e-14);
  }

  // ========================================================================
  // ActiveContraction midpoint fiber strain tests
  //
  // The series law is evaluated at e_1D^{n+1/2} when Tags::PreviousFiberStrain
  // is supplied (Chapelle, Le Tallec, Moireau, Sorine 2012, eq. 26). These
  // pin the evaluation strain and its chain-rule factor directly.
  // Equal evaluation strains give equal stresses at fixed active fields.
  // ========================================================================

  namespace
  {
    /**
     * @brief Builds a kinematic state with a prescribed fiber strain along x.
     *
     * For F = diag(1+a,1,1) and fiber (1,0,0): I4 = (1+a)^2, so the
     * Green-Lagrange fiber strain is ((1+a)^2 - 1)/2. Invert for a.
     */
    Solid::KinematicState makeUniaxialStateWithFiberStrain(Real strain)
    {
      const Real a = std::sqrt(1.0 + 2.0 * strain) - 1.0;
      Solid::KinematicState state(3);
      Math::SpatialMatrix<Real> H(3, 3);
      H.setZero();
      H(0, 0) = a;
      state.setDisplacementGradient(H);
      return state;
    }

    Solid::ActiveContraction<Solid::NeoHookean, Solid::ActiveFiberLaw>
    makeActiveContractionLaw()
    {
      Solid::NeoHookean passive(0.0, 0.0);
      Solid::ActiveFiberLaw::Parameters p;
      p.stiffness = 80.0;
      p.damping = 0.4;
      p.destructionRate = 0.5;
      p.crossBridgeStiffness = 60.0;
      p.contractility = 40.0;
      Solid::ActiveContraction law(passive, Solid::ActiveFiberLaw(p));
      return law;
    }

    void setDynamicTags(Solid::ConstitutivePoint& cp)
    {
      Math::SpatialVector<Real> fiber(3);
      fiber[0] = 1.0;
      fiber[1] = 0.0;
      fiber[2] = 0.0;
      cp.set<Solid::Tags::FiberDirection>(fiber);
      cp.set<Solid::Tags::TimeStep>(0.01);
      cp.set<Solid::Tags::PreviousActiveExtension>(-0.035);
      cp.set<Solid::Tags::PreviousActiveGamma>(1.7);
      cp.set<Solid::Tags::PreviousActiveBeta>(0.9);
      cp.set<Solid::Tags::ElectricalActivation>(1.1);
      cp.set<Solid::Tags::ActiveExtension>(-0.025);
      cp.set<Solid::Tags::ActiveGamma>(1.8);
      cp.set<Solid::Tags::ActiveBeta>(0.95);
    }
  }

  /**
   * @brief Verifies the series law is evaluated at the midpoint fiber strain:
   *        supplying a previous strain must reproduce the active extension
   *        obtained by evaluating at the midpoint directly.
   */
  TEST(Rodin_Solid_ActiveContraction, MidpointStrainIsUsedForTheSeriesLaw)
  {
    auto law = makeActiveContractionLaw();

    const Real ePrev = -0.02;
    const Real eCurr = 0.06;
    const Real eMid = 0.5 * (eCurr + ePrev);

    // Midpoint via the tag: current strain eCurr, previous strain ePrev.
    auto stateA = makeUniaxialStateWithFiberStrain(eCurr);
    Solid::ConstitutivePoint cpA(stateA);
    setDynamicTags(cpA);
    cpA.set<Solid::Tags::PreviousFiberStrain>(ePrev);
    auto preparedA = law;
    preparedA.setConstitutivePoint(cpA);

    // The same evaluation strain, reached directly with no tag.
    auto stateB = makeUniaxialStateWithFiberStrain(eMid);
    Solid::ConstitutivePoint cpB(stateB);
    setDynamicTags(cpB);
    auto preparedB = law;
    preparedB.setConstitutivePoint(cpB);

    ASSERT_TRUE(preparedA.isDynamic());
    ASSERT_TRUE(preparedB.isDynamic());

    // prepared.getFiberStrain() is the *current* fiber strain, so these differ ...
    EXPECT_NEAR(preparedA.getFiberStrain(), eCurr, 1e-12);
    EXPECT_NEAR(preparedB.getFiberStrain(), eMid, 1e-12);

    // The independent fields are unchanged and the evaluated stresses agree.
    EXPECT_DOUBLE_EQ(preparedA.getActiveExtension(), preparedB.getActiveExtension());
    EXPECT_NEAR(
      preparedA.getActiveResponse().stress, preparedB.getActiveResponse().stress, 1e-12);
    EXPECT_NEAR(preparedA.getActiveResidual().value[0],
      preparedB.getActiveResidual().value[0], 1e-12);
  }

  /**
   * @brief Verifies the partial tangent carries the 1/2 chain-rule factor
   *        when the midpoint strain is used, since d(e^{n+1/2})/d(e^{n+1}) = 1/2.
   */
  TEST(Rodin_Solid_ActiveContraction, MidpointStrainHalvesThePartialTangent)
  {
    auto law = makeActiveContractionLaw();

    const Real ePrev = -0.02;
    const Real eCurr = 0.06;
    const Real eMid = 0.5 * (eCurr + ePrev);

    auto stateA = makeUniaxialStateWithFiberStrain(eCurr);
    Solid::ConstitutivePoint cpA(stateA);
    setDynamicTags(cpA);
    cpA.set<Solid::Tags::PreviousFiberStrain>(ePrev);
    auto preparedA = law;
    preparedA.setConstitutivePoint(cpA);

    auto stateB = makeUniaxialStateWithFiberStrain(eMid);
    Solid::ConstitutivePoint cpB(stateB);
    setDynamicTags(cpB);
    auto preparedB = law;
    preparedB.setConstitutivePoint(cpB);

    EXPECT_GT(std::abs(preparedB.getActiveResponse().tangent), 1e-8);
    EXPECT_NEAR(preparedA.getActiveResponse().tangent,
      0.5 * preparedB.getActiveResponse().tangent, 1e-9);
  }

  /**
   * @brief Regression: omitting Tags::PreviousFiberStrain must reproduce the
   *        pre-midpoint behaviour exactly -- the series law falls back to the
   *        current strain and the tangent carries no 1/2 factor. Existing
   *        callers that never set the tag must be unaffected.
   */
  TEST(Rodin_Solid_ActiveContraction, OmittingPreviousFiberStrainFallsBackToCurrentStrain)
  {
    auto law = makeActiveContractionLaw();

    const Real eCurr = 0.06;

    // No tag: the evaluation strain is the current strain.
    auto stateA = makeUniaxialStateWithFiberStrain(eCurr);
    Solid::ConstitutivePoint cpA(stateA);
    setDynamicTags(cpA);
    auto preparedA = law;
    preparedA.setConstitutivePoint(cpA);

    // Tag equal to the current strain: the midpoint *is* the current strain,
    // so the stress agrees while the partial tangent picks up the 1/2.
    auto stateB = makeUniaxialStateWithFiberStrain(eCurr);
    Solid::ConstitutivePoint cpB(stateB);
    setDynamicTags(cpB);
    cpB.set<Solid::Tags::PreviousFiberStrain>(eCurr);
    auto preparedB = law;
    preparedB.setConstitutivePoint(cpB);

    EXPECT_NEAR(preparedA.getActiveExtension(), preparedB.getActiveExtension(), 1e-12);
    // Guard: without this the tangent comparison would pass trivially if the
    // tangent were near zero.
    EXPECT_GT(std::abs(preparedA.getActiveResponse().tangent), 1e-8);
    EXPECT_NEAR(preparedB.getActiveResponse().tangent,
      0.5 * preparedA.getActiveResponse().tangent, 1e-9);
  }
}

namespace Rodin::Tests::Unit
{
  template <class Scalar>
  class DistinctSpaceElasticity : public ::testing::Test
  {};

  using ElasticityScalars = ::testing::Types<Real, Complex>;
  TYPED_TEST_SUITE(DistinctSpaceElasticity, ElasticityScalars);

  /// @brief Distinct spaces reproduce affine strain work without nodal-layout assumptions.
  TYPED_TEST(DistinctSpaceElasticity, AffineStrainWork)
  {
    using Scalar = TypeParam;
    auto value = [](Real real, Real imaginary = 0) -> Scalar {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(real, imaginary);
      else
        return real;
    };
    auto conjugate = [](Scalar scalar) -> Scalar {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return std::conj(scalar);
      else
        return scalar;
    };
    auto coefficient = [](Real scalar) {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Variational::ComplexFunction(Complex(scalar));
      else
        return Variational::RealFunction(scalar);
    };
    for (const auto geometry : {Geometry::Polytope::Type::Segment,
           Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Quadrilateral,
           Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Hexahedron,
           Geometry::Polytope::Type::Pyramid, Geometry::Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Geometry::Polytope::Traits(geometry).getDimension();
      auto mesh =
        Geometry::LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
      for (size_t d = 1; d <= dimension; ++d)
      {
        mesh.getConnectivity().compute(dimension, d);
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      }
      Math::SpatialMatrix<Scalar> first(dimension, dimension),
        second(dimension, dimension);
      Scalar divergenceFirst = 0, divergenceSecond = 0, strainWork = 0;
      for (size_t row = 0; row < dimension; ++row)
      {
        for (size_t col = 0; col < dimension; ++col)
        {
          first(row, col) =
            value(Real(0.2) * (row + 1) + Real(0.3) * col, Real(0.1) * (row + col + 1));
          second(row, col) = value(Real(0.4) * (col + 1) - Real(0.15) * row,
            Real(-0.2) * (row + 1) + Real(0.05) * col);
        }
        divergenceFirst += first(row, row);
        divergenceSecond += second(row, row);
      }
      for (size_t row = 0; row < dimension; ++row)
      {
        for (size_t col = 0; col < dimension; ++col)
        {
          const Scalar firstStrain = (first(row, col) + first(col, row)) / Real(2);
          const Scalar secondStrain = (second(row, col) + second(col, row)) / Real(2);
          strainWork += firstStrain * conjugate(secondStrain);
        }
      }
      const Scalar expected =
        Real(2) * divergenceFirst * conjugate(divergenceSecond) + Real(6) * strainWork;
      auto check = [&](auto makeSpace) {
        auto trialSpace = makeSpace();
        auto testSpace = makeSpace();
        ASSERT_FALSE(trialSpace == testSpace);
        Variational::TrialFunction trial(trialSpace);
        Variational::TestFunction test(testSpace);
        Variational::GridFunction testField(testSpace);
        trial.getSolution() = Variational::VectorFunction(
          dimension, [first, value](const Geometry::Point& point) {
            Math::SpatialVector<Scalar> result = first * point.getPhysicalCoordinates();
            for (size_t i = 0; i < result.size(); ++i)
              result[i] += value(Real(0.7), Real(-0.3));
            return result;
          });
        testField = Variational::VectorFunction(
          dimension, [second, value](const Geometry::Point& point) {
            Math::SpatialVector<Scalar> result = second * point.getPhysicalCoordinates();
            for (size_t i = 0; i < result.size(); ++i)
              result[i] += value(Real(-0.2), Real(0.4));
            return result;
          });
        Variational::BilinearForm form(trial, test);
        form = Variational::LinearElasticityIntegral(trial, test)(
          coefficient(Real(2)), coefficient(Real(3)));
        for (size_t repeat = 0; repeat < 2; ++repeat)
        {
          form.assemble();
          const Scalar actual = Math::dot(
            form.getOperator() * trial.getSolution().getData(), testField.getData());
          EXPECT_LT(std::abs(actual - expected), Real(1e-10));
        }
        Variational::TestFunction sharedTest(trialSpace);
        Variational::GridFunction sharedField(trialSpace);
        sharedField = Variational::VectorFunction(
          dimension, [second, value](const Geometry::Point& point) {
            Math::SpatialVector<Scalar> result = second * point.getPhysicalCoordinates();
            for (size_t i = 0; i < result.size(); ++i)
              result[i] += value(Real(-0.2), Real(0.4));
            return result;
          });
        Variational::BilinearForm sharedForm(trial, sharedTest);
        sharedForm = Variational::LinearElasticityIntegral(trial, sharedTest)(
          coefficient(Real(2)), coefficient(Real(3)));
        sharedForm.assemble();
        const Scalar sharedAction =
          Math::dot(sharedForm.getOperator() * trial.getSolution().getData(),
            sharedField.getData());
        EXPECT_LT(std::abs(sharedAction - expected), Real(1e-10));
      };
      check(
        [&] { return Variational::P1<Math::SpatialVector<Scalar>>(mesh, dimension); });
      check([&] {
        return Variational::H1<2, Math::SpatialVector<Scalar>>(
          std::integral_constant<size_t, 2>{}, mesh, dimension);
      });
    }
  }
  TEST(Rodin_Solid_ActiveFiberLaw, CoupledResidualAllDerivativeBlocks)
  {
    Solid::ActiveFiberLaw::Parameters parameters;
    parameters.stiffness = 80;
    parameters.damping = 0.4;
    parameters.destructionRate = 0.5;
    parameters.crossBridgeStiffness = 60;
    parameters.contractility = 40;
    const Solid::ActiveFiberLaw law(parameters);
    const Solid::ActiveFiberLaw::State oldState{1.7, 0.9};
    constexpr Real Step = 1e-7;
    for (Real extension : {Real(-0.015), Real(-0.055)})
    {
      const Solid::ActiveFiberLaw::State state{1.8, 0.95};
      Solid::ActiveFiberLaw::Residual reference;
      law.getResidual(
        reference, 0.01, oldState, state, 0.06, -0.035, extension, 1.1, 0.5);
      for (size_t column = 0; column < 4; ++column)
      {
        auto plusState = state;
        auto minusState = state;
        if (column == 1)
        {
          plusState.gamma += Step;
          minusState.gamma -= Step;
        }
        if (column == 2)
        {
          plusState.beta += Step;
          minusState.beta -= Step;
        }
        Solid::ActiveFiberLaw::Residual plus, minus;
        law.getResidual(plus, 0.01, oldState, plusState,
          0.06 + (column == 3 ? 0.5 * Step : 0), -0.035,
          extension + (column == 0 ? Step : 0), 1.1, 0.5);
        law.getResidual(minus, 0.01, oldState, minusState,
          0.06 - (column == 3 ? 0.5 * Step : 0), -0.035,
          extension - (column == 0 ? Step : 0), 1.1, 0.5);
        for (size_t row = 0; row < 3; ++row)
        {
          const Real analytic = column == 3 ? reference.strainDerivative[row]
                                            : reference.derivative(row, column);
          EXPECT_NEAR((plus.value[row] - minus.value[row]) / (2 * Step), analytic,
            1e-7 * (1 + std::abs(analytic)))
            << row << "," << column;
        }
      }
    }
  }

  TEST(Rodin_Solid_ActiveFiberLaw, CoupledResidualPreservesRegularizedUpdate)
  {
    const Solid::ActiveFiberLaw law;
    const Solid::ActiveFiberLaw::State previous{0, 0};
    const auto state = law.update(0.01, previous, 0, 0, -1);
    Solid::ActiveFiberLaw::Residual residual;
    law.getResidual(residual, 0.01, previous, state, 0, 0, 0, -1);
    EXPECT_NEAR(residual.value[1], 0, 1e-30);
    EXPECT_EQ(residual.value[2], 0);
    EXPECT_EQ(residual.derivative(1, 0), 0);
    EXPECT_EQ(state.gamma, law.getMinimumGamma());
    EXPECT_THROW(
      law.getResidual(residual, 0, previous, state, 0, 0, 0, -1), Alert::Exception);
  }

  TEST(Rodin_Solid_ActiveContraction, GlobalDiscontinuousStateNewton)
  {
    using namespace Variational;
    using namespace Geometry;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
    P0 space(mesh);
    GridFunction extension(space), gamma(space), beta(space);
    extension = Real(-0.025);
    gamma = Real(1.8);
    beta = Real(0.95);
    const Solid::ActiveFiberLaw::State previous{1.7, 0.9};
    Solid::ActiveFiberLaw::Parameters parameters;
    parameters.stiffness = 80;
    parameters.damping = 0.4;
    parameters.destructionRate = 0.5;
    parameters.crossBridgeStiffness = 60;
    parameters.contractility = 40;
    const Solid::ActiveFiberLaw law(parameters);
    const auto evaluate = [&](const Point& point) {
      Solid::ActiveFiberLaw::Residual output;
      const Solid::ActiveFiberLaw::State state{gamma(point), beta(point)};
      law.getResidual(output, 0.01, previous, state, 0.06 + 0.02 * point(0), -0.035,
        extension(point), 1.1);
      return output;
    };
    const auto residual = [&](size_t row) {
      return RealFunction(
        [&, row](const Point& point) { return evaluate(point).value[row]; });
    };
    const auto derivative = [&](size_t row, size_t column) {
      return RealFunction([&, row, column](const Point& point) {
        return evaluate(point).derivative(row, column);
      });
    };
    TrialFunction dc(space), dg(space), db(space);
    TestFunction qc(space), qg(space), qb(space);
    const std::array increments{std::ref(dc), std::ref(dg), std::ref(db)};
    const std::array tests{std::ref(qc), std::ref(qg), std::ref(qb)};
    constexpr size_t QuadratureOrder = 4; // Integrates the quadratic balance exactly.
    Problem problem(dc, dg, db, qc, qg, qb);
    auto firstLoad = Integral(residual(0), qc);
    firstLoad.setOrder(QuadratureOrder);
    auto firstBlock = Integral(derivative(0, 0) * dc, qc);
    firstBlock.setOrder(QuadratureOrder);
    problem = firstBlock + firstLoad;
    for (size_t row = 0; row < 3; ++row)
    {
      if (row != 0)
      {
        auto load = Integral(residual(row), tests[row].get());
        load.setOrder(QuadratureOrder);
        problem += load;
      }
      for (size_t column = 0; column < 3; ++column)
      {
        if (row == 0 && column == 0)
          continue;
        auto block =
          Integral(derivative(row, column) * increments[column].get(), tests[row].get());
        block.setOrder(QuadratureOrder);
        problem += block;
      }
    }
    const Index count = extension.getData().size();
    Math::Vector<Real> iterate(3 * count);
    iterate << extension.getData(), gamma.getData(), beta.getData();
    Solver::SparseLU linearSolver(problem);
    Solver::NewtonSolver solver(linearSolver);
    solver.setAbsoluteTolerance(1e-11).setRelativeTolerance(1e-11).setMaxIterations(30);
    solver.setStepPolicy([&](auto& x, const auto& system, const auto&) {
      const Real norm = system.getSolution().norm();
      x += system.getSolution();
      extension.getData() = x.head(count);
      gamma.getData() = x.segment(count, count);
      beta.getData() = x.tail(count);
      return typename decltype(solver)::StepResult{.accepted = true, .stepNorm = norm};
    });
    solver.solve(iterate);
    ASSERT_TRUE(solver.getReport().converged)
      << "iterations=" << solver.getReport().iterations
      << " residual=" << solver.getReport().finalResidual;
    problem.assemble();
    EXPECT_LT(problem.getLinearSystem().getVector().norm(), 1e-10);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      const auto expected = law.update(0.01, previous, -0.035, extension(point), 1.1);
      EXPECT_NEAR(gamma(point), expected.gamma, 1e-10);
      EXPECT_NEAR(beta(point), expected.beta, 1e-10);
    }
    const RealFunction balance([&](const Point& point) {
      const Real c = extension(point);
      const Real midpoint = 0.5 * (c - 0.035);
      const Real d = 1 + 2 * midpoint;
      const Real strain = 0.06 + 0.02 * point(0);
      return d * d * d * (gamma(point) * beta(point) + 0.4 * (c + 0.035) / 0.01) -
        80 * (strain - midpoint) * (1 + 2 * strain);
    });
    const auto& formula = QF::PolytopeQuadratureFormula::get(4, Polytope::Type::Segment);
    Real integratedBalance = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const auto& quadrature = cell->getQuadrature(formula);
      for (size_t sample = 0; sample < formula.getSize(); ++sample)
      {
        const auto& point = quadrature.getPoint(sample);
        integratedBalance +=
          formula.getWeight(sample) * point.getDistortion() * balance(point);
      }
    }
    EXPECT_NEAR(integratedBalance, 0, 1e-10);
  }

  TEST(Rodin_Solid_ActiveContraction, DynamicInputRequiresIndependentFields)
  {
    auto law = makeActiveContractionLaw();
    const auto state = makeUniaxialStateWithFiberStrain(0.06);
    Solid::ConstitutivePoint point(state);
    point.set<Solid::Tags::TimeStep>(0.01);
    auto prepared = law;
    EXPECT_THROW(prepared.setConstitutivePoint(point), Alert::Exception);
    setDynamicTags(point);
    EXPECT_NO_THROW(prepared.setConstitutivePoint(point));
  }

}
