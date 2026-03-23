/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>

#include <Rodin/Solid.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

namespace Rodin::Tests::Unit
{
  // ========================================================================
  // KinematicState tests
  // ========================================================================

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
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(F(i, j), (i == j) ? 1.0 : 0.0, 1e-14);

    const auto& C = state.getRightCauchyGreenTensor();
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(C(i, j), (i == j) ? 1.0 : 0.0, 1e-14);
  }

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
    EXPECT_NEAR(C(0, 0), 1.0,  1e-14);
    EXPECT_NEAR(C(0, 1), 0.1,  1e-14);
    EXPECT_NEAR(C(1, 0), 0.1,  1e-14);
    EXPECT_NEAR(C(1, 1), 1.01, 1e-14);
  }

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
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(product(i, j), (i == j) ? 1.0 : 0.0, 1e-12);
  }

  // ========================================================================
  // Invariants tests
  // ========================================================================

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

  TEST(Rodin_Solid_FiberInvariants, AlignedFiber)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0,0)=0.2; H(0,1)=0.0; H(1,0)=0.0; H(1,1)=0.0;
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

  // ========================================================================
  // NeoHookean tests
  // ========================================================================

  TEST(Rodin_Solid_NeoHookean, ZeroDeformation)
  {
    const Real lambda = 1.0, mu = 0.5;
    Solid::NeoHookean law(lambda, mu);
    Solid::NeoHookean::Cache cache;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);
    law.setCache(cache, state);

    // At zero deformation, P should be zero
    Math::SpatialMatrix<Real> P;
    law.getFirstPiolaKirchhoffStress(P, cache, state);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(P(i, j), 0.0, 1e-14);

    // Energy should also be zero
    EXPECT_NEAR(law.getStrainEnergyDensity(cache, state), 0.0, 1e-14);
  }

  TEST(Rodin_Solid_NeoHookean, StressSymmetryUnderPureStretch)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);
    Solid::NeoHookean::Cache cache;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0,0)=0.1; H(0,1)=0.0; H(1,0)=0.0; H(1,1)=0.1;
    state.setDisplacementGradient(H);
    law.setCache(cache, state);

    Math::SpatialMatrix<Real> P;
    law.getFirstPiolaKirchhoffStress(P, cache, state);

    // Under uniform stretch, P should be diagonal
    EXPECT_NEAR(P(0, 1), 0.0, 1e-13);
    EXPECT_NEAR(P(1, 0), 0.0, 1e-13);
    // And P(0,0) == P(1,1) by symmetry
    EXPECT_NEAR(P(0, 0), P(1, 1), 1e-13);
  }

  TEST(Rodin_Solid_NeoHookean, TangentFiniteDifference)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0,0)=0.1; H(0,1)=0.05; H(1,0)=-0.02; H(1,1)=0.15;
    state.setDisplacementGradient(H);

    Solid::NeoHookean::Cache cache;
    law.setCache(cache, state);

    // Test tangent via finite differences
    Math::SpatialMatrix<Real> dF(2, 2);
    dF(0,0)=0.3; dF(0,1)=-0.2; dF(1,0)=0.1; dF(1,1)=0.4;

    Math::SpatialMatrix<Real> dP_analytical;
    law.getMaterialTangent(dP_analytical, cache, state, dF);

    // Finite difference approximation
    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(2);
    state_plus.setDisplacementGradient(H_plus);
    Solid::NeoHookean::Cache cache_plus;
    law.setCache(cache_plus, state_plus);
    Math::SpatialMatrix<Real> P_plus;
    law.getFirstPiolaKirchhoffStress(P_plus, cache_plus, state_plus);

    Math::SpatialMatrix<Real> P;
    law.getFirstPiolaKirchhoffStress(P, cache, state);

    // SpatialMatrix does not implement direct matrix subtraction.
    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-5);
  }

  // ========================================================================
  // SaintVenantKirchhoff tests
  // ========================================================================

  TEST(Rodin_Solid_SaintVenantKirchhoff, ZeroDeformation)
  {
    const Real lambda = 1.0, mu = 0.5;
    Solid::SaintVenantKirchhoff law(lambda, mu);
    Solid::SaintVenantKirchhoff::Cache cache;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);
    law.setCache(cache, state);

    // At zero deformation, P should be zero (E = 0, S = 0)
    Math::SpatialMatrix<Real> P;
    law.getFirstPiolaKirchhoffStress(P, cache, state);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(P(i, j), 0.0, 1e-14);

    EXPECT_NEAR(law.getStrainEnergyDensity(cache, state), 0.0, 1e-14);
  }

  TEST(Rodin_Solid_SaintVenantKirchhoff, TangentFiniteDifference)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::SaintVenantKirchhoff law(lambda, mu);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0,0)=0.1; H(0,1)=0.05; H(1,0)=-0.02; H(1,1)=0.15;
    state.setDisplacementGradient(H);

    Solid::SaintVenantKirchhoff::Cache cache;
    law.setCache(cache, state);

    Math::SpatialMatrix<Real> dF(2, 2);
    dF(0,0)=0.3; dF(0,1)=-0.2; dF(1,0)=0.1; dF(1,1)=0.4;

    Math::SpatialMatrix<Real> dP_analytical;
    law.getMaterialTangent(dP_analytical, cache, state, dF);

    // Finite difference
    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(2);
    state_plus.setDisplacementGradient(H_plus);
    Solid::SaintVenantKirchhoff::Cache cache_plus;
    law.setCache(cache_plus, state_plus);
    Math::SpatialMatrix<Real> P_plus;
    law.getFirstPiolaKirchhoffStress(P_plus, cache_plus, state_plus);

    Math::SpatialMatrix<Real> P;
    law.getFirstPiolaKirchhoffStress(P, cache, state);

    // SpatialMatrix does not implement direct matrix subtraction.
    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-5);
  }

  TEST(Rodin_Solid_SaintVenantKirchhoff, ReducesToLinearElasticity)
  {
    // For small deformations, SVK should behave like linear elasticity
    const Real lambda = 1.0, mu = 0.5;
    Solid::SaintVenantKirchhoff law(lambda, mu);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    const Real eps = 1e-6;
    H(0,0)=eps; H(0,1)=0.5*eps; H(1,0)=0.5*eps; H(1,1)=2.0*eps;
    state.setDisplacementGradient(H);

    Solid::SaintVenantKirchhoff::Cache cache;
    law.setCache(cache, state);

    Math::SpatialMatrix<Real> P;
    law.getFirstPiolaKirchhoffStress(P, cache, state);

    // For infinitesimal strain, P ~ sigma = lambda tr(epsilon) I + 2 mu epsilon
    // where epsilon = 0.5 (H + H^T)
    Math::SpatialMatrix<Real> epsilon = 0.5 * (H + H.transpose());
    Math::SpatialMatrix<Real> I2(2,2); I2.setIdentity();
    Math::SpatialMatrix<Real> sigma = lambda * epsilon.trace() * I2 + 2.0 * mu * epsilon;

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(P(i, j), sigma(i, j), 1e-10);
  }

  // ========================================================================
  // MooneyRivlin tests
  // ========================================================================

  TEST(Rodin_Solid_MooneyRivlin, ZeroDeformation)
  {
    const Real c1 = 0.5, c2 = 0.1, kappa = 10.0;
    Solid::MooneyRivlin law(c1, c2, kappa);
    Solid::MooneyRivlin::Cache cache;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H.setZero();
    state.setDisplacementGradient(H);
    law.setCache(cache, state);

    // At zero deformation, P should be zero
    Math::SpatialMatrix<Real> P;
    law.getFirstPiolaKirchhoffStress(P, cache, state);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(P(i, j), 0.0, 1e-12);
  }

  TEST(Rodin_Solid_MooneyRivlin, TangentFiniteDifference)
  {
    const Real c1 = 0.5, c2 = 0.2, kappa = 5.0;
    Solid::MooneyRivlin law(c1, c2, kappa);

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0,0)=0.1; H(0,1)=0.05; H(1,0)=-0.02; H(1,1)=0.15;
    state.setDisplacementGradient(H);

    Solid::MooneyRivlin::Cache cache;
    law.setCache(cache, state);

    Math::SpatialMatrix<Real> dF(2, 2);
    dF(0,0)=0.3; dF(0,1)=-0.2; dF(1,0)=0.1; dF(1,1)=0.4;

    Math::SpatialMatrix<Real> dP_analytical;
    law.getMaterialTangent(dP_analytical, cache, state, dF);

    // Finite difference
    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(2);
    state_plus.setDisplacementGradient(H_plus);
    Solid::MooneyRivlin::Cache cache_plus;
    law.setCache(cache_plus, state_plus);
    Math::SpatialMatrix<Real> P_plus;
    law.getFirstPiolaKirchhoffStress(P_plus, cache_plus, state_plus);

    Math::SpatialMatrix<Real> P;
    law.getFirstPiolaKirchhoffStress(P, cache, state);

    // SpatialMatrix does not implement direct matrix subtraction.
    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-4);
  }

  // ========================================================================
  // Hooke tests
  // ========================================================================

  TEST(Rodin_Solid_Hooke, IsotropicStress)
  {
    Solid::Hooke hooke(1.0, 0.5);
    Math::SpatialMatrix<Real> epsilon(2, 2);
    epsilon(0,0)=0.1; epsilon(0,1)=0.05; epsilon(1,0)=0.05; epsilon(1,1)=0.2;

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
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(E(i, j), 0.0, 1e-14);
  }

  TEST(Rodin_Solid_GreenLagrangeStrain, PureExtension)
  {
    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0,0)=0.1; H(0,1)=0.0; H(1,0)=0.0; H(1,1)=0.0;
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

  TEST(Rodin_Solid_CauchyStress, SymmetryCheck)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);
    Solid::NeoHookean::Cache cache;

    Solid::KinematicState state(2);
    Math::SpatialMatrix<Real> H(2, 2);
    H(0,0)=0.1; H(0,1)=0.05; H(1,0)=-0.02; H(1,1)=0.15;
    state.setDisplacementGradient(H);
    law.setCache(cache, state);

    Solid::CauchyStress<Solid::NeoHookean> cauchy(law);
    Math::SpatialMatrix<Real> sigma;
    cauchy.getCauchyStress(sigma, cache, state);

    // Cauchy stress must be symmetric
    EXPECT_NEAR(sigma(0, 1), sigma(1, 0), 1e-12);
  }

  // ========================================================================
  // InternalForce and MaterialTangent tests
  // ========================================================================

  TEST(Rodin_Solid_InternalForce, ZeroDisplacementZeroForce)
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

    Solid::InternalForce force(law, v);

    Math::Vector<Real> data = Math::Vector<Real>::Zero(Vh.getSize());
    force.setLinearizationPoint(data);

    auto cellIt = mesh.getCell(0);
    force.setPolytope(*cellIt);

    const size_t ndof = 3 * vdim;
    for (size_t i = 0; i < ndof; ++i)
      EXPECT_NEAR(force.integrate(i), 0.0, 1e-14);
  }

  TEST(Rodin_Solid_MaterialTangent, ZeroDisplacementSymmetry)
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

    Solid::MaterialTangent tangent(law, u, v);

    Math::Vector<Real> data = Math::Vector<Real>::Zero(Vh.getSize());
    tangent.setLinearizationPoint(data);

    auto cellIt = mesh.getCell(0);
    tangent.setPolytope(*cellIt);

    const size_t ndof = 3 * vdim;
    for (size_t i = 0; i < ndof; ++i)
      for (size_t j = 0; j < ndof; ++j)
        EXPECT_NEAR(tangent.integrate(i, j), tangent.integrate(j, i), 1e-12);
  }

  // ========================================================================
  // 3D tests
  // ========================================================================

  TEST(Rodin_Solid_NeoHookean, TangentFiniteDifference3D)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::NeoHookean law(lambda, mu);

    Solid::KinematicState state(3);
    Math::SpatialMatrix<Real> H(3, 3);
    H(0,0)=0.1; H(0,1)=0.05; H(0,2)=-0.02;
    H(1,0)=-0.02; H(1,1)=0.15; H(1,2)=0.03;
    H(2,0)=0.01; H(2,1)=-0.01; H(2,2)=0.08;
    state.setDisplacementGradient(H);

    Solid::NeoHookean::Cache cache;
    law.setCache(cache, state);

    Math::SpatialMatrix<Real> dF(3, 3);
    dF(0,0)=0.3; dF(0,1)=-0.2; dF(0,2)=0.1;
    dF(1,0)=0.1; dF(1,1)=0.4; dF(1,2)=-0.15;
    dF(2,0)=-0.05; dF(2,1)=0.2; dF(2,2)=0.25;

    Math::SpatialMatrix<Real> dP_analytical;
    law.getMaterialTangent(dP_analytical, cache, state, dF);

    // Finite difference
    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(3);
    state_plus.setDisplacementGradient(H_plus);
    Solid::NeoHookean::Cache cache_plus;
    law.setCache(cache_plus, state_plus);
    Math::SpatialMatrix<Real> P_plus;
    law.getFirstPiolaKirchhoffStress(P_plus, cache_plus, state_plus);

    Math::SpatialMatrix<Real> P;
    law.getFirstPiolaKirchhoffStress(P, cache, state);

    // SpatialMatrix does not implement direct matrix subtraction.
    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-5);
  }

  TEST(Rodin_Solid_SaintVenantKirchhoff, TangentFiniteDifference3D)
  {
    const Real lambda = 2.0, mu = 1.0;
    Solid::SaintVenantKirchhoff law(lambda, mu);

    Solid::KinematicState state(3);
    Math::SpatialMatrix<Real> H(3, 3);
    H(0,0)=0.1; H(0,1)=0.05; H(0,2)=-0.02;
    H(1,0)=-0.02; H(1,1)=0.15; H(1,2)=0.03;
    H(2,0)=0.01; H(2,1)=-0.01; H(2,2)=0.08;
    state.setDisplacementGradient(H);

    Solid::SaintVenantKirchhoff::Cache cache;
    law.setCache(cache, state);

    Math::SpatialMatrix<Real> dF(3, 3);
    dF(0,0)=0.3; dF(0,1)=-0.2; dF(0,2)=0.1;
    dF(1,0)=0.1; dF(1,1)=0.4; dF(1,2)=-0.15;
    dF(2,0)=-0.05; dF(2,1)=0.2; dF(2,2)=0.25;

    Math::SpatialMatrix<Real> dP_analytical;
    law.getMaterialTangent(dP_analytical, cache, state, dF);

    const Real eps = 1e-7;
    Math::SpatialMatrix<Real> H_plus = H + eps * dF;
    Solid::KinematicState state_plus(3);
    state_plus.setDisplacementGradient(H_plus);
    Solid::SaintVenantKirchhoff::Cache cache_plus;
    law.setCache(cache_plus, state_plus);
    Math::SpatialMatrix<Real> P_plus;
    law.getFirstPiolaKirchhoffStress(P_plus, cache_plus, state_plus);

    Math::SpatialMatrix<Real> P;
    law.getFirstPiolaKirchhoffStress(P, cache, state);

    // SpatialMatrix does not implement direct matrix subtraction.
    Math::SpatialMatrix<Real> dP_fd = (1.0 / eps) * P_plus + (-1.0 / eps) * P;

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(dP_analytical(i, j), dP_fd(i, j), 1e-5);
  }
}
