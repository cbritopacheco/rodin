/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <gtest/gtest.h>

#include <Rodin/Types.h>
#include <Rodin/Math/SpatialVector.h>
#include <Rodin/Math/SpatialMatrix.h>

#include <Rodin/Fluid/Local/FlowPoint.h>
#include <Rodin/Fluid/Local/Tags.h>
#include <Rodin/Fluid/Local/Input.h>
#include <Rodin/Fluid/Fields/StrainRate.h>
#include <Rodin/Fluid/Fields/ShearRate.h>
#include <Rodin/Fluid/Fields/Vorticity.h>
#include <Rodin/Fluid/Fields/CauchyStress.h>
#include <Rodin/Fluid/Constitutive/Newtonian.h>
#include <Rodin/Fluid/Constitutive/PowerLaw.h>
#include <Rodin/Fluid/Constitutive/CarreauYasuda.h>
#include <Rodin/Fluid/Constitutive/Bingham.h>
#include <Rodin/Fluid/Constitutive/HerschelBulkley.h>

using namespace Rodin;
using namespace Rodin::Fluid;

namespace Rodin::Tests::Unit
{
  // =========================================================================
  // Helper: build FlowPoint from velocity gradient (velocity set to zero)
  // =========================================================================
  static FlowPoint makeFlowPoint(const Math::SpatialMatrix<Real>& gradU)
  {
    const size_t d = static_cast<size_t>(gradU.rows());
    Math::SpatialVector<Real> vel(d);
    vel.setZero();
    return FlowPoint(vel, gradU);
  }

  // =========================================================================
  // FlowPoint
  // =========================================================================

  TEST(Rodin_Fluid_FlowPoint, ConstructWithVelocityAndGradient2D)
  {
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    Math::SpatialVector<Real> vel{{1.0, 2.0}};
    FlowPoint fp(vel, gradU);

    EXPECT_EQ(fp.getDimension(), 2u);
    EXPECT_NEAR(fp.getVelocity()(0), 1.0, 1e-15);
    EXPECT_NEAR(fp.getVelocity()(1), 2.0, 1e-15);
    EXPECT_FALSE(fp.hasPressure());
    EXPECT_FALSE(fp.hasPressureGradient());
    EXPECT_FALSE(fp.getPoint().has_value());
  }

  TEST(Rodin_Fluid_FlowPoint, SetAndGetPressure)
  {
    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(2, 2));
    EXPECT_FALSE(fp.hasPressure());
    fp.setPressure(3.14);
    EXPECT_TRUE(fp.hasPressure());
    EXPECT_NEAR(fp.getPressure().value(), 3.14, 1e-15);
  }

  TEST(Rodin_Fluid_FlowPoint, SetAndGetPressureGradient2D)
  {
    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(2, 2));
    Math::SpatialVector<Real> gradP{{0.5, -0.5}};
    fp.setPressureGradient(gradP);
    EXPECT_TRUE(fp.hasPressureGradient());
    EXPECT_NEAR(fp.getPressureGradient().value()(0),  0.5, 1e-15);
    EXPECT_NEAR(fp.getPressureGradient().value()(1), -0.5, 1e-15);
  }

  TEST(Rodin_Fluid_FlowPoint, SetAndGetAuxTag)
  {
    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(2, 2));
    EXPECT_FALSE(fp.has<Tags::Density>());
    fp.set<Tags::Density>(1050.0);
    EXPECT_TRUE(fp.has<Tags::Density>());
    EXPECT_NEAR(fp.get<Tags::Density>(), 1050.0, 1e-15);
  }

  TEST(Rodin_Fluid_FlowPoint, MultipleTagsIndependent)
  {
    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(3, 3));
    fp.set<Tags::Density>(1000.0);
    fp.set<Tags::DynamicViscosity>(0.001);
    fp.set<Tags::YieldStress>(5.0);

    EXPECT_NEAR(fp.get<Tags::Density>(), 1000.0, 1e-15);
    EXPECT_NEAR(fp.get<Tags::DynamicViscosity>(), 0.001, 1e-15);
    EXPECT_NEAR(fp.get<Tags::YieldStress>(), 5.0, 1e-15);
    EXPECT_FALSE(fp.has<Tags::WallDistance>());
  }

  TEST(Rodin_Fluid_FlowPoint, GetDimension3D)
  {
    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(3, 3));
    EXPECT_EQ(fp.getDimension(), 3u);
  }

  // =========================================================================
  // StrainRate
  // =========================================================================

  TEST(Rodin_Fluid_StrainRate, ZeroGradient)
  {
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    auto fp = makeFlowPoint(gradU);

    StrainRate sr;
    Math::SpatialMatrix<Real> D(2, 2);
    sr.getStrainRate(D, fp);

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(D(i, j), 0.0, 1e-15);
  }

  TEST(Rodin_Fluid_StrainRate, PureShear2D)
  {
    // gradU = [[0, gamma], [0, 0]]  =>  D = [[0, gamma/2], [gamma/2, 0]]
    const Real gamma = 4.0;
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = gamma;

    auto fp = makeFlowPoint(gradU);
    StrainRate sr;
    Math::SpatialMatrix<Real> D(2, 2);
    sr.getStrainRate(D, fp);

    EXPECT_NEAR(D(0, 0), 0.0,       1e-15);
    EXPECT_NEAR(D(1, 1), 0.0,       1e-15);
    EXPECT_NEAR(D(0, 1), gamma / 2, 1e-15);
    EXPECT_NEAR(D(1, 0), gamma / 2, 1e-15);
  }

  TEST(Rodin_Fluid_StrainRate, Symmetry)
  {
    // Arbitrary asymmetric gradient – D must be symmetric
    Math::SpatialMatrix<Real> gradU(3, 3);
    gradU(0, 0) = 1.0; gradU(0, 1) = 2.0; gradU(0, 2) = 3.0;
    gradU(1, 0) = 4.0; gradU(1, 1) = 5.0; gradU(1, 2) = 6.0;
    gradU(2, 0) = 7.0; gradU(2, 1) = 8.0; gradU(2, 2) = 9.0;

    StrainRate sr;
    Math::SpatialMatrix<Real> D(3, 3);
    sr.getStrainRate(D, gradU);

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(D(i, j), D(j, i), 1e-14);
  }

  TEST(Rodin_Fluid_StrainRate, DiagonalGradient)
  {
    // gradU = diag(1, 2) => D = gradU (already symmetric)
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 0) = 1.0;
    gradU(1, 1) = 2.0;
    auto fp = makeFlowPoint(gradU);

    StrainRate sr;
    Math::SpatialMatrix<Real> D(2, 2);
    sr.getStrainRate(D, fp);

    EXPECT_NEAR(D(0, 0), 1.0, 1e-15);
    EXPECT_NEAR(D(1, 1), 2.0, 1e-15);
    EXPECT_NEAR(D(0, 1), 0.0, 1e-15);
    EXPECT_NEAR(D(1, 0), 0.0, 1e-15);
  }

  // =========================================================================
  // ShearRate
  // =========================================================================

  TEST(Rodin_Fluid_ShearRate, ZeroGradient)
  {
    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(2, 2));
    ShearRate shearRate;
    EXPECT_NEAR(shearRate.getShearRate(fp), 0.0, 1e-15);
  }

  TEST(Rodin_Fluid_ShearRate, PureShear2DAnalytical)
  {
    // gradU = [[0, gamma], [0, 0]]  =>  D = [[0, gamma/2], [gamma/2, 0]]
    // ||D||_F^2 = 2*(gamma/2)^2 = gamma^2/2
    // dot(D,D) = gamma^2/2  =>  gamma_dot = sqrt(2 * gamma^2/2) = |gamma|
    const Real gamma = 3.0;
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = gamma;
    auto fp = makeFlowPoint(gradU);

    ShearRate shearRate;
    EXPECT_NEAR(shearRate.getShearRate(fp), gamma, 1e-14);
  }

  TEST(Rodin_Fluid_ShearRate, PureExtension2D)
  {
    // gradU = [[e, 0], [0, -e]]  (incompressible)
    // D = gradU => dot(D,D) = e^2 + e^2 = 2*e^2 => gamma_dot = 2|e|
    const Real e = 1.5;
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 0) =  e;
    gradU(1, 1) = -e;
    auto fp = makeFlowPoint(gradU);

    ShearRate shearRate;
    EXPECT_NEAR(shearRate.getShearRate(fp), 2.0 * e, 1e-14);
  }

  // =========================================================================
  // Vorticity
  // =========================================================================

  TEST(Rodin_Fluid_Vorticity, ZeroGradient2D)
  {
    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(2, 2));
    Vorticity vort;

    Math::SpatialMatrix<Real> W(2, 2);
    vort.getVorticityTensor(W, fp);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(W(i, j), 0.0, 1e-15);

    Math::SpatialVector<Real> omega;
    vort.getVorticityVector(omega, fp);
    EXPECT_EQ(omega.size(), 1);
    EXPECT_NEAR(omega(0), 0.0, 1e-15);
  }

  TEST(Rodin_Fluid_Vorticity, PureRotation2D)
  {
    // gradU = [[0, -omega], [omega, 0]]  =>  W = gradU, D = 0
    const Real w = 2.0;
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = -w;
    gradU(1, 0) =  w;
    auto fp = makeFlowPoint(gradU);

    Vorticity vort;
    Math::SpatialMatrix<Real> W(2, 2);
    vort.getVorticityTensor(W, fp);

    EXPECT_NEAR(W(0, 0),  0.0, 1e-15);
    EXPECT_NEAR(W(1, 1),  0.0, 1e-15);
    EXPECT_NEAR(W(0, 1), -w,   1e-15);
    EXPECT_NEAR(W(1, 0),  w,   1e-15);

    Math::SpatialVector<Real> omega;
    vort.getVorticityVector(omega, fp);
    EXPECT_EQ(omega.size(), 1);
    // omega_z = du_y/dx - du_x/dy = gradU(1,0) - gradU(0,1) = w - (-w) = 2w
    EXPECT_NEAR(omega(0), 2.0 * w, 1e-14);
  }

  TEST(Rodin_Fluid_Vorticity, Antisymmetry2D)
  {
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU(0, 0) = 1.0; gradU(0, 1) = 2.0;
    gradU(1, 0) = 3.0; gradU(1, 1) = 4.0;
    auto fp = makeFlowPoint(gradU);

    Vorticity vort;
    Math::SpatialMatrix<Real> W(2, 2);
    vort.getVorticityTensor(W, fp);

    // W must be anti-symmetric
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(W(i, j), -W(j, i), 1e-14);
  }

  TEST(Rodin_Fluid_Vorticity, VorticityVector3D)
  {
    // gradU such that curl(u) = (1, 2, 3)
    // curl_x = gradU(2,1) - gradU(1,2) = 1
    // curl_y = gradU(0,2) - gradU(2,0) = 2
    // curl_z = gradU(1,0) - gradU(0,1) = 3
    Math::SpatialMatrix<Real> gradU(3, 3);
    gradU.setZero();
    gradU(2, 1) = 0.5;  gradU(1, 2) = -0.5;  // diff = 1
    gradU(0, 2) = 1.0;  gradU(2, 0) = -1.0;  // diff = 2
    gradU(1, 0) = 1.5;  gradU(0, 1) = -1.5;  // diff = 3
    auto fp = makeFlowPoint(gradU);

    Vorticity vort;
    Math::SpatialVector<Real> omega;
    vort.getVorticityVector(omega, fp);

    ASSERT_EQ(omega.size(), 3);
    EXPECT_NEAR(omega(0), 1.0, 1e-14);
    EXPECT_NEAR(omega(1), 2.0, 1e-14);
    EXPECT_NEAR(omega(2), 3.0, 1e-14);
  }

  // =========================================================================
  // Newtonian rheology
  // =========================================================================

  TEST(Rodin_Fluid_Newtonian, ZeroShear)
  {
    Newtonian law(1.0);
    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(2, 2));

    Newtonian::Cache cache;
    law.setCache(cache, fp);

    Math::SpatialMatrix<Real> tau(2, 2);
    law.getDeviatoricStress(tau, cache, fp);

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(tau(i, j), 0.0, 1e-15);
  }

  TEST(Rodin_Fluid_Newtonian, PureShear2D)
  {
    // tau = 2 mu D, D(0,1)=D(1,0)=gamma/2  =>  tau(0,1)=tau(1,0)=mu*gamma
    const Real mu = 2.5;
    const Real gamma = 4.0;
    Newtonian law(mu);

    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = gamma;
    auto fp = makeFlowPoint(gradU);

    Newtonian::Cache cache;
    law.setCache(cache, fp);

    Math::SpatialMatrix<Real> tau(2, 2);
    law.getDeviatoricStress(tau, cache, fp);

    EXPECT_NEAR(tau(0, 0), 0.0,       1e-14);
    EXPECT_NEAR(tau(1, 1), 0.0,       1e-14);
    EXPECT_NEAR(tau(0, 1), mu * gamma, 1e-14);
    EXPECT_NEAR(tau(1, 0), mu * gamma, 1e-14);
  }

  TEST(Rodin_Fluid_Newtonian, TangentConsistency2D)
  {
    // For Newtonian: dTau = 2 mu dD, compare tangent vs finite difference
    const Real mu = 1.0;
    const Real gamma = 2.0;
    Newtonian law(mu);

    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = gamma;
    auto fp = makeFlowPoint(gradU);

    Newtonian::Cache cache;
    law.setCache(cache, fp);

    // Perturbation dGradU
    Math::SpatialMatrix<Real> dGradU(2, 2);
    dGradU.setZero();
    dGradU(0, 0) = 0.1;
    dGradU(1, 0) = 0.2;

    Math::SpatialMatrix<Real> dtau(2, 2);
    law.getTangent(dtau, cache, fp, dGradU);

    // Expected: 2 mu * 0.5*(dGradU + dGradU^T)
    const Math::SpatialMatrix<Real> dD = 0.5 * (dGradU + dGradU.transpose());
    const Math::SpatialMatrix<Real> expected = 2.0 * mu * dD;

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dtau(i, j), expected(i, j), 1e-14);
  }

  TEST(Rodin_Fluid_Newtonian, GetDynamicViscosity)
  {
    Newtonian law(0.003);
    EXPECT_NEAR(law.getDynamicViscosity(), 0.003, 1e-15);
  }

  // =========================================================================
  // PowerLaw rheology
  // =========================================================================

  TEST(Rodin_Fluid_PowerLaw, NewtonianLimitN1)
  {
    // n=1 => mu_eff = K, independent of shear rate => equivalent to Newtonian
    const Real K = 2.0;
    const Real gamma = 3.0;
    PowerLaw law(K, 1.0);

    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = gamma;
    auto fp = makeFlowPoint(gradU);

    PowerLaw::Cache cache;
    law.setCache(cache, fp);
    EXPECT_NEAR(cache.effectiveViscosity, K, 1e-10);

    Math::SpatialMatrix<Real> tau(2, 2);
    law.getDeviatoricStress(tau, cache, fp);
    // tau(0,1) = 2 * K * gamma/2 = K * gamma
    EXPECT_NEAR(tau(0, 1), K * gamma, 1e-10);
    EXPECT_NEAR(tau(1, 0), K * gamma, 1e-10);
  }

  TEST(Rodin_Fluid_PowerLaw, ShearThinning)
  {
    // n < 1 => shear-thinning: higher shear rate => lower effective viscosity
    const Real K  = 1.0;
    const Real n  = 0.5;
    PowerLaw law(K, n);

    auto makeGradU = [](Real gamma) {
      Math::SpatialMatrix<Real> g(2, 2);
      g.setZero();
      g(0, 1) = gamma;
      return g;
    };

    auto fp1 = makeFlowPoint(makeGradU(1.0));
    auto fp2 = makeFlowPoint(makeGradU(10.0));

    PowerLaw::Cache c1, c2;
    law.setCache(c1, fp1);
    law.setCache(c2, fp2);

    EXPECT_GT(c1.effectiveViscosity, c2.effectiveViscosity);
  }

  TEST(Rodin_Fluid_PowerLaw, ShearThickening)
  {
    // n > 1 => shear-thickening: higher shear rate => higher effective viscosity
    const Real K = 1.0;
    const Real n = 1.5;
    PowerLaw law(K, n);

    auto makeGradU = [](Real gamma) {
      Math::SpatialMatrix<Real> g(2, 2);
      g.setZero();
      g(0, 1) = gamma;
      return g;
    };

    auto fp1 = makeFlowPoint(makeGradU(1.0));
    auto fp2 = makeFlowPoint(makeGradU(10.0));

    PowerLaw::Cache c1, c2;
    law.setCache(c1, fp1);
    law.setCache(c2, fp2);

    EXPECT_LT(c1.effectiveViscosity, c2.effectiveViscosity);
  }

  TEST(Rodin_Fluid_PowerLaw, TangentConsistencyFiniteDiff)
  {
    // Check tangent by finite difference
    const Real K = 0.5;
    const Real n = 0.6;
    const Real eps = 1e-12;
    PowerLaw law(K, n, eps);

    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = 2.0;
    gradU(1, 0) = 0.5;
    auto fp = makeFlowPoint(gradU);

    PowerLaw::Cache cache;
    law.setCache(cache, fp);

    Math::SpatialMatrix<Real> dGradU(2, 2);
    dGradU.setZero();
    dGradU(0, 0) = 0.01;
    dGradU(1, 1) = -0.01;

    Math::SpatialMatrix<Real> dtau_analytic(2, 2);
    law.getTangent(dtau_analytic, cache, fp, dGradU);

    // Finite difference
    const Real h = 1e-6;
    auto fpPlus = makeFlowPoint(gradU + h * dGradU);
    auto fpMinus = makeFlowPoint(gradU - h * dGradU);

    PowerLaw::Cache cPlus, cMinus;
    law.setCache(cPlus, fpPlus);
    law.setCache(cMinus, fpMinus);

    Math::SpatialMatrix<Real> tauPlus(2, 2), tauMinus(2, 2);
    law.getDeviatoricStress(tauPlus,  cPlus,  fpPlus);
    law.getDeviatoricStress(tauMinus, cMinus, fpMinus);

    const Math::SpatialMatrix<Real> dtau_fd = (tauPlus - tauMinus) / (2.0 * h);

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dtau_analytic(i, j), dtau_fd(i, j), 1e-8);
  }

  // =========================================================================
  // CarreauYasuda rheology
  // =========================================================================

  TEST(Rodin_Fluid_CarreauYasuda, LowShearLimitApproachesMu0)
  {
    // At zero shear rate the viscosity should be close to mu_0
    const Real mu0   = 0.056;
    const Real muInf = 0.0035;
    const Real lambda = 3.313;
    const Real n = 0.3568;
    const Real a = 2.0;
    CarreauYasuda law(mu0, muInf, lambda, n, a);

    // Very small shear gradient
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = 1e-8;
    auto fp = makeFlowPoint(gradU);

    CarreauYasuda::Cache cache;
    law.setCache(cache, fp);

    EXPECT_NEAR(cache.effectiveViscosity, mu0, 1e-5);
  }

  TEST(Rodin_Fluid_CarreauYasuda, HighShearLimitApproachesMuInf)
  {
    const Real mu0   = 0.056;
    const Real muInf = 0.0035;
    const Real lambda = 3.313;
    const Real n = 0.3568;
    const Real a = 2.0;
    CarreauYasuda law(mu0, muInf, lambda, n, a);

    // Very large shear gradient
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = 1e6;
    auto fp = makeFlowPoint(gradU);

    CarreauYasuda::Cache cache;
    law.setCache(cache, fp);

    EXPECT_NEAR(cache.effectiveViscosity, muInf, 1e-4);
  }

  TEST(Rodin_Fluid_CarreauYasuda, TangentConsistencyFiniteDiff)
  {
    CarreauYasuda law(0.056, 0.0035, 3.313, 0.3568, 2.0);

    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = 0.5;
    auto fp = makeFlowPoint(gradU);

    CarreauYasuda::Cache cache;
    law.setCache(cache, fp);

    Math::SpatialMatrix<Real> dGradU(2, 2);
    dGradU.setZero();
    dGradU(0, 0) = 0.01;
    dGradU(1, 1) = -0.01;

    Math::SpatialMatrix<Real> dtau_analytic(2, 2);
    law.getTangent(dtau_analytic, cache, fp, dGradU);

    const Real h = 1e-6;
    auto fpPlus  = makeFlowPoint(gradU + h * dGradU);
    auto fpMinus = makeFlowPoint(gradU - h * dGradU);

    CarreauYasuda::Cache cPlus, cMinus;
    law.setCache(cPlus, fpPlus);
    law.setCache(cMinus, fpMinus);

    Math::SpatialMatrix<Real> tauPlus(2, 2), tauMinus(2, 2);
    law.getDeviatoricStress(tauPlus,  cPlus,  fpPlus);
    law.getDeviatoricStress(tauMinus, cMinus, fpMinus);

    const Math::SpatialMatrix<Real> dtau_fd = (tauPlus - tauMinus) / (2.0 * h);

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dtau_analytic(i, j), dtau_fd(i, j), 1e-8);
  }

  // =========================================================================
  // Bingham rheology
  // =========================================================================

  TEST(Rodin_Fluid_Bingham, ViscousLimitLargeM)
  {
    // For large m and tauY=0 the law reduces to Newtonian (mu_eff -> mu_p)
    const Real muP = 1.0;
    const Real tauY = 0.0;
    const Real m = 1e4;
    Bingham law(muP, tauY, m);

    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = 2.0;
    auto fp = makeFlowPoint(gradU);

    Bingham::Cache cache;
    law.setCache(cache, fp);

    EXPECT_NEAR(cache.effectiveViscosity, muP, 1e-8);
  }

  TEST(Rodin_Fluid_Bingham, StressIncreasesWithYieldStress)
  {
    // Higher yield stress => higher effective viscosity => larger deviatoric stress
    const Real gamma = 1.0;
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = gamma;
    auto fp = makeFlowPoint(gradU);

    Bingham law1(1.0, 1.0, 100.0);
    Bingham law2(1.0, 5.0, 100.0);

    Bingham::Cache c1, c2;
    law1.setCache(c1, fp);
    law2.setCache(c2, fp);

    Math::SpatialMatrix<Real> tau1(2, 2), tau2(2, 2);
    law1.getDeviatoricStress(tau1, c1, fp);
    law2.getDeviatoricStress(tau2, c2, fp);

    EXPECT_GT(std::abs(tau2(0, 1)), std::abs(tau1(0, 1)));
  }

  TEST(Rodin_Fluid_Bingham, TangentConsistencyFiniteDiff)
  {
    Bingham law(1.0, 2.0, 50.0, 1e-12);

    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = 1.0;
    auto fp = makeFlowPoint(gradU);

    Bingham::Cache cache;
    law.setCache(cache, fp);

    Math::SpatialMatrix<Real> dGradU(2, 2);
    dGradU.setZero();
    dGradU(0, 1) = 0.01;
    dGradU(1, 0) = 0.01;

    Math::SpatialMatrix<Real> dtau_analytic(2, 2);
    law.getTangent(dtau_analytic, cache, fp, dGradU);

    const Real h = 1e-6;
    auto fpPlus  = makeFlowPoint(gradU + h * dGradU);
    auto fpMinus = makeFlowPoint(gradU - h * dGradU);

    Bingham::Cache cPlus, cMinus;
    law.setCache(cPlus, fpPlus);
    law.setCache(cMinus, fpMinus);

    Math::SpatialMatrix<Real> tauPlus(2, 2), tauMinus(2, 2);
    law.getDeviatoricStress(tauPlus,  cPlus,  fpPlus);
    law.getDeviatoricStress(tauMinus, cMinus, fpMinus);

    const Math::SpatialMatrix<Real> dtau_fd = (tauPlus - tauMinus) / (2.0 * h);

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dtau_analytic(i, j), dtau_fd(i, j), 1e-8);
  }

  // =========================================================================
  // HerschelBulkley rheology
  // =========================================================================

  TEST(Rodin_Fluid_HerschelBulkley, ReducesToPowerLawWhenYieldStressZero)
  {
    // tauY=0 and n=0.5 => mu_eff = K * gdot^{n-1} (same as PowerLaw)
    const Real K = 1.0;
    const Real n = 0.5;
    const Real tauY = 0.0;
    const Real m = 1e4;
    HerschelBulkley hbLaw(K, n, tauY, m);
    PowerLaw pwLaw(K, n);

    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = 2.0;
    auto fp = makeFlowPoint(gradU);

    HerschelBulkley::Cache hbCache;
    hbLaw.setCache(hbCache, fp);

    PowerLaw::Cache pwCache;
    pwLaw.setCache(pwCache, fp);

    EXPECT_NEAR(hbCache.effectiveViscosity, pwCache.effectiveViscosity, 1e-6);
  }

  TEST(Rodin_Fluid_HerschelBulkley, TangentConsistencyFiniteDiff)
  {
    HerschelBulkley law(0.5, 0.6, 1.0, 50.0, 1e-12);

    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = 1.5;
    gradU(1, 0) = 0.3;
    auto fp = makeFlowPoint(gradU);

    HerschelBulkley::Cache cache;
    law.setCache(cache, fp);

    Math::SpatialMatrix<Real> dGradU(2, 2);
    dGradU.setZero();
    dGradU(0, 0) = 0.01;
    dGradU(1, 1) = -0.01;

    Math::SpatialMatrix<Real> dtau_analytic(2, 2);
    law.getTangent(dtau_analytic, cache, fp, dGradU);

    const Real h = 1e-6;
    auto fpPlus  = makeFlowPoint(gradU + h * dGradU);
    auto fpMinus = makeFlowPoint(gradU - h * dGradU);

    HerschelBulkley::Cache cPlus, cMinus;
    law.setCache(cPlus, fpPlus);
    law.setCache(cMinus, fpMinus);

    Math::SpatialMatrix<Real> tauPlus(2, 2), tauMinus(2, 2);
    law.getDeviatoricStress(tauPlus,  cPlus,  fpPlus);
    law.getDeviatoricStress(tauMinus, cMinus, fpMinus);

    const Math::SpatialMatrix<Real> dtau_fd = (tauPlus - tauMinus) / (2.0 * h);

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        EXPECT_NEAR(dtau_analytic(i, j), dtau_fd(i, j), 1e-8);
  }

  // =========================================================================
  // CauchyStress
  // =========================================================================

  TEST(Rodin_Fluid_CauchyStress, ZeroPressure)
  {
    // sigma = tau  when p=0
    const Real mu = 1.5;
    const Real gamma = 2.0;
    Newtonian law(mu);

    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    gradU(0, 1) = gamma;
    auto fp = makeFlowPoint(gradU);

    Newtonian::Cache cache;
    law.setCache(cache, fp);

    CauchyStress<Newtonian> cs(law);
    Math::SpatialMatrix<Real> sigma(2, 2);
    cs.getCauchyStress(sigma, cache, fp);

    // Deviatoric: tau(0,1) = tau(1,0) = mu * gamma, diagonals = 0
    EXPECT_NEAR(sigma(0, 1), mu * gamma, 1e-14);
    EXPECT_NEAR(sigma(1, 0), mu * gamma, 1e-14);
    EXPECT_NEAR(sigma(0, 0), 0.0,        1e-14);
    EXPECT_NEAR(sigma(1, 1), 0.0,        1e-14);
  }

  TEST(Rodin_Fluid_CauchyStress, WithPressure)
  {
    // sigma = tau - p*I
    const Real mu = 1.0;
    const Real gamma = 0.0;
    const Real p = 5.0;
    Newtonian law(mu);

    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(2, 2));
    fp.setPressure(p);

    Newtonian::Cache cache;
    law.setCache(cache, fp);

    CauchyStress<Newtonian> cs(law);
    Math::SpatialMatrix<Real> sigma(2, 2);
    cs.getCauchyStress(sigma, cache, fp);

    EXPECT_NEAR(sigma(0, 0), -p,  1e-14);
    EXPECT_NEAR(sigma(1, 1), -p,  1e-14);
    EXPECT_NEAR(sigma(0, 1),  0.0, 1e-14);
    EXPECT_NEAR(sigma(1, 0),  0.0, 1e-14);
  }

  // =========================================================================
  // Tags / Input helpers
  // =========================================================================

  TEST(Rodin_Fluid_Tags, AllStandardTagsSetAndGet)
  {
    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(3, 3));
    fp.set<Tags::Density>(1060.0);
    fp.set<Tags::DynamicViscosity>(0.0027);
    fp.set<Tags::WallDistance>(0.001);
    fp.set<Tags::YieldStress>(0.5);
    fp.set<Tags::ConsistencyIndex>(0.35);
    fp.set<Tags::FlowIndex>(0.61);
    fp.set<Tags::Regularization>(1e-4);
    fp.set<Tags::TurbulentKineticEnergy>(0.01);
    fp.set<Tags::DissipationRate>(0.001);
    fp.set<Tags::SpecificDissipationRate>(1.5);

    EXPECT_NEAR(fp.get<Tags::Density>(),                 1060.0, 1e-14);
    EXPECT_NEAR(fp.get<Tags::DynamicViscosity>(),        0.0027, 1e-15);
    EXPECT_NEAR(fp.get<Tags::WallDistance>(),             0.001, 1e-15);
    EXPECT_NEAR(fp.get<Tags::YieldStress>(),               0.5,  1e-15);
    EXPECT_NEAR(fp.get<Tags::ConsistencyIndex>(),          0.35, 1e-15);
    EXPECT_NEAR(fp.get<Tags::FlowIndex>(),                 0.61, 1e-15);
    EXPECT_NEAR(fp.get<Tags::Regularization>(),           1e-4,  1e-18);
    EXPECT_NEAR(fp.get<Tags::TurbulentKineticEnergy>(),   0.01,  1e-15);
    EXPECT_NEAR(fp.get<Tags::DissipationRate>(),          0.001, 1e-15);
    EXPECT_NEAR(fp.get<Tags::SpecificDissipationRate>(),   1.5,  1e-14);
  }

  TEST(Rodin_Fluid_Input, InputFunctionPopulates)
  {
    auto fp = makeFlowPoint(Math::SpatialMatrix<Real>::Zero(2, 2));

    InputFunction inject = [](FlowPoint& p) {
      p.set<Tags::Density>(1000.0);
      p.setPressure(101325.0);
    };
    inject(fp);

    EXPECT_TRUE(fp.has<Tags::Density>());
    EXPECT_NEAR(fp.get<Tags::Density>(), 1000.0, 1e-14);
    EXPECT_TRUE(fp.hasPressure());
    EXPECT_NEAR(fp.getPressure().value(), 101325.0, 1e-9);
  }
}
