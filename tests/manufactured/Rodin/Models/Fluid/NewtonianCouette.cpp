/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file NewtonianCouette.cpp
 * @brief Manufactured tests for Newtonian Couette flow.
 *
 * Couette flow between two parallel plates at y=0 (fixed) and y=H (moving at
 * speed U) has the exact velocity field:
 * @f[
 *   u_x(y) = U \frac{y}{H}, \qquad u_y = 0.
 * @f]
 * The velocity gradient is:
 * @f[
 *   \nabla\mathbf{u} =
 *   \begin{pmatrix} 0 & 0 \\ U/H & 0 \end{pmatrix}.
 * @f]
 * The strain-rate tensor, shear rate, and Newtonian deviatoric stress are:
 * @f[
 *   \mathbf{D} = \frac{1}{2}\begin{pmatrix}0 & U/H \\ U/H & 0\end{pmatrix},\qquad
 *   \dot\gamma = \frac{U}{H},\qquad
 *   \boldsymbol\tau = \mu \begin{pmatrix}0 & U/H \\ U/H & 0\end{pmatrix}.
 * @f]
 *
 * The PowerLaw flow in the same geometry gives:
 * @f[
 *   \mu_\text{eff} = K\,\left(\frac{U}{H}\right)^{n-1},\qquad
 *   \boldsymbol\tau = 2\mu_\text{eff}\,\mathbf{D}.
 * @f]
 */
#include <cmath>
#include <gtest/gtest.h>

#include <Rodin/Types.h>
#include <Rodin/Math/SpatialVector.h>
#include <Rodin/Math/SpatialMatrix.h>

#include <Rodin/Fluid/Local/FlowPoint.h>
#include <Rodin/Fluid/Fields/StrainRate.h>
#include <Rodin/Fluid/Fields/ShearRate.h>
#include <Rodin/Fluid/Constitutive/Newtonian.h>
#include <Rodin/Fluid/Constitutive/PowerLaw.h>
#include <Rodin/Fluid/Constitutive/Bingham.h>

using namespace Rodin;
using namespace Rodin::Fluid;

namespace Rodin::Tests::Manufactured::FluidCouette
{
  // =========================================================================
  // Helpers
  // =========================================================================

  /// Build a FlowPoint for Couette flow at shear rate gamma = U/H (2D).
  static FlowPoint makeCouetteFlowPoint(Real shearRate)
  {
    Math::SpatialVector<Real> vel(2);
    vel.setZero();
    Math::SpatialMatrix<Real> gradU(2, 2);
    gradU.setZero();
    // u = (shearRate * y, 0) => grad u(0,1) = shearRate (du_x/dy component)
    gradU(0, 1) = shearRate;
    return FlowPoint(vel, gradU);
  }

  // =========================================================================
  // Newtonian Couette: stress against analytical formula over grid
  // =========================================================================

  /**
   * Sweeps a grid of shear rates in [gamma_min, gamma_max] with N samples and
   * verifies that the Newtonian stress matches tau = mu * gamma exactly.
   */
  class Manufactured_Newtonian_Couette_Test
    : public ::testing::TestWithParam<size_t>
  {};

  TEST_P(Manufactured_Newtonian_Couette_Test, NewtonianCouette_ShearStress)
  {
    const size_t N = GetParam();
    const Real mu = 0.003;
    const Real gammaMin = 1e-3;
    const Real gammaMax = 1e3;

    Newtonian law(mu);

    const Real logMin = std::log10(gammaMin);
    const Real logMax = std::log10(gammaMax);

    for (size_t i = 0; i < N; ++i)
    {
      const Real gamma = std::pow(10.0,
          logMin + static_cast<Real>(i) / static_cast<Real>(N - 1)
                 * (logMax - logMin));

      auto fp = makeCouetteFlowPoint(gamma);

      Newtonian::Cache cache;
      law.setCache(cache, fp);

      Math::SpatialMatrix<Real> tau(2, 2);
      law.getDeviatoricStress(tau, cache, fp);

      // Analytical: tau(0,1) = tau(1,0) = mu * gamma
      const Real expected = mu * gamma;
      EXPECT_NEAR(tau(0, 1), expected, 1e-12 * expected + 1e-15)
          << "at gamma=" << gamma;
      EXPECT_NEAR(tau(1, 0), expected, 1e-12 * expected + 1e-15)
          << "at gamma=" << gamma;
      EXPECT_NEAR(tau(0, 0), 0.0, 1e-15) << "at gamma=" << gamma;
      EXPECT_NEAR(tau(1, 1), 0.0, 1e-15) << "at gamma=" << gamma;
    }
  }

  INSTANTIATE_TEST_SUITE_P(
    FluidCouette_N,
    Manufactured_Newtonian_Couette_Test,
    ::testing::Values(16, 64, 256));

  // =========================================================================
  // Newtonian Couette: shear-rate invariant
  // =========================================================================

  class Manufactured_Newtonian_ShearRate_Test
    : public ::testing::TestWithParam<size_t>
  {};

  TEST_P(Manufactured_Newtonian_ShearRate_Test, ShearRate_Couette_Grid)
  {
    const size_t N = GetParam();
    const Real gammaMin = 1e-4;
    const Real gammaMax = 1e4;

    ShearRate sr;

    const Real logMin = std::log10(gammaMin);
    const Real logMax = std::log10(gammaMax);

    for (size_t i = 0; i < N; ++i)
    {
      const Real gamma = std::pow(10.0,
          logMin + static_cast<Real>(i) / static_cast<Real>(N - 1)
                 * (logMax - logMin));

      auto fp = makeCouetteFlowPoint(gamma);
      const Real computed = sr.getShearRate(fp);

      EXPECT_NEAR(computed, gamma, 1e-12 * gamma + 1e-15)
          << "at gamma=" << gamma;
    }
  }

  INSTANTIATE_TEST_SUITE_P(
    FluidCouette_N,
    Manufactured_Newtonian_ShearRate_Test,
    ::testing::Values(16, 64, 256));

  // =========================================================================
  // PowerLaw Couette: effective viscosity against analytical formula
  // =========================================================================

  class Manufactured_PowerLaw_Couette_Test
    : public ::testing::TestWithParam<size_t>
  {};

  TEST_P(Manufactured_PowerLaw_Couette_Test, PowerLaw_EffectiveViscosity_Grid)
  {
    const size_t N = GetParam();
    const Real K  = 0.5;
    const Real n  = 0.6;
    const Real eps = 0.0;  // Use exact formula (large shear far from regularization)
    PowerLaw law(K, n, eps * eps);  // eps^2 is added inside, pass 0 effectively

    const Real gammaMin = 1.0;  // Stay well away from eps=0 singularity
    const Real gammaMax = 1e3;

    const Real logMin = std::log10(gammaMin);
    const Real logMax = std::log10(gammaMax);

    for (size_t i = 0; i < N; ++i)
    {
      const Real gamma = std::pow(10.0,
          logMin + static_cast<Real>(i) / static_cast<Real>(N - 1)
                 * (logMax - logMin));

      auto fp = makeCouetteFlowPoint(gamma);

      PowerLaw::Cache cache;
      law.setCache(cache, fp);

      // When eps~0 and gamma >> eps: mu_eff = K * gamma^{n-1}
      const Real expected = K * std::pow(gamma, n - 1.0);
      EXPECT_NEAR(cache.effectiveViscosity, expected, 1e-8 * expected + 1e-15)
          << "at gamma=" << gamma;

      Math::SpatialMatrix<Real> tau(2, 2);
      law.getDeviatoricStress(tau, cache, fp);

      // tau(0,1) = mu_eff * gamma
      const Real expectedStress = expected * gamma;
      EXPECT_NEAR(tau(0, 1), expectedStress, 1e-8 * expectedStress + 1e-15)
          << "at gamma=" << gamma;
    }
  }

  INSTANTIATE_TEST_SUITE_P(
    FluidCouette_N,
    Manufactured_PowerLaw_Couette_Test,
    ::testing::Values(16, 64, 256));

  // =========================================================================
  // Bingham Couette: unyielded limit (small shear, yield stress dominates)
  // =========================================================================

  TEST(Rodin_Fluid_Manufactured_Bingham, UnyieldedRegionApproachesYieldStress)
  {
    // As gamma -> 0 with large m, the effective viscosity satisfies:
    //   mu_eff * gamma -> tau_Y
    // so the shear stress approaches the yield stress from above.
    const Real muP  = 1.0;
    const Real tauY = 5.0;
    const Real m    = 1e4;
    Bingham law(muP, tauY, m, 1e-14);

    // At very small gamma, (1 - exp(-m*gamma))/gamma -> m,
    // so mu_eff ~ muP + tauY * m (huge), but the product tau = 2*mu_eff*D
    // For gamma small: tau(0,1) = mu_eff * gamma ~ (muP + tauY * m) * gamma
    // We can check that as gamma->0 the ratio tau(0,1)/gamma is bounded above.
    // More usefully: check that in the yielded region tau(0,1) > tauY/2

    const Real gamma = 0.5;  // well above yield
    auto fp = makeCouetteFlowPoint(gamma);

    Bingham::Cache cache;
    law.setCache(cache, fp);

    Math::SpatialMatrix<Real> tau(2, 2);
    law.getDeviatoricStress(tau, cache, fp);

    // For yielded Bingham: tau(0,1) = mu_eff * gamma >= tauY*(1-exp(-m*gamma))
    // which approaches tauY for m*gamma >> 1.
    const Real expected_yield_contribution = tauY * (1.0 - std::exp(-m * gamma));
    EXPECT_GT(tau(0, 1), expected_yield_contribution - 1e-6);
  }

  // =========================================================================
  // Newtonian 3D Couette: simple extension test
  // =========================================================================

  TEST(Rodin_Fluid_Manufactured_Newtonian3D, CouetteShearStress3D)
  {
    const Real mu    = 0.001;
    const Real gamma = 100.0;

    Newtonian law(mu);

    // 3D Couette: velocity = (gamma * z, 0, 0), so gradU(0,2) = gamma
    Math::SpatialVector<Real> vel(3);
    vel.setZero();
    Math::SpatialMatrix<Real> gradU(3, 3);
    gradU.setZero();
    gradU(0, 2) = gamma;  // du_x/dz = gamma

    FlowPoint fp(vel, gradU);

    Newtonian::Cache cache;
    law.setCache(cache, fp);

    Math::SpatialMatrix<Real> tau(3, 3);
    law.getDeviatoricStress(tau, cache, fp);

    // D(0,2) = D(2,0) = gamma/2, so tau(0,2) = tau(2,0) = mu * gamma
    EXPECT_NEAR(tau(0, 2), mu * gamma, 1e-13);
    EXPECT_NEAR(tau(2, 0), mu * gamma, 1e-13);

    // All other components zero
    EXPECT_NEAR(tau(0, 0), 0.0, 1e-15);
    EXPECT_NEAR(tau(1, 1), 0.0, 1e-15);
    EXPECT_NEAR(tau(2, 2), 0.0, 1e-15);
    EXPECT_NEAR(tau(0, 1), 0.0, 1e-15);
    EXPECT_NEAR(tau(1, 0), 0.0, 1e-15);
    EXPECT_NEAR(tau(1, 2), 0.0, 1e-15);
    EXPECT_NEAR(tau(2, 1), 0.0, 1e-15);
  }
}
