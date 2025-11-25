/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational/H1/JacobiPolynomial.h"

using namespace Rodin;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  //==========================================================================
  // Jacobi Polynomial with α = β = 0 (Legendre polynomials)
  //==========================================================================

  TEST(JacobiPolynomial, P00_K0_IsLegendre)
  {
    // P_0^{0,0}(x) = 1 (same as Legendre)
    Real P, dP;

    JacobiPolynomial<0>::getValue(P, dP, 0.0, 0.0, 0.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    JacobiPolynomial<0>::getValue(P, dP, 0.0, 0.0, 0.5);
    EXPECT_NEAR(P, 1.0, 1e-14);

    JacobiPolynomial<0>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);
  }

  TEST(JacobiPolynomial, P00_K1_IsLegendre)
  {
    // P_1^{0,0}(x) = x (same as Legendre)
    Real P, dP;

    JacobiPolynomial<1>::getValue(P, dP, 0.0, 0.0, 0.0);
    EXPECT_NEAR(P, 0.0, 1e-14);

    JacobiPolynomial<1>::getValue(P, dP, 0.0, 0.0, 0.5);
    EXPECT_NEAR(P, 0.5, 1e-14);

    JacobiPolynomial<1>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    JacobiPolynomial<1>::getValue(P, dP, 0.0, 0.0, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-14);
  }

  TEST(JacobiPolynomial, P00_K2_IsLegendre)
  {
    // P_2^{0,0}(x) = (3x^2 - 1)/2 (same as Legendre)
    Real P, dP;

    JacobiPolynomial<2>::getValue(P, dP, 0.0, 0.0, 0.0);
    EXPECT_NEAR(P, -0.5, 1e-14);

    Real x = 0.5;
    Real expected = (3.0 * x * x - 1.0) / 2.0;
    JacobiPolynomial<2>::getValue(P, dP, 0.0, 0.0, x);
    EXPECT_NEAR(P, expected, 1e-14);
  }

  //==========================================================================
  // Jacobi Polynomial Derivative Tests
  //==========================================================================

  TEST(JacobiPolynomial, Derivative_K0)
  {
    // d/dx P_0^{α,β}(x) = 0 for all α, β
    Real P, dP;

    JacobiPolynomial<0>::getValue(P, dP, 0.0, 0.0, 0.5);
    EXPECT_NEAR(dP, 0.0, 1e-14);

    JacobiPolynomial<0>::getValue(P, dP, 1.0, 0.0, 0.5);
    EXPECT_NEAR(dP, 0.0, 1e-14);

    JacobiPolynomial<0>::getValue(P, dP, 2.0, 1.0, 0.5);
    EXPECT_NEAR(dP, 0.0, 1e-14);
  }

  TEST(JacobiPolynomial, Derivative_K1_Alpha0_Beta0)
  {
    // d/dx P_1^{0,0}(x) = 1
    Real P, dP;

    JacobiPolynomial<1>::getValue(P, dP, 0.0, 0.0, 0.0);
    EXPECT_NEAR(dP, 1.0, 1e-14);

    JacobiPolynomial<1>::getValue(P, dP, 0.0, 0.0, 0.5);
    EXPECT_NEAR(dP, 1.0, 1e-14);
  }

  TEST(JacobiPolynomial, Derivative_K2_Alpha0_Beta0)
  {
    // d/dx P_2^{0,0}(x) = 3x (Legendre)
    Real P, dP;

    JacobiPolynomial<2>::getValue(P, dP, 0.0, 0.0, 0.0);
    EXPECT_NEAR(dP, 0.0, 1e-14);

    JacobiPolynomial<2>::getValue(P, dP, 0.0, 0.0, 0.5);
    EXPECT_NEAR(dP, 1.5, 1e-14);

    JacobiPolynomial<2>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(dP, 3.0, 1e-14);
  }

  //==========================================================================
  // Jacobi Polynomial with α = 1, β = 0 (used in Dubiner basis)
  //==========================================================================

  TEST(JacobiPolynomial, P10_K0)
  {
    // P_0^{1,0}(x) = 1
    Real P, dP;

    JacobiPolynomial<0>::getValue(P, dP, 1.0, 0.0, 0.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    JacobiPolynomial<0>::getValue(P, dP, 1.0, 0.0, 0.5);
    EXPECT_NEAR(P, 1.0, 1e-14);
  }

  TEST(JacobiPolynomial, P10_K1)
  {
    // P_1^{1,0}(x) = (α - β + (α + β + 2)x) / 2 = (1 + 3x) / 2
    Real P, dP;

    JacobiPolynomial<1>::getValue(P, dP, 1.0, 0.0, 0.0);
    EXPECT_NEAR(P, 0.5, 1e-14);  // (1 + 0)/2 = 0.5

    JacobiPolynomial<1>::getValue(P, dP, 1.0, 0.0, 1.0);
    EXPECT_NEAR(P, 2.0, 1e-14);  // (1 + 3)/2 = 2

    JacobiPolynomial<1>::getValue(P, dP, 1.0, 0.0, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-14);  // (1 - 3)/2 = -1
  }

  //==========================================================================
  // Jacobi Polynomial with higher α, β values (used in Dubiner for triangles)
  //==========================================================================

  TEST(JacobiPolynomial, P30_K0)
  {
    // P_0^{3,0}(x) = 1
    Real P, dP;

    JacobiPolynomial<0>::getValue(P, dP, 3.0, 0.0, 0.0);
    EXPECT_NEAR(P, 1.0, 1e-14);
  }

  TEST(JacobiPolynomial, P30_K1)
  {
    // P_1^{3,0}(x) = (α - β + (α + β + 2)x) / 2 = (3 + 5x) / 2
    Real P, dP;

    JacobiPolynomial<1>::getValue(P, dP, 3.0, 0.0, 0.0);
    EXPECT_NEAR(P, 1.5, 1e-14);  // (3 + 0)/2 = 1.5

    JacobiPolynomial<1>::getValue(P, dP, 3.0, 0.0, 1.0);
    EXPECT_NEAR(P, 4.0, 1e-14);  // (3 + 5)/2 = 4
  }

  //==========================================================================
  // Jacobi Polynomial Endpoint Tests
  //==========================================================================

  TEST(JacobiPolynomial, EndpointValue_x1)
  {
    // P_K^{α,β}(1) = (K + α choose K) = Γ(K + α + 1) / (Γ(K + 1) Γ(α + 1))
    // For α = 0: P_K^{0,β}(1) = 1
    Real P, dP;

    JacobiPolynomial<0>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    JacobiPolynomial<1>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    JacobiPolynomial<2>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    JacobiPolynomial<5>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-13);
  }

  TEST(JacobiPolynomial, EndpointValue_xm1)
  {
    // P_K^{α,β}(-1) = (-1)^K (K + β choose K)
    // For β = 0: P_K^{α,0}(-1) = (-1)^K
    Real P, dP;

    JacobiPolynomial<0>::getValue(P, dP, 1.0, 0.0, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);  // (-1)^0 = 1

    JacobiPolynomial<1>::getValue(P, dP, 1.0, 0.0, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-14);  // (-1)^1 = -1

    JacobiPolynomial<2>::getValue(P, dP, 1.0, 0.0, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);  // (-1)^2 = 1

    JacobiPolynomial<3>::getValue(P, dP, 1.0, 0.0, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-14);  // (-1)^3 = -1
  }

  //==========================================================================
  // Jacobi Polynomial Derivative Identity
  //==========================================================================

  TEST(JacobiPolynomial, DerivativeIdentity)
  {
    // d/dx P_K^{α,β}(x) = (K + α + β + 1)/2 * P_{K-1}^{α+1,β+1}(x)
    Real P, dP;

    // Test for K=2, α=0, β=0
    // d/dx P_2^{0,0}(x) = (2 + 0 + 0 + 1)/2 * P_1^{1,1}(x) = 1.5 * P_1^{1,1}(x)
    Real x = 0.3;
    JacobiPolynomial<2>::getValue(P, dP, 0.0, 0.0, x);

    Real P1_11, dP1_11;
    JacobiPolynomial<1>::getValue(P1_11, dP1_11, 1.0, 1.0, x);

    Real expected_deriv = 1.5 * P1_11;
    EXPECT_NEAR(dP, expected_deriv, 1e-13);
  }

  //==========================================================================
  // Jacobi Polynomial Recurrence Consistency
  //==========================================================================

  TEST(JacobiPolynomial, RecurrenceConsistency_K5)
  {
    // Verify consistency at multiple points for K=5
    Real P, dP;
    const double tol = 1e-12;

    // Test at several points
    std::vector<Real> test_points = {-0.9, -0.5, 0.0, 0.5, 0.9};

    for (Real x : test_points)
    {
      JacobiPolynomial<5>::getValue(P, dP, 0.0, 0.0, x);

      // Verify P_5(x) at x=0 is 0 (odd polynomial)
      // and P_5(±1) = (±1)^5 = ±1
      if (std::abs(x) < 1e-10)
      {
        EXPECT_NEAR(P, 0.0, tol);
      }
    }

    JacobiPolynomial<5>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, tol);

    JacobiPolynomial<5>::getValue(P, dP, 0.0, 0.0, -1.0);
    EXPECT_NEAR(P, -1.0, tol);
  }

  //==========================================================================
  // Jacobi Polynomial Higher Orders (used in Dubiner tetrahedron)
  //==========================================================================

  TEST(JacobiPolynomial, HigherAlpha_K2)
  {
    // Test P_2^{5,0}(x) - used in Dubiner tetrahedron for (p,q,r) modes
    Real P, dP;

    JacobiPolynomial<2>::getValue(P, dP, 5.0, 0.0, 0.0);
    // This should be computable via recurrence
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));

    JacobiPolynomial<2>::getValue(P, dP, 5.0, 0.0, 1.0);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));
  }

  TEST(JacobiPolynomial, VeryHighAlpha_K3)
  {
    // Test P_3^{10,0}(x) - verifying stability
    Real P, dP;

    JacobiPolynomial<3>::getValue(P, dP, 10.0, 0.0, 0.0);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));
  }

  //==========================================================================
  // Higher Order Tests (K = 5, 6)
  //==========================================================================

  TEST(JacobiPolynomial, P00_K5_EndpointValues)
  {
    // P_5^{0,0}(1) = 1, P_5^{0,0}(-1) = -1 (Legendre)
    Real P, dP;

    JacobiPolynomial<5>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-13);

    JacobiPolynomial<5>::getValue(P, dP, 0.0, 0.0, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-13);
  }

  TEST(JacobiPolynomial, P00_K6_EndpointValues)
  {
    // P_6^{0,0}(1) = 1, P_6^{0,0}(-1) = 1 (Legendre)
    Real P, dP;

    JacobiPolynomial<6>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-13);

    JacobiPolynomial<6>::getValue(P, dP, 0.0, 0.0, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-13);  // Even polynomial
  }

  TEST(JacobiPolynomial, HigherAlpha_K5)
  {
    // Test P_5^{7,0}(x) - used in high-order Dubiner tetrahedron
    Real P, dP;

    JacobiPolynomial<5>::getValue(P, dP, 7.0, 0.0, 0.0);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));

    JacobiPolynomial<5>::getValue(P, dP, 7.0, 0.0, 0.5);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));

    JacobiPolynomial<5>::getValue(P, dP, 7.0, 0.0, 1.0);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));
  }

  TEST(JacobiPolynomial, HigherAlpha_K6)
  {
    // Test P_6^{9,0}(x) - used in high-order Dubiner tetrahedron
    Real P, dP;

    JacobiPolynomial<6>::getValue(P, dP, 9.0, 0.0, 0.0);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));

    JacobiPolynomial<6>::getValue(P, dP, 9.0, 0.0, 0.5);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));
  }

  TEST(JacobiPolynomial, DerivativeIdentity_K5)
  {
    // d/dx P_5^{α,β}(x) = (5 + α + β + 1)/2 * P_4^{α+1,β+1}(x)
    Real P, dP;
    Real x = 0.3;

    JacobiPolynomial<5>::getValue(P, dP, 0.0, 0.0, x);

    Real P4_11, dP4_11;
    JacobiPolynomial<4>::getValue(P4_11, dP4_11, 1.0, 1.0, x);

    Real expected_deriv = 3.0 * P4_11;  // (5 + 0 + 0 + 1)/2 = 3
    EXPECT_NEAR(dP, expected_deriv, 1e-12);
  }

  TEST(JacobiPolynomial, DerivativeIdentity_K6)
  {
    // d/dx P_6^{α,β}(x) = (6 + α + β + 1)/2 * P_5^{α+1,β+1}(x)
    Real P, dP;
    Real x = 0.4;

    JacobiPolynomial<6>::getValue(P, dP, 0.0, 0.0, x);

    Real P5_11, dP5_11;
    JacobiPolynomial<5>::getValue(P5_11, dP5_11, 1.0, 1.0, x);

    Real expected_deriv = 3.5 * P5_11;  // (6 + 0 + 0 + 1)/2 = 3.5
    EXPECT_NEAR(dP, expected_deriv, 1e-11);
  }

  TEST(JacobiPolynomial, RecurrenceConsistency_K6)
  {
    // Verify consistency at multiple points for K=6
    Real P, dP;
    const double tol = 1e-12;

    // Test at several points
    std::vector<Real> test_points = {-0.9, -0.5, 0.0, 0.5, 0.9};

    for (Real x : test_points)
    {
      JacobiPolynomial<6>::getValue(P, dP, 0.0, 0.0, x);
      // Just verify no NaN/Inf
      EXPECT_FALSE(std::isnan(P));
      EXPECT_FALSE(std::isinf(P));
      EXPECT_FALSE(std::isnan(dP));
      EXPECT_FALSE(std::isinf(dP));
    }

    // P_6(0) for Legendre (α=β=0)
    JacobiPolynomial<6>::getValue(P, dP, 0.0, 0.0, 0.0);
    Real expected = -5.0 / 16.0;  // Known value for P_6(0)
    EXPECT_NEAR(P, expected, tol);
  }
}
