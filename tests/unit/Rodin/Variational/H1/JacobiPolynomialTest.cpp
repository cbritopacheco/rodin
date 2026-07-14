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
  // Jacobi Polynomial with \alpha = \beta = 0 (Legendre polynomials)
  //==========================================================================

  /// @brief P_0^{0,0}(x) = 1 (same as Legendre).
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

  /// @brief P_1^{0,0}(x) = x (same as Legendre).
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

  /// @brief P_2^{0,0}(x) = (3x^2 - 1)/2 (same as Legendre).
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

  /// @brief Verifies @f$ \frac{d}{dx} P_0^{\alpha,\beta}(x) = 0 @f$ for all tested @f$ \alpha @f$ and @f$ \beta @f$.
  TEST(JacobiPolynomial, Derivative_K0)
  {
    // d/dx P_0^{\alpha,\beta}(x) = 0 for all \alpha, \beta
    Real P, dP;

    JacobiPolynomial<0>::getValue(P, dP, 0.0, 0.0, 0.5);
    EXPECT_NEAR(dP, 0.0, 1e-14);

    JacobiPolynomial<0>::getValue(P, dP, 1.0, 0.0, 0.5);
    EXPECT_NEAR(dP, 0.0, 1e-14);

    JacobiPolynomial<0>::getValue(P, dP, 2.0, 1.0, 0.5);
    EXPECT_NEAR(dP, 0.0, 1e-14);
  }

  /// @brief D/dx P_1^{0,0}(x) = 1.
  TEST(JacobiPolynomial, Derivative_K1_Alpha0_Beta0)
  {
    // d/dx P_1^{0,0}(x) = 1
    Real P, dP;

    JacobiPolynomial<1>::getValue(P, dP, 0.0, 0.0, 0.0);
    EXPECT_NEAR(dP, 1.0, 1e-14);

    JacobiPolynomial<1>::getValue(P, dP, 0.0, 0.0, 0.5);
    EXPECT_NEAR(dP, 1.0, 1e-14);
  }

  /// @brief D/dx P_2^{0,0}(x) = 3x (Legendre).
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
  // Jacobi Polynomial with \alpha = 1, \beta = 0 (used in Dubiner basis)
  //==========================================================================

  /// @brief P_0^{1,0}(x) = 1.
  TEST(JacobiPolynomial, P10_K0)
  {
    // P_0^{1,0}(x) = 1
    Real P, dP;

    JacobiPolynomial<0>::getValue(P, dP, 1.0, 0.0, 0.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    JacobiPolynomial<0>::getValue(P, dP, 1.0, 0.0, 0.5);
    EXPECT_NEAR(P, 1.0, 1e-14);
  }

  /// @brief Verifies @f$ P_1^{1,0}(x) = (1 + 3x) / 2 @f$.
  TEST(JacobiPolynomial, P10_K1)
  {
    // P_1^{1,0}(x) = (\alpha - \beta + (\alpha + \beta + 2)x) / 2 = (1 + 3x) / 2
    Real P, dP;

    JacobiPolynomial<1>::getValue(P, dP, 1.0, 0.0, 0.0);
    EXPECT_NEAR(P, 0.5, 1e-14);  // (1 + 0)/2 = 0.5

    JacobiPolynomial<1>::getValue(P, dP, 1.0, 0.0, 1.0);
    EXPECT_NEAR(P, 2.0, 1e-14);  // (1 + 3)/2 = 2

    JacobiPolynomial<1>::getValue(P, dP, 1.0, 0.0, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-14);  // (1 - 3)/2 = -1
  }

  //==========================================================================
  // Jacobi Polynomial with higher \alpha, \beta values (used in Dubiner for triangles)
  //==========================================================================

  /// @brief P_0^{3,0}(x) = 1.
  TEST(JacobiPolynomial, P30_K0)
  {
    // P_0^{3,0}(x) = 1
    Real P, dP;

    JacobiPolynomial<0>::getValue(P, dP, 3.0, 0.0, 0.0);
    EXPECT_NEAR(P, 1.0, 1e-14);
  }

  /// @brief Verifies @f$ P_1^{3,0}(x) = (3 + 5x) / 2 @f$.
  TEST(JacobiPolynomial, P30_K1)
  {
    // P_1^{3,0}(x) = (\alpha - \beta + (\alpha + \beta + 2)x) / 2 = (3 + 5x) / 2
    Real P, dP;

    JacobiPolynomial<1>::getValue(P, dP, 3.0, 0.0, 0.0);
    EXPECT_NEAR(P, 1.5, 1e-14);  // (3 + 0)/2 = 1.5

    JacobiPolynomial<1>::getValue(P, dP, 3.0, 0.0, 1.0);
    EXPECT_NEAR(P, 4.0, 1e-14);  // (3 + 5)/2 = 4
  }

  //==========================================================================
  // Jacobi Polynomial Endpoint Tests
  //==========================================================================

  /// @brief Verifies the endpoint identity @f$ P_K^{\alpha,\beta}(1) = \Gamma(K+\alpha+1)/(\Gamma(K+1)\Gamma(\alpha+1)) @f$.
  TEST(JacobiPolynomial, EndpointValue_x1)
  {
    // P_K^{\alpha,\beta}(1) = (K + \alpha choose K) = \Gamma(K + \alpha + 1) / (\Gamma(K + 1) \Gamma(\alpha + 1))
    // For \alpha = 0: P_K^{0,\beta}(1) = 1
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

  /// @brief Verifies the endpoint identity @f$ P_K^{\alpha,\beta}(-1) = (-1)^K {K+\beta \choose K} @f$.
  TEST(JacobiPolynomial, EndpointValue_xm1)
  {
    // P_K^{\alpha,\beta}(-1) = (-1)^K (K + \beta choose K)
    // For \beta = 0: P_K^{\alpha,0}(-1) = (-1)^K
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

  /// @brief Verifies @f$ \frac{d}{dx}P_K^{\alpha,\beta}(x) = \frac{K+\alpha+\beta+1}{2} P_{K-1}^{\alpha+1,\beta+1}(x) @f$.
  TEST(JacobiPolynomial, DerivativeIdentity)
  {
    // d/dx P_K^{\alpha,\beta}(x) = (K + \alpha + \beta + 1)/2 * P_{K-1}^{\alpha+1,\beta+1}(x)
    Real P, dP;

    // Test for K=2, \alpha=0, \beta=0
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

  /// @brief Verify consistency at multiple points for K=5.
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

  /// @brief Test P_2^{5,0}(x) - used in Dubiner tetrahedron for (p,q,r) modes.
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

  /// @brief Test P_3^{10,0}(x) - verifying stability.
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

  /// @brief P_5^{0,0}(1) = 1, P_5^{0,0}(-1) = -1 (Legendre).
  TEST(JacobiPolynomial, P00_K5_EndpointValues)
  {
    // P_5^{0,0}(1) = 1, P_5^{0,0}(-1) = -1 (Legendre)
    Real P, dP;

    JacobiPolynomial<5>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-13);

    JacobiPolynomial<5>::getValue(P, dP, 0.0, 0.0, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-13);
  }

  /// @brief P_6^{0,0}(1) = 1, P_6^{0,0}(-1) = 1 (Legendre).
  TEST(JacobiPolynomial, P00_K6_EndpointValues)
  {
    // P_6^{0,0}(1) = 1, P_6^{0,0}(-1) = 1 (Legendre)
    Real P, dP;

    JacobiPolynomial<6>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-13);

    JacobiPolynomial<6>::getValue(P, dP, 0.0, 0.0, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-13);  // Even polynomial
  }

  /// @brief Test P_5^{7,0}(x) - used in high-order Dubiner tetrahedron.
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

  /// @brief Test P_6^{9,0}(x) - used in high-order Dubiner tetrahedron.
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

  /// @brief Verifies @f$ \frac{d}{dx}P_5^{\alpha,\beta}(x) = \frac{6+\alpha+\beta}{2} P_4^{\alpha+1,\beta+1}(x) @f$.
  TEST(JacobiPolynomial, DerivativeIdentity_K5)
  {
    // d/dx P_5^{\alpha,\beta}(x) = (5 + \alpha + \beta + 1)/2 * P_4^{\alpha+1,\beta+1}(x)
    Real P, dP;
    Real x = 0.3;

    JacobiPolynomial<5>::getValue(P, dP, 0.0, 0.0, x);

    Real P4_11, dP4_11;
    JacobiPolynomial<4>::getValue(P4_11, dP4_11, 1.0, 1.0, x);

    Real expected_deriv = 3.0 * P4_11;  // (5 + 0 + 0 + 1)/2 = 3
    EXPECT_NEAR(dP, expected_deriv, 1e-12);
  }

  /// @brief Verifies @f$ \frac{d}{dx}P_6^{\alpha,\beta}(x) = \frac{7+\alpha+\beta}{2} P_5^{\alpha+1,\beta+1}(x) @f$.
  TEST(JacobiPolynomial, DerivativeIdentity_K6)
  {
    // d/dx P_6^{\alpha,\beta}(x) = (6 + \alpha + \beta + 1)/2 * P_5^{\alpha+1,\beta+1}(x)
    Real P, dP;
    Real x = 0.4;

    JacobiPolynomial<6>::getValue(P, dP, 0.0, 0.0, x);

    Real P5_11, dP5_11;
    JacobiPolynomial<5>::getValue(P5_11, dP5_11, 1.0, 1.0, x);

    Real expected_deriv = 3.5 * P5_11;  // (6 + 0 + 0 + 1)/2 = 3.5
    EXPECT_NEAR(dP, expected_deriv, 1e-11);
  }

  /// @brief Verify consistency at multiple points for K=6.
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

    // P_6(0) for Legendre (\alpha=\beta=0)
    JacobiPolynomial<6>::getValue(P, dP, 0.0, 0.0, 0.0);
    Real expected = -5.0 / 16.0;  // Known value for P_6(0)
    EXPECT_NEAR(P, expected, tol);
  }

  //==========================================================================
  // Very High Order Tests (K = 15)
  //==========================================================================

  /// @brief P_15^{0,0}(1) = 1, P_15^{0,0}(-1) = -1 (Legendre, odd degree).
  TEST(JacobiPolynomial, P00_K15_EndpointValues)
  {
    // P_15^{0,0}(1) = 1, P_15^{0,0}(-1) = -1 (Legendre, odd degree)
    Real P, dP;

    JacobiPolynomial<15>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-11);

    JacobiPolynomial<15>::getValue(P, dP, 0.0, 0.0, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-11);  // Odd polynomial
  }

  /// @brief P_15(0) = 0 for odd degree Legendre polynomial.
  TEST(JacobiPolynomial, P00_K15_AtZero)
  {
    // P_15(0) = 0 for odd degree Legendre polynomial
    Real P, dP;

    JacobiPolynomial<15>::getValue(P, dP, 0.0, 0.0, 0.0);
    EXPECT_NEAR(P, 0.0, 1e-12);  // Odd polynomial
  }

  /// @brief Test P_15^{10,0}(x) - used in high-order Dubiner tetrahedron.
  TEST(JacobiPolynomial, HigherAlpha_K15)
  {
    // Test P_15^{10,0}(x) - used in high-order Dubiner tetrahedron
    Real P, dP;

    JacobiPolynomial<15>::getValue(P, dP, 10.0, 0.0, 0.0);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));

    JacobiPolynomial<15>::getValue(P, dP, 10.0, 0.0, 0.5);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));

    JacobiPolynomial<15>::getValue(P, dP, 10.0, 0.0, 1.0);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));
  }

  /// @brief Test P_15^{20,0}(x) - verifying stability for very high alpha.
  TEST(JacobiPolynomial, VeryHighAlpha_K15)
  {
    // Test P_15^{20,0}(x) - verifying stability for very high alpha
    Real P, dP;

    JacobiPolynomial<15>::getValue(P, dP, 20.0, 0.0, 0.0);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));

    JacobiPolynomial<15>::getValue(P, dP, 20.0, 0.0, 0.5);
    EXPECT_FALSE(std::isnan(P));
    EXPECT_FALSE(std::isinf(P));
  }

  /// @brief Verify consistency at multiple points for K=15.
  TEST(JacobiPolynomial, RecurrenceConsistency_K15)
  {
    // Verify consistency at multiple points for K=15
    Real P, dP;

    // Test at several points
    std::vector<Real> test_points = {-0.9, -0.5, 0.0, 0.5, 0.9};

    for (Real x : test_points)
    {
      JacobiPolynomial<15>::getValue(P, dP, 0.0, 0.0, x);
      // Just verify no NaN/Inf
      EXPECT_FALSE(std::isnan(P));
      EXPECT_FALSE(std::isinf(P));
      EXPECT_FALSE(std::isnan(dP));
      EXPECT_FALSE(std::isinf(dP));
    }

    JacobiPolynomial<15>::getValue(P, dP, 0.0, 0.0, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-11);

    JacobiPolynomial<15>::getValue(P, dP, 0.0, 0.0, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-11);
  }

  /// @brief Verifies @f$ \frac{d}{dx}P_{15}^{\alpha,\beta}(x) = \frac{16+\alpha+\beta}{2} P_{14}^{\alpha+1,\beta+1}(x) @f$.
  TEST(JacobiPolynomial, DerivativeIdentity_K15)
  {
    // d/dx P_15^{\alpha,\beta}(x) = (15 + \alpha + \beta + 1)/2 * P_14^{\alpha+1,\beta+1}(x)
    Real P, dP;
    Real x = 0.3;

    JacobiPolynomial<15>::getValue(P, dP, 0.0, 0.0, x);

    Real P14_11, dP14_11;
    JacobiPolynomial<14>::getValue(P14_11, dP14_11, 1.0, 1.0, x);

    Real expected_deriv = 8.0 * P14_11;  // (15 + 0 + 0 + 1)/2 = 8
    EXPECT_NEAR(dP, expected_deriv, 1e-9);
  }
}
