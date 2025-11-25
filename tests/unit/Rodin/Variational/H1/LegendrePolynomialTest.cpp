/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational/H1/LegendrePolynomial.h"

using namespace Rodin;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  //==========================================================================
  // Legendre Polynomial Value Tests
  //==========================================================================

  TEST(LegendrePolynomial, P0_Value)
  {
    // P_0(x) = 1 for all x
    Real P, dP;

    LegendrePolynomial<0>::getValue(P, dP, 0.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    LegendrePolynomial<0>::getValue(P, dP, 0.5);
    EXPECT_NEAR(P, 1.0, 1e-14);

    LegendrePolynomial<0>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    LegendrePolynomial<0>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);
  }

  TEST(LegendrePolynomial, P1_Value)
  {
    // P_1(x) = x
    Real P, dP;

    LegendrePolynomial<1>::getValue(P, dP, 0.0);
    EXPECT_NEAR(P, 0.0, 1e-14);

    LegendrePolynomial<1>::getValue(P, dP, 0.5);
    EXPECT_NEAR(P, 0.5, 1e-14);

    LegendrePolynomial<1>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-14);

    LegendrePolynomial<1>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);
  }

  TEST(LegendrePolynomial, P2_Value)
  {
    // P_2(x) = (3x^2 - 1) / 2
    Real P, dP;

    LegendrePolynomial<2>::getValue(P, dP, 0.0);
    EXPECT_NEAR(P, -0.5, 1e-14);  // (3*0 - 1)/2 = -0.5

    LegendrePolynomial<2>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);   // (3*1 - 1)/2 = 1

    LegendrePolynomial<2>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);   // (3*1 - 1)/2 = 1

    Real x = 0.5;
    Real expected = (3.0 * x * x - 1.0) / 2.0;
    LegendrePolynomial<2>::getValue(P, dP, x);
    EXPECT_NEAR(P, expected, 1e-14);
  }

  TEST(LegendrePolynomial, P3_Value)
  {
    // P_3(x) = (5x^3 - 3x) / 2
    Real P, dP;

    LegendrePolynomial<3>::getValue(P, dP, 0.0);
    EXPECT_NEAR(P, 0.0, 1e-14);

    LegendrePolynomial<3>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);  // (5 - 3)/2 = 1

    LegendrePolynomial<3>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-14);

    Real x = 0.5;
    Real expected = (5.0 * x * x * x - 3.0 * x) / 2.0;
    LegendrePolynomial<3>::getValue(P, dP, x);
    EXPECT_NEAR(P, expected, 1e-14);
  }

  TEST(LegendrePolynomial, P4_Value)
  {
    // P_4(x) = (35x^4 - 30x^2 + 3) / 8
    Real P, dP;

    LegendrePolynomial<4>::getValue(P, dP, 0.0);
    EXPECT_NEAR(P, 3.0/8.0, 1e-14);

    LegendrePolynomial<4>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    LegendrePolynomial<4>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);
  }

  //==========================================================================
  // Legendre Polynomial Derivative Tests
  //==========================================================================

  TEST(LegendrePolynomial, P0_Derivative)
  {
    // P'_0(x) = 0
    Real P, dP;

    LegendrePolynomial<0>::getValue(P, dP, 0.0);
    EXPECT_NEAR(dP, 0.0, 1e-14);

    LegendrePolynomial<0>::getValue(P, dP, 0.5);
    EXPECT_NEAR(dP, 0.0, 1e-14);
  }

  TEST(LegendrePolynomial, P1_Derivative)
  {
    // P'_1(x) = 1
    Real P, dP;

    LegendrePolynomial<1>::getValue(P, dP, 0.0);
    EXPECT_NEAR(dP, 1.0, 1e-14);

    LegendrePolynomial<1>::getValue(P, dP, 0.5);
    EXPECT_NEAR(dP, 1.0, 1e-14);

    LegendrePolynomial<1>::getValue(P, dP, -1.0);
    EXPECT_NEAR(dP, 1.0, 1e-14);
  }

  TEST(LegendrePolynomial, P2_Derivative)
  {
    // P'_2(x) = 3x
    Real P, dP;

    LegendrePolynomial<2>::getValue(P, dP, 0.0);
    EXPECT_NEAR(dP, 0.0, 1e-14);

    LegendrePolynomial<2>::getValue(P, dP, 0.5);
    EXPECT_NEAR(dP, 1.5, 1e-14);

    LegendrePolynomial<2>::getValue(P, dP, 1.0);
    EXPECT_NEAR(dP, 3.0, 1e-14);

    LegendrePolynomial<2>::getValue(P, dP, -1.0);
    EXPECT_NEAR(dP, -3.0, 1e-14);
  }

  TEST(LegendrePolynomial, P3_Derivative)
  {
    // P'_3(x) = (15x^2 - 3) / 2
    Real P, dP;

    LegendrePolynomial<3>::getValue(P, dP, 0.0);
    EXPECT_NEAR(dP, -1.5, 1e-14);  // (0 - 3)/2 = -1.5

    Real x = 0.5;
    Real expected = (15.0 * x * x - 3.0) / 2.0;
    LegendrePolynomial<3>::getValue(P, dP, x);
    EXPECT_NEAR(dP, expected, 1e-13);
  }

  //==========================================================================
  // Legendre Polynomial Endpoint Tests
  //==========================================================================

  TEST(LegendrePolynomial, EndpointValues_K0_to_K10)
  {
    // P_K(1) = 1 for all K
    // P_K(-1) = (-1)^K for all K
    Real P, dP;

    LegendrePolynomial<0>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);
    LegendrePolynomial<0>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    LegendrePolynomial<1>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);
    LegendrePolynomial<1>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-14);

    LegendrePolynomial<2>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);
    LegendrePolynomial<2>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    LegendrePolynomial<5>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);
    LegendrePolynomial<5>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-14);

    LegendrePolynomial<10>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-13);
    LegendrePolynomial<10>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-13);
  }

  //==========================================================================
  // Legendre Polynomial Orthogonality (Sanity Check)
  //==========================================================================

  TEST(LegendrePolynomial, Orthogonality_P0_P1)
  {
    // ∫_{-1}^{1} P_0(x) P_1(x) dx = 0
    // Numerical integration via Gauss-Legendre quadrature
    const int N = 20;
    Real sum = 0.0;
    for (int i = 0; i < N; ++i)
    {
      Real x = -1.0 + 2.0 * (i + 0.5) / N;
      Real P0, dP0, P1, dP1;
      LegendrePolynomial<0>::getValue(P0, dP0, x);
      LegendrePolynomial<1>::getValue(P1, dP1, x);
      sum += P0 * P1 * (2.0 / N);
    }
    EXPECT_NEAR(sum, 0.0, 1e-12);
  }

  TEST(LegendrePolynomial, Orthogonality_P1_P2)
  {
    // ∫_{-1}^{1} P_1(x) P_2(x) dx = 0
    const int N = 20;
    Real sum = 0.0;
    for (int i = 0; i < N; ++i)
    {
      Real x = -1.0 + 2.0 * (i + 0.5) / N;
      Real P1, dP1, P2, dP2;
      LegendrePolynomial<1>::getValue(P1, dP1, x);
      LegendrePolynomial<2>::getValue(P2, dP2, x);
      sum += P1 * P2 * (2.0 / N);
    }
    EXPECT_NEAR(sum, 0.0, 1e-12);
  }

  //==========================================================================
  // Legendre Polynomial Symmetry Tests
  //==========================================================================

  TEST(LegendrePolynomial, Symmetry_EvenDegree)
  {
    // P_K(-x) = P_K(x) for even K
    Real P_pos, dP_pos, P_neg, dP_neg;
    Real x = 0.3;

    LegendrePolynomial<2>::getValue(P_pos, dP_pos, x);
    LegendrePolynomial<2>::getValue(P_neg, dP_neg, -x);
    EXPECT_NEAR(P_pos, P_neg, 1e-14);

    LegendrePolynomial<4>::getValue(P_pos, dP_pos, x);
    LegendrePolynomial<4>::getValue(P_neg, dP_neg, -x);
    EXPECT_NEAR(P_pos, P_neg, 1e-14);

    LegendrePolynomial<6>::getValue(P_pos, dP_pos, x);
    LegendrePolynomial<6>::getValue(P_neg, dP_neg, -x);
    EXPECT_NEAR(P_pos, P_neg, 1e-14);
  }

  TEST(LegendrePolynomial, Symmetry_OddDegree)
  {
    // P_K(-x) = -P_K(x) for odd K
    Real P_pos, dP_pos, P_neg, dP_neg;
    Real x = 0.3;

    LegendrePolynomial<1>::getValue(P_pos, dP_pos, x);
    LegendrePolynomial<1>::getValue(P_neg, dP_neg, -x);
    EXPECT_NEAR(P_pos, -P_neg, 1e-14);

    LegendrePolynomial<3>::getValue(P_pos, dP_pos, x);
    LegendrePolynomial<3>::getValue(P_neg, dP_neg, -x);
    EXPECT_NEAR(P_pos, -P_neg, 1e-14);

    LegendrePolynomial<5>::getValue(P_pos, dP_pos, x);
    LegendrePolynomial<5>::getValue(P_neg, dP_neg, -x);
    EXPECT_NEAR(P_pos, -P_neg, 1e-14);
  }

  //==========================================================================
  // Higher Order Tests (K = 5, 6)
  //==========================================================================

  TEST(LegendrePolynomial, P5_Value)
  {
    // P_5(x) = (63x^5 - 70x^3 + 15x) / 8
    Real P, dP;

    LegendrePolynomial<5>::getValue(P, dP, 0.0);
    EXPECT_NEAR(P, 0.0, 1e-14);  // Odd polynomial

    LegendrePolynomial<5>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    LegendrePolynomial<5>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-14);

    Real x = 0.5;
    Real expected = (63.0 * std::pow(x, 5) - 70.0 * std::pow(x, 3) + 15.0 * x) / 8.0;
    LegendrePolynomial<5>::getValue(P, dP, x);
    EXPECT_NEAR(P, expected, 1e-14);
  }

  TEST(LegendrePolynomial, P6_Value)
  {
    // P_6(x) = (231x^6 - 315x^4 + 105x^2 - 5) / 16
    Real P, dP;

    LegendrePolynomial<6>::getValue(P, dP, 0.0);
    Real expected_at_0 = -5.0 / 16.0;
    EXPECT_NEAR(P, expected_at_0, 1e-14);

    LegendrePolynomial<6>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);

    LegendrePolynomial<6>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, 1.0, 1e-14);  // Even polynomial

    Real x = 0.5;
    Real expected = (231.0 * std::pow(x, 6) - 315.0 * std::pow(x, 4) + 105.0 * x * x - 5.0) / 16.0;
    LegendrePolynomial<6>::getValue(P, dP, x);
    EXPECT_NEAR(P, expected, 1e-13);
  }

  TEST(LegendrePolynomial, P5_Derivative)
  {
    // P'_5(x) = (315x^4 - 210x^2 + 15) / 8
    Real P, dP;

    Real x = 0.5;
    Real expected = (315.0 * std::pow(x, 4) - 210.0 * x * x + 15.0) / 8.0;
    LegendrePolynomial<5>::getValue(P, dP, x);
    EXPECT_NEAR(dP, expected, 1e-13);

    LegendrePolynomial<5>::getValue(P, dP, 0.0);
    EXPECT_NEAR(dP, 15.0 / 8.0, 1e-14);
  }

  TEST(LegendrePolynomial, P6_Derivative)
  {
    // P'_6(x) = (1386x^5 - 1260x^3 + 210x) / 16
    Real P, dP;

    Real x = 0.5;
    Real expected = (1386.0 * std::pow(x, 5) - 1260.0 * std::pow(x, 3) + 210.0 * x) / 16.0;
    LegendrePolynomial<6>::getValue(P, dP, x);
    EXPECT_NEAR(dP, expected, 1e-12);

    LegendrePolynomial<6>::getValue(P, dP, 0.0);
    EXPECT_NEAR(dP, 0.0, 1e-14);  // Even polynomial has odd derivative
  }

  TEST(LegendrePolynomial, Orthogonality_P2_P5)
  {
    // ∫_{-1}^{1} P_2(x) P_5(x) dx = 0
    const int N = 30;
    Real sum = 0.0;
    for (int i = 0; i < N; ++i)
    {
      Real x = -1.0 + 2.0 * (i + 0.5) / N;
      Real P2, dP2, P5, dP5;
      LegendrePolynomial<2>::getValue(P2, dP2, x);
      LegendrePolynomial<5>::getValue(P5, dP5, x);
      sum += P2 * P5 * (2.0 / N);
    }
    EXPECT_NEAR(sum, 0.0, 1e-12);
  }

  //==========================================================================
  // Very High Order Tests (K = 15)
  //==========================================================================

  TEST(LegendrePolynomial, P15_EndpointValues)
  {
    // P_15(1) = 1, P_15(-1) = -1 (odd degree)
    Real P, dP;

    LegendrePolynomial<15>::getValue(P, dP, 1.0);
    EXPECT_NEAR(P, 1.0, 1e-11);

    LegendrePolynomial<15>::getValue(P, dP, -1.0);
    EXPECT_NEAR(P, -1.0, 1e-11);
  }

  TEST(LegendrePolynomial, P15_AtZero)
  {
    // P_15(0) = 0 (odd polynomial)
    Real P, dP;

    LegendrePolynomial<15>::getValue(P, dP, 0.0);
    EXPECT_NEAR(P, 0.0, 1e-12);
  }

  TEST(LegendrePolynomial, P15_Symmetry_OddDegree)
  {
    // P_15(-x) = -P_15(x) for odd K
    Real P_pos, dP_pos, P_neg, dP_neg;

    std::vector<Real> test_points = {0.1, 0.3, 0.5, 0.7, 0.9};

    for (Real x : test_points)
    {
      LegendrePolynomial<15>::getValue(P_pos, dP_pos, x);
      LegendrePolynomial<15>::getValue(P_neg, dP_neg, -x);
      EXPECT_NEAR(P_pos, -P_neg, 1e-12);
    }
  }

  TEST(LegendrePolynomial, P15_Value_NoNaN)
  {
    // Verify no NaN/Inf at various points
    Real P, dP;

    std::vector<Real> test_points = {-0.9, -0.5, 0.0, 0.5, 0.9};

    for (Real x : test_points)
    {
      LegendrePolynomial<15>::getValue(P, dP, x);
      EXPECT_FALSE(std::isnan(P));
      EXPECT_FALSE(std::isinf(P));
      EXPECT_FALSE(std::isnan(dP));
      EXPECT_FALSE(std::isinf(dP));
    }
  }

  TEST(LegendrePolynomial, P15_Derivative_NoNaN)
  {
    Real P, dP;

    std::vector<Real> test_points = {-0.9, -0.5, 0.0, 0.5, 0.9};

    for (Real x : test_points)
    {
      LegendrePolynomial<15>::getValue(P, dP, x);
      EXPECT_FALSE(std::isnan(dP));
      EXPECT_FALSE(std::isinf(dP));
    }
  }

  TEST(LegendrePolynomial, Orthogonality_P5_P15)
  {
    // ∫_{-1}^{1} P_5(x) P_15(x) dx = 0
    // Use Gauss-Legendre quadrature with 11 points (exact for deg <= 21)
    // Nodes and weights for 11-point Gauss-Legendre quadrature on [-1, 1]
    const std::array<Real, 11> nodes = {
      -0.9782286581460570,
      -0.8870625997680953,
      -0.7301520055740494,
      -0.5190961292068118,
      -0.2695431559523450,
      0.0,
      0.2695431559523450,
      0.5190961292068118,
      0.7301520055740494,
      0.8870625997680953,
      0.9782286581460570
    };
    const std::array<Real, 11> weights = {
      0.0556685671161737,
      0.1255803694649046,
      0.1862902109277343,
      0.2331937645919905,
      0.2628045445102467,
      0.2729250867779006,
      0.2628045445102467,
      0.2331937645919905,
      0.1862902109277343,
      0.1255803694649046,
      0.0556685671161737
    };

    Real sum = 0.0;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
      Real P5, dP5, P15, dP15;
      LegendrePolynomial<5>::getValue(P5, dP5, nodes[i]);
      LegendrePolynomial<15>::getValue(P15, dP15, nodes[i]);
      sum += weights[i] * P5 * P15;
    }
    EXPECT_NEAR(sum, 0.0, 1e-12);
  }
}
