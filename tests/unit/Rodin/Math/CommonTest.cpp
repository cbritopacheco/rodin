/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/Common.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

class CommonTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test abs function for various types
TEST_F(CommonTest, AbsFunction)
{
  // Test abs for Real numbers
  EXPECT_DOUBLE_EQ(abs(5.0), 5.0);
  EXPECT_DOUBLE_EQ(abs(-5.0), 5.0);
  EXPECT_DOUBLE_EQ(abs(0.0), 0.0);
  EXPECT_DOUBLE_EQ(abs(-3.14), 3.14);

  // Test abs for integers
  EXPECT_EQ(abs(10), 10);
  EXPECT_EQ(abs(-10), 10);
  EXPECT_EQ(abs(0), 0);

  // Test abs for float
  EXPECT_FLOAT_EQ(abs(2.5f), 2.5f);
  EXPECT_FLOAT_EQ(abs(-2.5f), 2.5f);
}

// Test exp function
TEST_F(CommonTest, ExpFunction)
{
  // Test exp for various values
  EXPECT_DOUBLE_EQ(exp(0.0), 1.0);
  EXPECT_NEAR(exp(1.0), M_E, 1e-10);
  EXPECT_NEAR(exp(2.0), M_E * M_E, 1e-10);
  EXPECT_NEAR(exp(-1.0), 1.0 / M_E, 1e-10);

  // Test exp for float
  EXPECT_FLOAT_EQ(exp(0.0f), 1.0f);
  EXPECT_NEAR(exp(1.0f), static_cast<float>(M_E), 1e-6f);
}

// Test conj function for complex numbers
TEST_F(CommonTest, ConjFunction)
{
  // Test complex conjugate
  Complex z1(3.0, 4.0);  // 3 + 4i
  Complex conj_z1 = conj(z1);
  EXPECT_DOUBLE_EQ(conj_z1.real(), 3.0);
  EXPECT_DOUBLE_EQ(conj_z1.imag(), -4.0);

  // Test conjugate of real number (imaginary part is 0)
  Complex z2(5.0, 0.0);  // 5 + 0i
  Complex conj_z2 = conj(z2);
  EXPECT_DOUBLE_EQ(conj_z2.real(), 5.0);
  EXPECT_DOUBLE_EQ(conj_z2.imag(), 0.0);

  // Test conjugate of purely imaginary number
  Complex z3(0.0, -2.0);  // 0 - 2i
  Complex conj_z3 = conj(z3);
  EXPECT_DOUBLE_EQ(conj_z3.real(), 0.0);
  EXPECT_DOUBLE_EQ(conj_z3.imag(), 2.0);

  // Test conjugate of zero
  Complex z4(0.0, 0.0);
  Complex conj_z4 = conj(z4);
  EXPECT_DOUBLE_EQ(conj_z4.real(), 0.0);
  EXPECT_DOUBLE_EQ(conj_z4.imag(), 0.0);
}

// Test conj function for Eigen matrices
TEST_F(CommonTest, ConjMatrixFunction)
{
  // Test conjugate of real matrix (should be unchanged)
  Eigen::Matrix<Real, 2, 2> realMat;
  realMat << 1.0, 2.0,
             3.0, 4.0;

  auto conjRealMat = conj(realMat);
  EXPECT_DOUBLE_EQ(conjRealMat(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(conjRealMat(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(conjRealMat(1, 0), 3.0);
  EXPECT_DOUBLE_EQ(conjRealMat(1, 1), 4.0);

  // Test conjugate of complex matrix
  Eigen::Matrix<Complex, 2, 2> complexMat;
  complexMat << Complex(1.0, 2.0), Complex(3.0, -1.0),
                Complex(-2.0, 1.0), Complex(0.0, 3.0);

  auto conjComplexMat = conj(complexMat);
  EXPECT_DOUBLE_EQ(conjComplexMat(0, 0).real(), 1.0);
  EXPECT_DOUBLE_EQ(conjComplexMat(0, 0).imag(), -2.0);
  EXPECT_DOUBLE_EQ(conjComplexMat(0, 1).real(), 3.0);
  EXPECT_DOUBLE_EQ(conjComplexMat(0, 1).imag(), 1.0);
  EXPECT_DOUBLE_EQ(conjComplexMat(1, 0).real(), -2.0);
  EXPECT_DOUBLE_EQ(conjComplexMat(1, 0).imag(), -1.0);
  EXPECT_DOUBLE_EQ(conjComplexMat(1, 1).real(), 0.0);
  EXPECT_DOUBLE_EQ(conjComplexMat(1, 1).imag(), -3.0);
}

// Test conj function for Eigen vectors
TEST_F(CommonTest, ConjVectorFunction)
{
  // Test conjugate of real vector
  Eigen::Vector<Real, 3> realVec;
  realVec << 1.0, 2.0, 3.0;

  auto conjRealVec = conj(realVec);
  EXPECT_DOUBLE_EQ(conjRealVec[0], 1.0);
  EXPECT_DOUBLE_EQ(conjRealVec[1], 2.0);
  EXPECT_DOUBLE_EQ(conjRealVec[2], 3.0);

  // Test conjugate of complex vector
  Eigen::Vector<Complex, 3> complexVec;
  complexVec << Complex(1.0, 2.0), Complex(3.0, -1.0), Complex(-2.0, 1.0);

  auto conjComplexVec = conj(complexVec);
  EXPECT_DOUBLE_EQ(conjComplexVec[0].real(), 1.0);
  EXPECT_DOUBLE_EQ(conjComplexVec[0].imag(), -2.0);
  EXPECT_DOUBLE_EQ(conjComplexVec[1].real(), 3.0);
  EXPECT_DOUBLE_EQ(conjComplexVec[1].imag(), 1.0);
  EXPECT_DOUBLE_EQ(conjComplexVec[2].real(), -2.0);
  EXPECT_DOUBLE_EQ(conjComplexVec[2].imag(), -1.0);
}

// Test edge cases and special values
TEST_F(CommonTest, EdgeCases)
{
  // Test with very small numbers
  EXPECT_NEAR(abs(1e-10), 1e-10, 1e-15);
  EXPECT_NEAR(abs(-1e-10), 1e-10, 1e-15);

  // Test with very large numbers
  EXPECT_DOUBLE_EQ(abs(1e10), 1e10);
  EXPECT_DOUBLE_EQ(abs(-1e10), 1e10);

  // Test exp with small values
  EXPECT_NEAR(exp(1e-10), 1.0 + 1e-10, 1e-15);

  // Test complex number with very small imaginary part
  Complex z(5.0, 1e-15);
  Complex conj_z = conj(z);
  EXPECT_DOUBLE_EQ(conj_z.real(), 5.0);
  EXPECT_DOUBLE_EQ(conj_z.imag(), -1e-15);
}

// Test const nature of functions where applicable
TEST_F(CommonTest, constFunctions)
{
  // Test that abs is const
  const Real abs_val = abs(-5.0);
  EXPECT_DOUBLE_EQ(abs_val, 5.0);

  // Test that exp is const
  const Real exp_val = exp(0.0);
  EXPECT_DOUBLE_EQ(exp_val, 1.0);

  // Test that conj is const for complex numbers
  const Complex z(3.0, 4.0);
  const Complex conj_z = conj(z);
  EXPECT_DOUBLE_EQ(conj_z.real(), 3.0);
  EXPECT_DOUBLE_EQ(conj_z.imag(), -4.0);
}

// Test mathematical properties
TEST_F(CommonTest, MathematicalProperties)
{
  // Test that |z|^2 = z * conj(z) for complex numbers
  Complex z(3.0, 4.0);
  Complex conj_z = conj(z);
  Real magnitude_squared = std::norm(z);  // |z|^2
  Complex product = z * conj_z;
  EXPECT_NEAR(magnitude_squared, product.real(), 1e-10);
  EXPECT_NEAR(0.0, product.imag(), 1e-10);

  // Test that exp(a + b) = exp(a) * exp(b)
  Real a = 1.0;
  Real b = 2.0;
  EXPECT_NEAR(exp(a + b), exp(a) * exp(b), 1e-10);

  // Test that exp(0) = 1
  EXPECT_DOUBLE_EQ(exp(0.0), 1.0);

  // Test that abs(-x) = abs(x)
  Real x = 3.14;
  EXPECT_DOUBLE_EQ(abs(-x), abs(x));
}

// Test function composition and chaining
TEST_F(CommonTest, FunctionComposition)
{
  // Test abs(exp(x)) for various x
  EXPECT_DOUBLE_EQ(abs(exp(1.0)), exp(1.0));  // exp is always positive
  EXPECT_DOUBLE_EQ(abs(exp(-1.0)), exp(-1.0));
  EXPECT_DOUBLE_EQ(abs(exp(0.0)), 1.0);

  // Test conj(conj(z)) = z
  Complex z(2.0, -3.0);
  Complex double_conj = conj(conj(z));
  EXPECT_DOUBLE_EQ(double_conj.real(), z.real());
  EXPECT_DOUBLE_EQ(double_conj.imag(), z.imag());

  // Test abs(conj(z)) = abs(z)
  EXPECT_DOUBLE_EQ(abs(conj(z)), abs(z));
}
