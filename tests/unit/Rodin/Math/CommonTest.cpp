/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/Common.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialMatrix.h"
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

TEST_F(CommonTest, Pow2Function)
{
  EXPECT_DOUBLE_EQ(pow2(3.0), 9.0);
  EXPECT_DOUBLE_EQ(pow2(-4.0), 16.0);
  EXPECT_DOUBLE_EQ(pow2(0.0), 0.0);
}

TEST_F(CommonTest, PowNFunction)
{
  EXPECT_DOUBLE_EQ((pow<0>(5.0)), 1.0);
  EXPECT_DOUBLE_EQ((pow<1>(3.0)), 3.0);
  EXPECT_DOUBLE_EQ((pow<2>(3.0)), 9.0);
  EXPECT_DOUBLE_EQ((pow<3>(2.0)), 8.0);
  EXPECT_DOUBLE_EQ((pow<4>(2.0)), 16.0);
  EXPECT_DOUBLE_EQ((pow<10>(2.0)), 1024.0);
}

TEST_F(CommonTest, PowRuntimeFunction)
{
  EXPECT_DOUBLE_EQ(pow(2.0, 3.0), 8.0);
  EXPECT_DOUBLE_EQ(pow(4.0, 0.5), 2.0);
  EXPECT_DOUBLE_EQ(pow(2.0, 0.0), 1.0);
}

TEST_F(CommonTest, SqrtFunction)
{
  EXPECT_DOUBLE_EQ(sqrt(4.0), 2.0);
  EXPECT_DOUBLE_EQ(sqrt(9.0), 3.0);
  EXPECT_DOUBLE_EQ(sqrt(0.0), 0.0);
  EXPECT_DOUBLE_EQ(sqrt(1.0), 1.0);
  EXPECT_NEAR(sqrt(2.0), 1.41421356, 1e-6);
}

TEST_F(CommonTest, IsNaNFunction)
{
  EXPECT_TRUE(isNaN(nan<Real>()));
  EXPECT_FALSE(isNaN(1.0));
  EXPECT_FALSE(isNaN(0.0));

  EXPECT_TRUE(isNaN(nan<Complex>()));
  EXPECT_TRUE(isNaN(Complex(1.0, nan<Real>())));
  EXPECT_FALSE(isNaN(Complex(1.0, 2.0)));
}

TEST_F(CommonTest, IsInfFunction)
{
  EXPECT_TRUE(isInf(std::numeric_limits<Real>::infinity()));
  EXPECT_TRUE(isInf(-std::numeric_limits<Real>::infinity()));
  EXPECT_FALSE(isInf(1.0));
  EXPECT_FALSE(isInf(0.0));
}

TEST_F(CommonTest, TrigFunctions)
{
  EXPECT_DOUBLE_EQ(cos(0.0), 1.0);
  EXPECT_DOUBLE_EQ(sin(0.0), 0.0);
  EXPECT_DOUBLE_EQ(tan(0.0), 0.0);
  EXPECT_NEAR(cos(M_PI), -1.0, 1e-10);
  EXPECT_NEAR(sin(M_PI / 2.0), 1.0, 1e-10);
  EXPECT_NEAR(tan(M_PI / 4.0), 1.0, 1e-10);
}

TEST_F(CommonTest, HyperbolicFunctions)
{
  EXPECT_DOUBLE_EQ(cosh(0.0), 1.0);
  EXPECT_DOUBLE_EQ(sinh(0.0), 0.0);
  EXPECT_DOUBLE_EQ(tanh(0.0), 0.0);
  EXPECT_NEAR(cosh(1.0), 1.5430806, 1e-6);
  EXPECT_NEAR(sinh(1.0), 1.1752012, 1e-6);
  EXPECT_NEAR(tanh(1.0), 0.7615942, 1e-6);
}

TEST_F(CommonTest, InverseTrigFunctions)
{
  EXPECT_DOUBLE_EQ(acos(1.0), 0.0);
  EXPECT_DOUBLE_EQ(asin(0.0), 0.0);
  EXPECT_DOUBLE_EQ(atan(0.0), 0.0);
  EXPECT_NEAR(acos(0.0), M_PI / 2.0, 1e-10);
  EXPECT_NEAR(asin(1.0), M_PI / 2.0, 1e-10);
  EXPECT_NEAR(atan(1.0), M_PI / 4.0, 1e-10);
  EXPECT_NEAR(atan2(1.0, 1.0), M_PI / 4.0, 1e-10);
  EXPECT_DOUBLE_EQ(atan2(0.0, 1.0), 0.0);
}

TEST_F(CommonTest, LogFunctions)
{
  EXPECT_DOUBLE_EQ(log(1.0), 0.0);
  EXPECT_NEAR(log(M_E), 1.0, 1e-10);
  EXPECT_DOUBLE_EQ(log2(1.0), 0.0);
  EXPECT_DOUBLE_EQ(log2(8.0), 3.0);
  EXPECT_DOUBLE_EQ(log10(1.0), 0.0);
  EXPECT_DOUBLE_EQ(log10(100.0), 2.0);
}

TEST_F(CommonTest, SgnFunction)
{
  EXPECT_DOUBLE_EQ(sgn(-5.0), -1.0);
  EXPECT_DOUBLE_EQ(sgn(5.0), 1.0);
  EXPECT_DOUBLE_EQ(sgn(0.0), 0.0);
  EXPECT_EQ(sgn(-3), -1);
  EXPECT_EQ(sgn(3), 1);
  EXPECT_EQ(sgn(0), 0);
}

TEST_F(CommonTest, BinomFunction)
{
  EXPECT_EQ(binom(Integer(5), Integer(2)), 10);
  EXPECT_EQ(binom(Integer(10), Integer(3)), 120);
  EXPECT_EQ(binom(Integer(0), Integer(0)), 1);
  EXPECT_EQ(binom(Integer(5), Integer(0)), 1);
  EXPECT_EQ(binom(Integer(5), Integer(5)), 1);
}

TEST_F(CommonTest, FactorialFunction)
{
  EXPECT_EQ(factorial(Integer(0)), 1);
  EXPECT_EQ(factorial(Integer(1)), 1);
  EXPECT_EQ(factorial(Integer(5)), 120);
  EXPECT_EQ(factorial(Integer(10)), 3628800);
}

TEST_F(CommonTest, PermutationFunction)
{
  EXPECT_EQ(permutation(Integer(5), Integer(2)), 20);
  EXPECT_EQ(permutation(Integer(5), Integer(0)), 1);
  EXPECT_EQ(permutation(Integer(5), Integer(5)), 120);
}

TEST_F(CommonTest, NanFactory)
{
  EXPECT_TRUE(isNaN(nan<Real>()));

  auto cn = nan<Complex>();
  EXPECT_TRUE(isNaN(cn));
  EXPECT_TRUE(std::isnan(cn.real()));
  EXPECT_TRUE(std::isnan(cn.imag()));
}

TEST_F(CommonTest, FormLanguageHelpers)
{
  EXPECT_DOUBLE_EQ(sum(2.0, 3.0), 5.0);
  EXPECT_DOUBLE_EQ(minus(5.0), -5.0);
  EXPECT_DOUBLE_EQ(minus(5.0, 3.0), 2.0);
  EXPECT_DOUBLE_EQ(mult(3.0, 4.0), 12.0);
  EXPECT_DOUBLE_EQ(division(10.0, 2.0), 5.0);
}

TEST_F(CommonTest, DotProductScalar)
{
  EXPECT_DOUBLE_EQ(dot(Real(2.0), Real(3.0)), 6.0);

  // dot(z1, z2) = z1 * conj(z2)
  // (1+2i) * conj(3+4i) = (1+2i) * (3-4i) = 3 - 4i + 6i - 8i^2 = 3 + 2i + 8 = 11 + 2i
  Complex result = dot(Complex(1, 2), Complex(3, 4));
  EXPECT_DOUBLE_EQ(result.real(), 11.0);
  EXPECT_DOUBLE_EQ(result.imag(), 2.0);
}

TEST_F(CommonTest, DotProductEigen)
{
  Eigen::Vector3d v1(1, 2, 3);
  Eigen::Vector3d v2(4, 5, 6);
  EXPECT_DOUBLE_EQ(dot(v1, v2), 32.0);
}

TEST_F(CommonTest, MinFunction)
{
  EXPECT_DOUBLE_EQ(min(3.0, 5.0), 3.0);
  EXPECT_DOUBLE_EQ(min(5.0, 3.0), 3.0);
  EXPECT_DOUBLE_EQ(min(-1.0, 1.0), -1.0);
  EXPECT_DOUBLE_EQ(min(0.0, 0.0), 0.0);
  EXPECT_EQ(min(2, 7), 2);
}

TEST_F(CommonTest, MaxFunction)
{
  EXPECT_DOUBLE_EQ(max(3.0, 5.0), 5.0);
  EXPECT_DOUBLE_EQ(max(5.0, 3.0), 5.0);
  EXPECT_DOUBLE_EQ(max(-1.0, 1.0), 1.0);
  EXPECT_DOUBLE_EQ(max(0.0, 0.0), 0.0);
  EXPECT_EQ(max(2, 7), 7);
}

TEST_F(CommonTest, ClampFunction)
{
  EXPECT_DOUBLE_EQ(clamp(5.0, 0.0, 10.0), 5.0);
  EXPECT_DOUBLE_EQ(clamp(-1.0, 0.0, 10.0), 0.0);
  EXPECT_DOUBLE_EQ(clamp(15.0, 0.0, 10.0), 10.0);
  EXPECT_DOUBLE_EQ(clamp(0.0, 0.0, 10.0), 0.0);
  EXPECT_DOUBLE_EQ(clamp(10.0, 0.0, 10.0), 10.0);
  EXPECT_EQ(clamp(3, 1, 5), 3);
  EXPECT_EQ(clamp(0, 1, 5), 1);
  EXPECT_EQ(clamp(7, 1, 5), 5);
}

TEST_F(CommonTest, DotProductSpatialVector)
{
  SpatialVector<Real> a({1.0, 2.0, 3.0});
  SpatialVector<Real> b({4.0, 5.0, 6.0});
  EXPECT_DOUBLE_EQ(dot(a, b), 32.0);
}

TEST_F(CommonTest, DotProductSpatialMatrix)
{
  SpatialMatrix<Real> a(2, 2);
  a(0, 0) = 1; a(0, 1) = 2;
  a(1, 0) = 3; a(1, 1) = 4;

  SpatialMatrix<Real> b(2, 2);
  b(0, 0) = 5; b(0, 1) = 6;
  b(1, 0) = 7; b(1, 1) = 8;

  EXPECT_DOUBLE_EQ(dot(a, b), 70.0);
}

TEST_F(CommonTest, DotProductMixedEigenSpatialVector)
{
  Eigen::Vector3d ev(1, 2, 3);
  SpatialVector<Real> sv({4.0, 5.0, 6.0});

  EXPECT_DOUBLE_EQ(dot(ev, sv), 32.0);
  EXPECT_DOUBLE_EQ(dot(sv, ev), 32.0);
}
