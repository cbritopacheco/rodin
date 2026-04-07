/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

class MatrixTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test basic matrix type aliases
TEST_F(MatrixTest, TypeAliases)
{
  // Test that Matrix<Real> is the same as RealMatrix
  static_assert(std::is_same_v<Matrix<Real>, RealMatrix>);

  // Test that Matrix<Complex> is the same as ComplexMatrix
  static_assert(std::is_same_v<Matrix<Complex>, ComplexMatrix>);
}

// Test matrix construction and basic operations
TEST_F(MatrixTest, Construction)
{
  // Test default construction
  RealMatrix m1;
  EXPECT_EQ(m1.rows(), 0);
  EXPECT_EQ(m1.cols(), 0);

  // Test sized construction
  RealMatrix m2(3, 4);
  EXPECT_EQ(m2.rows(), 3);
  EXPECT_EQ(m2.cols(), 4);

  // Test zero initialization
  RealMatrix m3 = RealMatrix::Zero(2, 3);
  EXPECT_EQ(m3.rows(), 2);
  EXPECT_EQ(m3.cols(), 3);
  for (int i = 0; i < m3.rows(); ++i)
    for (int j = 0; j < m3.cols(); ++j)
      EXPECT_DOUBLE_EQ(m3(i, j), 0.0);
}

// Test matrix assignment and access
TEST_F(MatrixTest, AssignmentAndAccess)
{
  RealMatrix m(2, 2);

  // Test element assignment
  m(0, 0) = 1.0;
  m(0, 1) = 2.0;
  m(1, 0) = 3.0;
  m(1, 1) = 4.0;

  // Test element access
  EXPECT_DOUBLE_EQ(m(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(m(1, 0), 3.0);
  EXPECT_DOUBLE_EQ(m(1, 1), 4.0);
}

// Test matrix arithmetic operations
TEST_F(MatrixTest, ArithmeticOperations)
{
  RealMatrix m1(2, 2);
  m1 << 1, 2,
        3, 4;

  RealMatrix m2(2, 2);
  m2 << 5, 6,
        7, 8;

  // Test matrix addition
  RealMatrix sum = m1 + m2;
  EXPECT_DOUBLE_EQ(sum(0, 0), 6.0);
  EXPECT_DOUBLE_EQ(sum(0, 1), 8.0);
  EXPECT_DOUBLE_EQ(sum(1, 0), 10.0);
  EXPECT_DOUBLE_EQ(sum(1, 1), 12.0);

  // Test matrix subtraction
  RealMatrix diff = m2 - m1;
  EXPECT_DOUBLE_EQ(diff(0, 0), 4.0);
  EXPECT_DOUBLE_EQ(diff(0, 1), 4.0);
  EXPECT_DOUBLE_EQ(diff(1, 0), 4.0);
  EXPECT_DOUBLE_EQ(diff(1, 1), 4.0);

  // Test scalar multiplication
  RealMatrix scaled = 2.0 * m1;
  EXPECT_DOUBLE_EQ(scaled(0, 0), 2.0);
  EXPECT_DOUBLE_EQ(scaled(0, 1), 4.0);
  EXPECT_DOUBLE_EQ(scaled(1, 0), 6.0);
  EXPECT_DOUBLE_EQ(scaled(1, 1), 8.0);
}

// Test matrix multiplication
TEST_F(MatrixTest, MatrixMultiplication)
{
  RealMatrix m1(2, 3);
  m1 << 1, 2, 3,
        4, 5, 6;

  RealMatrix m2(3, 2);
  m2 << 7,  8,
        9,  10,
        11, 12;

  RealMatrix result = m1 * m2;
  EXPECT_EQ(result.rows(), 2);
  EXPECT_EQ(result.cols(), 2);

  // Expected result: [58, 64; 139, 154]
  EXPECT_DOUBLE_EQ(result(0, 0), 58.0);
  EXPECT_DOUBLE_EQ(result(0, 1), 64.0);
  EXPECT_DOUBLE_EQ(result(1, 0), 139.0);
  EXPECT_DOUBLE_EQ(result(1, 1), 154.0);
}

// Test FixedSizeMatrix
TEST_F(MatrixTest, FixedSizeMatrix)
{
  using Matrix2x2 = FixedSizeMatrix<Real, 2, 2>;
  using Matrix3x3 = FixedSizeMatrix<Real, 3, 3>;

  // Test compile-time sizes
  static_assert(Matrix2x2::RowsAtCompileTime == 2);
  static_assert(Matrix2x2::ColsAtCompileTime == 2);
  static_assert(Matrix3x3::RowsAtCompileTime == 3);
  static_assert(Matrix3x3::ColsAtCompileTime == 3);

  // Test construction
  Matrix2x2 m;
  m << 1, 2,
       3, 4;

  EXPECT_EQ(m.rows(), 2);
  EXPECT_EQ(m.cols(), 2);
  EXPECT_DOUBLE_EQ(m(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m(1, 1), 4.0);
}

// Test SpatialMatrix functionality
TEST_F(MatrixTest, SpatialMatrix)
{
  SpatialMatrix<Real> sm(3, 3);
  sm.setIdentity();

  EXPECT_EQ(sm.rows(), 3);
  EXPECT_EQ(sm.cols(), 3);

  // Check identity matrix
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      if (i == j)
        EXPECT_DOUBLE_EQ(sm(i, j), 1.0);
      else
        EXPECT_DOUBLE_EQ(sm(i, j), 0.0);
    }
  }
}

// Test complex matrix operations
TEST_F(MatrixTest, ComplexMatrix)
{
  ComplexMatrix cm(2, 2);
  cm(0, 0) = Complex(1.0, 2.0);  // 1 + 2i
  cm(0, 1) = Complex(3.0, -1.0); // 3 - i
  cm(1, 0) = Complex(-2.0, 1.0); // -2 + i
  cm(1, 1) = Complex(0.0, 3.0);  // 3i

  // Test conjugate transpose
  ComplexMatrix hermitian = cm.adjoint();
  EXPECT_EQ(hermitian.rows(), 2);
  EXPECT_EQ(hermitian.cols(), 2);

  // Check conjugate transpose values
  EXPECT_DOUBLE_EQ(hermitian(0, 0).real(), 1.0);
  EXPECT_DOUBLE_EQ(hermitian(0, 0).imag(), -2.0);
  EXPECT_DOUBLE_EQ(hermitian(1, 0).real(), 3.0);
  EXPECT_DOUBLE_EQ(hermitian(1, 0).imag(), 1.0);
}

// Test matrix resizing
TEST_F(MatrixTest, Resizing)
{
  RealMatrix m(2, 2);
  m << 1, 2,
       3, 4;

  // Test resize
  m.resize(3, 3);
  EXPECT_EQ(m.rows(), 3);
  EXPECT_EQ(m.cols(), 3);

  // Test conservativeResize
  RealMatrix m2(2, 2);
  m2 << 1, 2,
        3, 4;
  m2.conservativeResize(3, 3);
  EXPECT_EQ(m2.rows(), 3);
  EXPECT_EQ(m2.cols(), 3);
  EXPECT_DOUBLE_EQ(m2(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m2(1, 1), 4.0);
}

// Test matrix norms and properties
TEST_F(MatrixTest, MatrixProperties)
{
  RealMatrix m(2, 2);
  m << 3, 4,
       0, 0;

  // Test Frobenius norm
  Real frobNorm = m.norm();
  EXPECT_DOUBLE_EQ(frobNorm, 5.0);  // sqrt(3^2 + 4^2) = 5

  // Test trace
  Real trace = m.trace();
  EXPECT_DOUBLE_EQ(trace, 3.0);  // 3 + 0 = 3

  // Test determinant for square matrix
  RealMatrix square(2, 2);
  square << 1, 2,
            3, 4;
  Real det = square.determinant();
  EXPECT_DOUBLE_EQ(det, -2.0);  // 1*4 - 2*3 = -2
}

// --- Comprehensive SpatialMatrix Tests ---

TEST_F(MatrixTest, SpatialMatrixConstruction)
{
  SpatialMatrix<Real> sm;
  EXPECT_EQ(sm.rows(), 0);
  EXPECT_EQ(sm.cols(), 0);

  SpatialMatrix<Real> sm2(2, 3);
  EXPECT_EQ(sm2.rows(), 2);
  EXPECT_EQ(sm2.cols(), 3);

  auto sm3 = sm2;
  EXPECT_EQ(sm3.rows(), 2);
  EXPECT_EQ(sm3.cols(), 3);
}

TEST_F(MatrixTest, SpatialMatrixElementAccess)
{
  SpatialMatrix<Real> sm(2, 2);
  sm(0, 0) = 1;
  sm(0, 1) = 2;
  sm(1, 0) = 3;
  sm(1, 1) = 4;

  EXPECT_DOUBLE_EQ(sm(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(sm(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(sm(1, 0), 3.0);
  EXPECT_DOUBLE_EQ(sm(1, 1), 4.0);
}

TEST_F(MatrixTest, SpatialMatrixSetZeroConstantIdentity)
{
  SpatialMatrix<Real> sm(3, 3);

  sm.setZero();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_DOUBLE_EQ(sm(i, j), 0.0);

  sm.setConstant(5.0);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_DOUBLE_EQ(sm(i, j), 5.0);

  sm.setIdentity();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      if (i == j)
        EXPECT_DOUBLE_EQ(sm(i, j), 1.0);
      else
        EXPECT_DOUBLE_EQ(sm(i, j), 0.0);
    }
}

TEST_F(MatrixTest, SpatialMatrixResize)
{
  SpatialMatrix<Real> sm(2, 2);
  sm(0, 0) = 1;
  sm(0, 1) = 2;
  sm(1, 0) = 3;
  sm(1, 1) = 4;

  sm.resize(3, 3);
  EXPECT_EQ(sm.rows(), 3);
  EXPECT_EQ(sm.cols(), 3);
}

TEST_F(MatrixTest, SpatialMatrixNorms)
{
  SpatialMatrix<Real> sm(2, 2);
  sm(0, 0) = 3;
  sm(0, 1) = 0;
  sm(1, 0) = 0;
  sm(1, 1) = 4;

  EXPECT_DOUBLE_EQ(sm.squaredNorm(), 25.0);
  EXPECT_DOUBLE_EQ(sm.norm(), 5.0);
}

TEST_F(MatrixTest, SpatialMatrixDot)
{
  SpatialMatrix<Real> a(2, 2);
  a(0, 0) = 1; a(0, 1) = 2;
  a(1, 0) = 3; a(1, 1) = 4;

  SpatialMatrix<Real> b(2, 2);
  b(0, 0) = 5; b(0, 1) = 6;
  b(1, 0) = 7; b(1, 1) = 8;

  EXPECT_DOUBLE_EQ(a.dot(b), 70.0);
}

TEST_F(MatrixTest, SpatialMatrixTranspose)
{
  SpatialMatrix<Real> sm(2, 3);
  sm(0, 0) = 1; sm(0, 1) = 2; sm(0, 2) = 3;
  sm(1, 0) = 4; sm(1, 1) = 5; sm(1, 2) = 6;

  auto t = sm.transpose();
  EXPECT_EQ(t.rows(), 3);
  EXPECT_EQ(t.cols(), 2);
  EXPECT_DOUBLE_EQ(t(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(t(1, 0), 2.0);
  EXPECT_DOUBLE_EQ(t(2, 0), 3.0);
  EXPECT_DOUBLE_EQ(t(0, 1), 4.0);
  EXPECT_DOUBLE_EQ(t(1, 1), 5.0);
  EXPECT_DOUBLE_EQ(t(2, 1), 6.0);
}

TEST_F(MatrixTest, SpatialMatrixTrace)
{
  SpatialMatrix<Real> sm1(1, 1);
  sm1(0, 0) = 5;
  EXPECT_DOUBLE_EQ(sm1.trace(), 5.0);

  SpatialMatrix<Real> sm2(2, 2);
  sm2(0, 0) = 1; sm2(0, 1) = 2;
  sm2(1, 0) = 3; sm2(1, 1) = 4;
  EXPECT_DOUBLE_EQ(sm2.trace(), 5.0);

  SpatialMatrix<Real> sm3(3, 3);
  sm3.setIdentity();
  EXPECT_DOUBLE_EQ(sm3.trace(), 3.0);
}

TEST_F(MatrixTest, SpatialMatrixValue)
{
  SpatialMatrix<Real> sm(1, 1);
  sm(0, 0) = 7.5;
  EXPECT_DOUBLE_EQ(sm.value(), 7.5);
}

TEST_F(MatrixTest, SpatialMatrixDeterminant1x1)
{
  SpatialMatrix<Real> sm(1, 1);
  sm(0, 0) = 5;
  EXPECT_DOUBLE_EQ(sm.determinant(), 5.0);
}

TEST_F(MatrixTest, SpatialMatrixDeterminant2x2)
{
  SpatialMatrix<Real> sm(2, 2);
  sm(0, 0) = 1; sm(0, 1) = 2;
  sm(1, 0) = 3; sm(1, 1) = 4;
  EXPECT_DOUBLE_EQ(sm.determinant(), -2.0);
}

TEST_F(MatrixTest, SpatialMatrixDeterminant3x3)
{
  SpatialMatrix<Real> sm(3, 3);
  sm(0, 0) = 1; sm(0, 1) = 2; sm(0, 2) = 3;
  sm(1, 0) = 0; sm(1, 1) = 1; sm(1, 2) = 4;
  sm(2, 0) = 5; sm(2, 1) = 6; sm(2, 2) = 0;
  EXPECT_NEAR(sm.determinant(), 1.0, 1e-10);
}

TEST_F(MatrixTest, SpatialMatrixInverse1x1)
{
  SpatialMatrix<Real> sm(1, 1);
  sm(0, 0) = 4;
  auto inv = sm.inverse();
  EXPECT_DOUBLE_EQ(inv(0, 0), 0.25);
}

TEST_F(MatrixTest, SpatialMatrixInverse2x2)
{
  SpatialMatrix<Real> sm(2, 2);
  sm(0, 0) = 4; sm(0, 1) = 7;
  sm(1, 0) = 2; sm(1, 1) = 6;
  auto inv = sm.inverse();

  EXPECT_NEAR(inv(0, 0), 0.6, 1e-10);
  EXPECT_NEAR(inv(0, 1), -0.7, 1e-10);
  EXPECT_NEAR(inv(1, 0), -0.2, 1e-10);
  EXPECT_NEAR(inv(1, 1), 0.4, 1e-10);
}

TEST_F(MatrixTest, SpatialMatrixInverse3x3)
{
  SpatialMatrix<Real> sm(3, 3);
  sm(0, 0) = 1; sm(0, 1) = 0; sm(0, 2) = 0;
  sm(1, 0) = 0; sm(1, 1) = 2; sm(1, 2) = 0;
  sm(2, 0) = 0; sm(2, 1) = 0; sm(2, 2) = 3;

  auto inv = sm.inverse();
  EXPECT_NEAR(inv(0, 0), 1.0, 1e-10);
  EXPECT_NEAR(inv(0, 1), 0.0, 1e-10);
  EXPECT_NEAR(inv(0, 2), 0.0, 1e-10);
  EXPECT_NEAR(inv(1, 0), 0.0, 1e-10);
  EXPECT_NEAR(inv(1, 1), 0.5, 1e-10);
  EXPECT_NEAR(inv(1, 2), 0.0, 1e-10);
  EXPECT_NEAR(inv(2, 0), 0.0, 1e-10);
  EXPECT_NEAR(inv(2, 1), 0.0, 1e-10);
  EXPECT_NEAR(inv(2, 2), 1.0 / 3.0, 1e-10);
}

TEST_F(MatrixTest, SpatialMatrixSolve)
{
  SpatialMatrix<Real> A(2, 2);
  A(0, 0) = 2; A(0, 1) = 1;
  A(1, 0) = 1; A(1, 1) = 3;

  SpatialVector<Real> b(2);
  b(0) = 5;
  b(1) = 7;

  auto x = A.solve(b);
  EXPECT_NEAR(x(0), 1.6, 1e-10);
  EXPECT_NEAR(x(1), 1.8, 1e-10);
}

TEST_F(MatrixTest, SpatialMatrixRowCol)
{
  SpatialMatrix<Real> sm(2, 3);
  sm(0, 0) = 1; sm(0, 1) = 2; sm(0, 2) = 3;
  sm(1, 0) = 4; sm(1, 1) = 5; sm(1, 2) = 6;

  auto r = sm.row(0);
  EXPECT_EQ(r.size(), 3);
  EXPECT_DOUBLE_EQ(r(0), 1.0);
  EXPECT_DOUBLE_EQ(r(1), 2.0);
  EXPECT_DOUBLE_EQ(r(2), 3.0);

  auto c = sm.col(1);
  EXPECT_EQ(c.size(), 2);
  EXPECT_DOUBLE_EQ(c(0), 2.0);
  EXPECT_DOUBLE_EQ(c(1), 5.0);
}

TEST_F(MatrixTest, SpatialMatrixScalarMultiply)
{
  SpatialMatrix<Real> sm(2, 2);
  sm(0, 0) = 1; sm(0, 1) = 2;
  sm(1, 0) = 3; sm(1, 1) = 4;

  sm *= 2.0;
  EXPECT_DOUBLE_EQ(sm(0, 0), 2.0);
  EXPECT_DOUBLE_EQ(sm(0, 1), 4.0);
  EXPECT_DOUBLE_EQ(sm(1, 0), 6.0);
  EXPECT_DOUBLE_EQ(sm(1, 1), 8.0);
}

TEST_F(MatrixTest, SpatialMatrixMatrixMultiply)
{
  SpatialMatrix<Real> a(2, 2);
  a(0, 0) = 1; a(0, 1) = 2;
  a(1, 0) = 3; a(1, 1) = 4;

  SpatialMatrix<Real> b(2, 2);
  b(0, 0) = 5; b(0, 1) = 6;
  b(1, 0) = 7; b(1, 1) = 8;

  a *= b;
  EXPECT_DOUBLE_EQ(a(0, 0), 19.0);
  EXPECT_DOUBLE_EQ(a(0, 1), 22.0);
  EXPECT_DOUBLE_EQ(a(1, 0), 43.0);
  EXPECT_DOUBLE_EQ(a(1, 1), 50.0);
}

TEST_F(MatrixTest, SpatialMatrixPseudoInverse)
{
  SpatialMatrix<Real> sm(2, 2);
  sm.setIdentity();

  auto pi = sm.pseudoInverse();
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
    {
      if (i == j)
        EXPECT_NEAR(pi(i, j), 1.0, 1e-10);
      else
        EXPECT_NEAR(pi(i, j), 0.0, 1e-10);
    }
}

TEST_F(MatrixTest, SpatialMatrixNonSquare)
{
  SpatialMatrix<Real> sm(2, 3);
  sm(0, 0) = 1; sm(0, 1) = 2; sm(0, 2) = 3;
  sm(1, 0) = 4; sm(1, 1) = 5; sm(1, 2) = 6;

  EXPECT_EQ(sm.rows(), 2);
  EXPECT_EQ(sm.cols(), 3);

  auto t = sm.transpose();
  EXPECT_EQ(t.rows(), 3);
  EXPECT_EQ(t.cols(), 2);
}

TEST_F(MatrixTest, SpatialMatrixBinaryAdd)
{
  SpatialMatrix<Real> a(2, 2);
  a(0, 0) = 1; a(0, 1) = 2;
  a(1, 0) = 3; a(1, 1) = 4;

  SpatialMatrix<Real> b(2, 2);
  b(0, 0) = 5; b(0, 1) = 6;
  b(1, 0) = 7; b(1, 1) = 8;

  auto c = a + b;
  EXPECT_DOUBLE_EQ(c(0, 0), 6.0);
  EXPECT_DOUBLE_EQ(c(0, 1), 8.0);
  EXPECT_DOUBLE_EQ(c(1, 0), 10.0);
  EXPECT_DOUBLE_EQ(c(1, 1), 12.0);
}

TEST_F(MatrixTest, SpatialMatrixBinarySubtract)
{
  SpatialMatrix<Real> a(2, 2);
  a(0, 0) = 5; a(0, 1) = 6;
  a(1, 0) = 7; a(1, 1) = 8;

  SpatialMatrix<Real> b(2, 2);
  b(0, 0) = 1; b(0, 1) = 2;
  b(1, 0) = 3; b(1, 1) = 4;

  auto c = a - b;
  EXPECT_DOUBLE_EQ(c(0, 0), 4.0);
  EXPECT_DOUBLE_EQ(c(0, 1), 4.0);
  EXPECT_DOUBLE_EQ(c(1, 0), 4.0);
  EXPECT_DOUBLE_EQ(c(1, 1), 4.0);
}

TEST_F(MatrixTest, SpatialMatrixBinaryMultiply)
{
  SpatialMatrix<Real> a(2, 2);
  a(0, 0) = 1; a(0, 1) = 2;
  a(1, 0) = 3; a(1, 1) = 4;

  SpatialMatrix<Real> b(2, 2);
  b(0, 0) = 5; b(0, 1) = 6;
  b(1, 0) = 7; b(1, 1) = 8;

  auto c = a * b;
  EXPECT_DOUBLE_EQ(c(0, 0), 19.0);
  EXPECT_DOUBLE_EQ(c(0, 1), 22.0);
  EXPECT_DOUBLE_EQ(c(1, 0), 43.0);
  EXPECT_DOUBLE_EQ(c(1, 1), 50.0);
}

TEST_F(MatrixTest, SpatialMatrixBinarySubtract3x3)
{
  SpatialMatrix<Real> a(3, 3);
  a.setIdentity();

  SpatialMatrix<Real> b(3, 3);
  b.setIdentity();

  auto c = a - b;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_DOUBLE_EQ(c(i, j), 0.0);
}

TEST_F(MatrixTest, SpatialMatrixBinarySubtractNonSquare)
{
  SpatialMatrix<Real> a(2, 3);
  a(0, 0) = 1; a(0, 1) = 2; a(0, 2) = 3;
  a(1, 0) = 4; a(1, 1) = 5; a(1, 2) = 6;

  SpatialMatrix<Real> b(2, 3);
  b(0, 0) = 6; b(0, 1) = 5; b(0, 2) = 4;
  b(1, 0) = 3; b(1, 1) = 2; b(1, 2) = 1;

  auto c = a - b;
  EXPECT_DOUBLE_EQ(c(0, 0), -5.0);
  EXPECT_DOUBLE_EQ(c(0, 1), -3.0);
  EXPECT_DOUBLE_EQ(c(0, 2), -1.0);
  EXPECT_DOUBLE_EQ(c(1, 0), 1.0);
  EXPECT_DOUBLE_EQ(c(1, 1), 3.0);
  EXPECT_DOUBLE_EQ(c(1, 2), 5.0);
}
