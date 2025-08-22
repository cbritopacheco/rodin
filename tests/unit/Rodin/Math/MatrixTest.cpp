/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/Matrix.h"
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
  
  // Test SpatialMatrix type
  static_assert(std::is_same_v<SpatialMatrix<Real>, 
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, 0, RODIN_MAXIMAL_SPACE_DIMENSION, RODIN_MAXIMAL_SPACE_DIMENSION>>);
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

// Test PointMatrix functionality  
TEST_F(MatrixTest, PointMatrix)
{
  PointMatrix pm(3, 4);  // 3D points, 4 of them
  pm.setZero();
  
  EXPECT_EQ(pm.rows(), 3);
  EXPECT_EQ(pm.cols(), 4);
  
  // Set some point coordinates
  pm(0, 0) = 1.0; // x-coord of first point
  pm(1, 0) = 2.0; // y-coord of first point
  pm(2, 0) = 3.0; // z-coord of first point
  
  EXPECT_DOUBLE_EQ(pm(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(pm(1, 0), 2.0);
  EXPECT_DOUBLE_EQ(pm(2, 0), 3.0);
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