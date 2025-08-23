/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

class SparseMatrixTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test basic sparse matrix type aliases
TEST_F(SparseMatrixTest, TypeAliases)
{
  // Test SparseMatrix type alias
  static_assert(std::is_same_v<SparseMatrix<Real>, Eigen::SparseMatrix<Real>>);
  static_assert(std::is_same_v<SparseMatrix<Complex>, Eigen::SparseMatrix<Complex>>);
}

// Test sparse matrix construction and basic properties
TEST_F(SparseMatrixTest, Construction)
{
  // Test default construction
  SparseMatrix<Real> sm1;
  EXPECT_EQ(sm1.rows(), 0);
  EXPECT_EQ(sm1.cols(), 0);
  EXPECT_EQ(sm1.nonZeros(), 0);

  // Test sized construction
  SparseMatrix<Real> sm2(5, 5);
  EXPECT_EQ(sm2.rows(), 5);
  EXPECT_EQ(sm2.cols(), 5);
  EXPECT_EQ(sm2.nonZeros(), 0);  // No elements inserted yet

  // Test construction with reserve
  SparseMatrix<Real> sm3(3, 3);
  sm3.reserve(10);  // Reserve space for 10 non-zeros
  EXPECT_EQ(sm3.rows(), 3);
  EXPECT_EQ(sm3.cols(), 3);
}

// Test sparse matrix element insertion and access
TEST_F(SparseMatrixTest, ElementInsertionAndAccess)
{
  SparseMatrix<Real> sm(3, 3);

  // Insert elements using triplets
  std::vector<Eigen::Triplet<Real>> triplets;
  triplets.push_back(Eigen::Triplet<Real>(0, 0, 1.0));
  triplets.push_back(Eigen::Triplet<Real>(0, 2, 2.0));
  triplets.push_back(Eigen::Triplet<Real>(1, 1, 3.0));
  triplets.push_back(Eigen::Triplet<Real>(2, 0, 4.0));
  triplets.push_back(Eigen::Triplet<Real>(2, 2, 5.0));

  sm.setFromTriplets(triplets.begin(), triplets.end());

  EXPECT_EQ(sm.nonZeros(), 5);

  // Test coefficient access
  EXPECT_DOUBLE_EQ(sm.coeff(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(sm.coeff(0, 2), 2.0);
  EXPECT_DOUBLE_EQ(sm.coeff(1, 1), 3.0);
  EXPECT_DOUBLE_EQ(sm.coeff(2, 0), 4.0);
  EXPECT_DOUBLE_EQ(sm.coeff(2, 2), 5.0);

  // Test zero elements
  EXPECT_DOUBLE_EQ(sm.coeff(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(sm.coeff(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(sm.coeff(1, 2), 0.0);
}

// Test sparse matrix arithmetic operations
TEST_F(SparseMatrixTest, ArithmeticOperations)
{
  SparseMatrix<Real> sm1(3, 3);
  SparseMatrix<Real> sm2(3, 3);

  // Initialize first matrix
  std::vector<Eigen::Triplet<Real>> triplets1;
  triplets1.push_back(Eigen::Triplet<Real>(0, 0, 1.0));
  triplets1.push_back(Eigen::Triplet<Real>(1, 1, 2.0));
  triplets1.push_back(Eigen::Triplet<Real>(2, 2, 3.0));
  sm1.setFromTriplets(triplets1.begin(), triplets1.end());

  // Initialize second matrix
  std::vector<Eigen::Triplet<Real>> triplets2;
  triplets2.push_back(Eigen::Triplet<Real>(0, 0, 4.0));
  triplets2.push_back(Eigen::Triplet<Real>(0, 1, 5.0));
  triplets2.push_back(Eigen::Triplet<Real>(1, 1, 6.0));
  sm2.setFromTriplets(triplets2.begin(), triplets2.end());

  // Test matrix addition
  SparseMatrix<Real> sum = sm1 + sm2;
  EXPECT_DOUBLE_EQ(sum.coeff(0, 0), 5.0);  // 1 + 4
  EXPECT_DOUBLE_EQ(sum.coeff(0, 1), 5.0);  // 0 + 5
  EXPECT_DOUBLE_EQ(sum.coeff(1, 1), 8.0);  // 2 + 6
  EXPECT_DOUBLE_EQ(sum.coeff(2, 2), 3.0);  // 3 + 0

  // Test matrix subtraction
  SparseMatrix<Real> diff = sm2 - sm1;
  EXPECT_DOUBLE_EQ(diff.coeff(0, 0), 3.0);  // 4 - 1
  EXPECT_DOUBLE_EQ(diff.coeff(0, 1), 5.0);  // 5 - 0
  EXPECT_DOUBLE_EQ(diff.coeff(1, 1), 4.0);  // 6 - 2
  EXPECT_DOUBLE_EQ(diff.coeff(2, 2), -3.0); // 0 - 3

  // Test scalar multiplication
  SparseMatrix<Real> scaled = 2.0 * sm1;
  EXPECT_DOUBLE_EQ(scaled.coeff(0, 0), 2.0);
  EXPECT_DOUBLE_EQ(scaled.coeff(1, 1), 4.0);
  EXPECT_DOUBLE_EQ(scaled.coeff(2, 2), 6.0);
}

// Test sparse matrix-vector multiplication
TEST_F(SparseMatrixTest, MatrixVectorMultiplication)
{
  SparseMatrix<Real> sm(3, 3);

  // Create a simple sparse matrix
  std::vector<Eigen::Triplet<Real>> triplets;
  triplets.push_back(Eigen::Triplet<Real>(0, 0, 2.0));
  triplets.push_back(Eigen::Triplet<Real>(0, 1, 1.0));
  triplets.push_back(Eigen::Triplet<Real>(1, 1, 3.0));
  triplets.push_back(Eigen::Triplet<Real>(2, 0, 1.0));
  triplets.push_back(Eigen::Triplet<Real>(2, 2, 4.0));
  sm.setFromTriplets(triplets.begin(), triplets.end());

  // Create a vector
  Vector<Real> v(3);
  v << 1, 2, 3;

  // Test matrix-vector multiplication
  Vector<Real> result = sm * v;
  EXPECT_EQ(result.size(), 3);
  EXPECT_DOUBLE_EQ(result[0], 4.0);  // 2*1 + 1*2 + 0*3 = 4
  EXPECT_DOUBLE_EQ(result[1], 6.0);  // 0*1 + 3*2 + 0*3 = 6
  EXPECT_DOUBLE_EQ(result[2], 13.0); // 1*1 + 0*2 + 4*3 = 13
}

// Test sparse matrix-matrix multiplication
TEST_F(SparseMatrixTest, MatrixMatrixMultiplication)
{
  SparseMatrix<Real> sm1(2, 3);
  SparseMatrix<Real> sm2(3, 2);

  // Initialize first matrix (2x3)
  std::vector<Eigen::Triplet<Real>> triplets1;
  triplets1.push_back(Eigen::Triplet<Real>(0, 0, 1.0));
  triplets1.push_back(Eigen::Triplet<Real>(0, 2, 2.0));
  triplets1.push_back(Eigen::Triplet<Real>(1, 1, 3.0));
  sm1.setFromTriplets(triplets1.begin(), triplets1.end());

  // Initialize second matrix (3x2)
  std::vector<Eigen::Triplet<Real>> triplets2;
  triplets2.push_back(Eigen::Triplet<Real>(0, 0, 4.0));
  triplets2.push_back(Eigen::Triplet<Real>(1, 1, 5.0));
  triplets2.push_back(Eigen::Triplet<Real>(2, 0, 6.0));
  sm2.setFromTriplets(triplets2.begin(), triplets2.end());

  // Test matrix-matrix multiplication
  SparseMatrix<Real> result = sm1 * sm2;
  EXPECT_EQ(result.rows(), 2);
  EXPECT_EQ(result.cols(), 2);

  // Expected: [1*4+0*0+2*6, 1*0+0*5+2*0] = [16, 0]
  //           [0*4+3*0+0*6, 0*0+3*5+0*0] = [0, 15]
  EXPECT_DOUBLE_EQ(result.coeff(0, 0), 16.0);
  EXPECT_DOUBLE_EQ(result.coeff(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(result.coeff(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(result.coeff(1, 1), 15.0);
}

// Test axpy function (y = alpha * x + y)
TEST_F(SparseMatrixTest, AxpyOperation)
{
  SparseMatrix<Real> x(3, 3);
  SparseMatrix<Real> y(3, 3);

  // Initialize x matrix
  std::vector<Eigen::Triplet<Real>> tripletsX;
  tripletsX.push_back(Eigen::Triplet<Real>(0, 0, 1.0));
  tripletsX.push_back(Eigen::Triplet<Real>(1, 1, 2.0));
  tripletsX.push_back(Eigen::Triplet<Real>(2, 2, 3.0));
  x.setFromTriplets(tripletsX.begin(), tripletsX.end());

  // Initialize y matrix
  std::vector<Eigen::Triplet<Real>> tripletsY;
  tripletsY.push_back(Eigen::Triplet<Real>(0, 0, 4.0));
  tripletsY.push_back(Eigen::Triplet<Real>(0, 1, 5.0));
  tripletsY.push_back(Eigen::Triplet<Real>(1, 1, 6.0));
  y.setFromTriplets(tripletsY.begin(), tripletsY.end());

  // Test axpy: y = 2.0 * x + y
  Real alpha = 2.0;
  axpy(y, alpha, x);

  // Expected results: y = 2*x + y_original
  EXPECT_DOUBLE_EQ(y.coeff(0, 0), 6.0);  // 2*1 + 4 = 6
  EXPECT_DOUBLE_EQ(y.coeff(0, 1), 5.0);  // 2*0 + 5 = 5
  EXPECT_DOUBLE_EQ(y.coeff(1, 1), 10.0); // 2*2 + 6 = 10
  EXPECT_DOUBLE_EQ(y.coeff(2, 2), 6.0);  // 2*3 + 0 = 6
}

// Test sparse matrix transpose
TEST_F(SparseMatrixTest, Transpose)
{
  SparseMatrix<Real> sm(2, 3);

  // Initialize matrix
  std::vector<Eigen::Triplet<Real>> triplets;
  triplets.push_back(Eigen::Triplet<Real>(0, 0, 1.0));
  triplets.push_back(Eigen::Triplet<Real>(0, 2, 2.0));
  triplets.push_back(Eigen::Triplet<Real>(1, 1, 3.0));
  sm.setFromTriplets(triplets.begin(), triplets.end());

  // Test transpose
  SparseMatrix<Real> transposed = sm.transpose();
  EXPECT_EQ(transposed.rows(), 3);
  EXPECT_EQ(transposed.cols(), 2);

  // Check transposed values
  EXPECT_DOUBLE_EQ(transposed.coeff(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(transposed.coeff(1, 1), 3.0);
  EXPECT_DOUBLE_EQ(transposed.coeff(2, 0), 2.0);
  EXPECT_DOUBLE_EQ(transposed.coeff(0, 1), 0.0);
}

// Test sparse matrix properties and operations
TEST_F(SparseMatrixTest, MatrixProperties)
{
  SparseMatrix<Real> sm(3, 3);

  // Create a matrix with some structure
  std::vector<Eigen::Triplet<Real>> triplets;
  triplets.push_back(Eigen::Triplet<Real>(0, 0, 4.0));
  triplets.push_back(Eigen::Triplet<Real>(0, 1, 1.0));
  triplets.push_back(Eigen::Triplet<Real>(1, 0, 1.0));
  triplets.push_back(Eigen::Triplet<Real>(1, 1, 3.0));
  triplets.push_back(Eigen::Triplet<Real>(2, 2, 2.0));
  sm.setFromTriplets(triplets.begin(), triplets.end());

  // Test number of non-zeros
  EXPECT_EQ(sm.nonZeros(), 5);

  // Test norm (Frobenius norm)
  Real norm = sm.norm();
  EXPECT_GT(norm, 0.0);

  // Test trace
  Real trace = 0.0;
  for (int i = 0; i < sm.rows(); ++i)
    trace += sm.coeff(i, i);
  EXPECT_DOUBLE_EQ(trace, 9.0);  // 4 + 3 + 2 = 9
}

// Test sparse matrix resizing and reserve
TEST_F(SparseMatrixTest, ResizingAndReserve)
{
  SparseMatrix<Real> sm(2, 2);

  // Add some elements
  std::vector<Eigen::Triplet<Real>> triplets;
  triplets.push_back(Eigen::Triplet<Real>(0, 0, 1.0));
  triplets.push_back(Eigen::Triplet<Real>(1, 1, 2.0));
  sm.setFromTriplets(triplets.begin(), triplets.end());

  // Test resize
  sm.conservativeResize(4, 4);
  EXPECT_EQ(sm.rows(), 4);
  EXPECT_EQ(sm.cols(), 4);

  // Previous values should be preserved
  EXPECT_DOUBLE_EQ(sm.coeff(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(sm.coeff(1, 1), 2.0);

  // Test reserve
  sm.reserve(20);
  // Reserve doesn't change the matrix structure, just internal storage
  EXPECT_EQ(sm.rows(), 4);
  EXPECT_EQ(sm.cols(), 4);
}

// Test complex sparse matrix
TEST_F(SparseMatrixTest, ComplexSparseMatrix)
{
  SparseMatrix<Complex> csm(2, 2);

  // Initialize with complex values
  std::vector<Eigen::Triplet<Complex>> triplets;
  triplets.push_back(Eigen::Triplet<Complex>(0, 0, Complex(1.0, 2.0)));
  triplets.push_back(Eigen::Triplet<Complex>(0, 1, Complex(3.0, -1.0)));
  triplets.push_back(Eigen::Triplet<Complex>(1, 0, Complex(-2.0, 1.0)));
  triplets.push_back(Eigen::Triplet<Complex>(1, 1, Complex(0.0, 3.0)));
  csm.setFromTriplets(triplets.begin(), triplets.end());

  // Test conjugate transpose
  SparseMatrix<Complex> hermitian = csm.adjoint();
  EXPECT_EQ(hermitian.rows(), 2);
  EXPECT_EQ(hermitian.cols(), 2);

  // Check conjugate transpose values
  EXPECT_DOUBLE_EQ(hermitian.coeff(0, 0).real(), 1.0);
  EXPECT_DOUBLE_EQ(hermitian.coeff(0, 0).imag(), -2.0);
  EXPECT_DOUBLE_EQ(hermitian.coeff(1, 0).real(), 3.0);
  EXPECT_DOUBLE_EQ(hermitian.coeff(1, 0).imag(), 1.0);
}

// Test sparse matrix iterators
TEST_F(SparseMatrixTest, Iterators)
{
  SparseMatrix<Real> sm(3, 3);

  // Initialize matrix
  std::vector<Eigen::Triplet<Real>> triplets;
  triplets.push_back(Eigen::Triplet<Real>(0, 0, 1.0));
  triplets.push_back(Eigen::Triplet<Real>(1, 1, 2.0));
  triplets.push_back(Eigen::Triplet<Real>(2, 2, 3.0));
  sm.setFromTriplets(triplets.begin(), triplets.end());

  // Test iteration over non-zero elements
  Real sum = 0.0;
  for (int k = 0; k < sm.outerSize(); ++k)
  {
    for (SparseMatrix<Real>::InnerIterator it(sm, k); it; ++it)
    {
      sum += it.value();
    }
  }
  EXPECT_DOUBLE_EQ(sum, 6.0);  // 1 + 2 + 3 = 6
}
