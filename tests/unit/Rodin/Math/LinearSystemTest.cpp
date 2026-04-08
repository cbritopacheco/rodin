/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

class LinearSystemTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test LinearSystem with SparseMatrix
TEST_F(LinearSystemTest, SparseMatrixLinearSystem)
{
  using LSType = LinearSystem<SparseMatrix<Real>, Vector<Real>>;

  // Test default construction
  LSType ls;

  // Test that we can get references to components
  auto& A = ls.getOperator();
  auto& x = ls.getSolution();
  auto& b = ls.getVector();

  // Test that components are initially empty/zero-sized
  EXPECT_EQ(A.rows(), 0);
  EXPECT_EQ(A.cols(), 0);
  EXPECT_EQ(x.size(), 0);
  EXPECT_EQ(b.size(), 0);
}

// Test LinearSystem with dense Matrix
TEST_F(LinearSystemTest, DenseMatrixLinearSystem)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;

  // Test default construction
  LSType ls;

  // Test that we can get references to components
  auto& A = ls.getOperator();
  auto& x = ls.getSolution();
  auto& b = ls.getVector();

  // Test that components are initially empty/zero-sized
  EXPECT_EQ(A.rows(), 0);
  EXPECT_EQ(A.cols(), 0);
  EXPECT_EQ(x.size(), 0);
  EXPECT_EQ(b.size(), 0);
}

// Test LinearSystem copy semantics
TEST_F(LinearSystemTest, CopySemantics)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;

  LSType ls1;

  // Resize components for testing
  ls1.getOperator().resize(3, 3);
  ls1.getSolution().resize(3);
  ls1.getVector().resize(3);

  // Set some values
  ls1.getOperator().setIdentity();
  ls1.getSolution().setOnes();
  ls1.getVector().setOnes();

  // Test copy constructor
  LSType ls2(ls1);
  EXPECT_EQ(ls2.getOperator().rows(), 3);
  EXPECT_EQ(ls2.getOperator().cols(), 3);
  EXPECT_EQ(ls2.getSolution().size(), 3);
  EXPECT_EQ(ls2.getVector().size(), 3);

  // Test copy assignment
  LSType ls3;
  ls3 = ls1;
  EXPECT_EQ(ls3.getOperator().rows(), 3);
  EXPECT_EQ(ls3.getOperator().cols(), 3);
  EXPECT_EQ(ls3.getSolution().size(), 3);
  EXPECT_EQ(ls3.getVector().size(), 3);
}

// Test LinearSystem move semantics
TEST_F(LinearSystemTest, MoveSemantics)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;

  LSType ls1;

  // Resize components for testing
  ls1.getOperator().resize(3, 3);
  ls1.getSolution().resize(3);
  ls1.getVector().resize(3);

  // Test move constructor
  LSType ls2(std::move(ls1));
  EXPECT_EQ(ls2.getOperator().rows(), 3);
  EXPECT_EQ(ls2.getOperator().cols(), 3);
  EXPECT_EQ(ls2.getSolution().size(), 3);
  EXPECT_EQ(ls2.getVector().size(), 3);

  // Test move assignment
  LSType ls3;
  LSType ls4;
  ls4.getOperator().resize(2, 2);
  ls4.getSolution().resize(2);
  ls4.getVector().resize(2);

  ls3 = std::move(ls4);
  EXPECT_EQ(ls3.getOperator().rows(), 2);
  EXPECT_EQ(ls3.getOperator().cols(), 2);
  EXPECT_EQ(ls3.getSolution().size(), 2);
  EXPECT_EQ(ls3.getVector().size(), 2);
}

// Test const access to LinearSystem components
TEST_F(LinearSystemTest, ConstAccess)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;

  LSType ls;
  ls.getOperator().resize(2, 2);
  ls.getSolution().resize(2);
  ls.getVector().resize(2);

  // Test const access
  const LSType& const_ls = ls;
  const auto& A = const_ls.getOperator();
  const auto& x = const_ls.getSolution();
  const auto& b = const_ls.getVector();

  EXPECT_EQ(A.rows(), 2);
  EXPECT_EQ(A.cols(), 2);
  EXPECT_EQ(x.size(), 2);
  EXPECT_EQ(b.size(), 2);
}

// Test LinearSystem eliminate function for dense matrices
TEST_F(LinearSystemTest, DenseMatrixEliminate)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;

  LSType ls;

  // Set up a simple 3x3 system
  ls.getOperator().resize(3, 3);
  ls.getVector().resize(3);

  // Set up matrix and vector
  ls.getOperator() << 2, 1, 0,
                       1, 3, 1,
                       0, 1, 2;

  ls.getVector() << 1, 2, 3;

  // Eliminate DOF 1 (middle row/column) with value 5.0
  IndexMap<Real> dofs;
  dofs[1] = 5.0;

  ls.eliminate(dofs);

  // Check that row 1 and column 1 have been zeroed except for diagonal
  EXPECT_DOUBLE_EQ(ls.getOperator()(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(2, 1), 0.0);

  // Check that diagonal element is 1
  EXPECT_DOUBLE_EQ(ls.getOperator()(1, 1), 1.0);

  // Check that RHS has been set correctly
  EXPECT_DOUBLE_EQ(ls.getVector()(1), 5.0);
}

// Test LinearSystem with complex numbers
TEST_F(LinearSystemTest, ComplexLinearSystem)
{
  using LSType = LinearSystem<Matrix<Complex>, Vector<Complex>>;

  LSType ls;

  // Test that we can get references to components
  auto& A = ls.getOperator();
  auto& x = ls.getSolution();
  auto& b = ls.getVector();

  // Resize and set some complex values
  A.resize(2, 2);
  x.resize(2);
  b.resize(2);

  Complex i(0, 1);  // imaginary unit
  A(0, 0) = Complex(1, 0);
  A(0, 1) = i;
  A(1, 0) = -i;
  A(1, 1) = Complex(1, 0);

  EXPECT_EQ(A(0, 1), i);
  EXPECT_EQ(A(1, 0), -i);
}

// Test tuple-like interface
TEST_F(LinearSystemTest, TupleInterface)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;

  LSType ls;
  ls.getOperator().resize(2, 2);
  ls.getSolution().resize(2);
  ls.getVector().resize(2);

  // Test std::tuple_size
  static_assert(std::tuple_size_v<LSType> == 3);

  // Test std::tuple_element
  static_assert(std::is_same_v<std::tuple_element_t<0, LSType>, Matrix<Real>>);
  static_assert(std::is_same_v<std::tuple_element_t<1, LSType>, Vector<Real>>);
  static_assert(std::is_same_v<std::tuple_element_t<2, LSType>, Vector<Real>>);

  auto& A = Math::get<0>(ls);
  auto& x = Math::get<1>(ls);
  auto& b = Math::get<2>(ls);

  EXPECT_EQ(A.rows(), 2);
  EXPECT_EQ(A.cols(), 2);
  EXPECT_EQ(x.size(), 2);
  EXPECT_EQ(b.size(), 2);
}

// Test Math::get interface
TEST_F(LinearSystemTest, MathGetInterface)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;

  LSType ls;
  ls.getOperator().resize(2, 2);
  ls.getSolution().resize(2);
  ls.getVector().resize(2);

  // Test Math::get<0> (operator)
  auto& A = Math::get<0>(ls);
  EXPECT_EQ(A.rows(), 2);
  EXPECT_EQ(A.cols(), 2);

  // Test Math::get<1> (solution)
  auto& x = Math::get<1>(ls);
  EXPECT_EQ(x.size(), 2);

  // Test Math::get<2> (vector)
  auto& b = Math::get<2>(ls);
  EXPECT_EQ(b.size(), 2);

  // Test const version
  const LSType& const_ls = ls;
  const auto& const_A = Math::get<0>(const_ls);
  const auto& const_x = Math::get<1>(const_ls);
  const auto& const_b = Math::get<2>(const_ls);

  EXPECT_EQ(const_A.rows(), 2);
  EXPECT_EQ(const_x.size(), 2);
  EXPECT_EQ(const_b.size(), 2);
}

// Test FormLanguage::Traits for LinearSystem
TEST_F(LinearSystemTest, FormLanguageTraits)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;
  using TraitsType = FormLanguage::Traits<LSType>;

  // Test that traits are properly defined
  static_assert(std::is_same_v<typename TraitsType::OperatorType, Matrix<Real>>);
  static_assert(std::is_same_v<typename TraitsType::VectorType, Vector<Real>>);
  static_assert(std::is_same_v<typename TraitsType::ScalarType, Real>);
}

// Test LinearSystem type properties
TEST_F(LinearSystemTest, TypeProperties)
{
  using DenseLSType = LinearSystem<Matrix<Real>, Vector<Real>>;
  using SparseLSType = LinearSystem<SparseMatrix<Real>, Vector<Real>>;

  // Test that LinearSystem types are properly constructible
  static_assert(std::is_default_constructible_v<DenseLSType>);
  static_assert(std::is_copy_constructible_v<DenseLSType>);
  static_assert(std::is_move_constructible_v<DenseLSType>);
  static_assert(std::is_copy_assignable_v<DenseLSType>);
  static_assert(std::is_move_assignable_v<DenseLSType>);

  static_assert(std::is_default_constructible_v<SparseLSType>);
  static_assert(std::is_copy_constructible_v<SparseLSType>);
  static_assert(std::is_move_constructible_v<SparseLSType>);
  static_assert(std::is_copy_assignable_v<SparseLSType>);
  static_assert(std::is_move_assignable_v<SparseLSType>);
}

// Test sparse matrix eliminate
TEST_F(LinearSystemTest, SparseMatrixEliminate)
{
  using LSType = LinearSystem<SparseMatrix<Real>, Vector<Real>>;

  LSType ls;

  // Build a simple 3x3 sparse system
  SparseMatrix<Real> A(3, 3);
  std::vector<Eigen::Triplet<Real>> trips;
  trips.emplace_back(0, 0, 2.0);
  trips.emplace_back(0, 1, 1.0);
  trips.emplace_back(1, 0, 1.0);
  trips.emplace_back(1, 1, 3.0);
  trips.emplace_back(1, 2, 1.0);
  trips.emplace_back(2, 1, 1.0);
  trips.emplace_back(2, 2, 2.0);
  A.setFromTriplets(trips.begin(), trips.end());

  ls.getOperator() = A;
  ls.getVector().resize(3);
  ls.getVector() << 1, 2, 3;
  ls.getSolution().resize(3);

  // Eliminate DOF 1 with value 5.0
  IndexMap<Real> dofs;
  dofs[1] = 5.0;
  ls.eliminate(dofs);

  // Check diagonal is 1
  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(1, 1), 1.0);

  // Check off-diagonals in row/col 1 are zeroed
  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(2, 1), 0.0);

  // Check RHS at eliminated DOF
  EXPECT_DOUBLE_EQ(ls.getVector()(1), 5.0);
}

// Test multiple DOF elimination
TEST_F(LinearSystemTest, MultipleDOFElimination)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;

  LSType ls;
  ls.getOperator().resize(4, 4);
  ls.getVector().resize(4);
  ls.getSolution().resize(4);

  ls.getOperator() << 4, 1, 0, 0,
                       1, 4, 1, 0,
                       0, 1, 4, 1,
                       0, 0, 1, 4;
  ls.getVector() << 1, 2, 3, 4;

  // Eliminate DOFs 0 and 3 simultaneously
  IndexMap<Real> dofs;
  dofs[0] = 1.0;
  dofs[3] = 2.0;
  ls.eliminate(dofs);

  // Check diagonals of eliminated DOFs
  EXPECT_DOUBLE_EQ(ls.getOperator()(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(3, 3), 1.0);

  // Check RHS at eliminated DOFs
  EXPECT_DOUBLE_EQ(ls.getVector()(0), 1.0);
  EXPECT_DOUBLE_EQ(ls.getVector()(3), 2.0);

  // Check rows/cols of eliminated DOFs are zeroed (except diagonal)
  EXPECT_DOUBLE_EQ(ls.getOperator()(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(3, 2), 0.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(2, 3), 0.0);

  // Interior DOFs should still have their original stiffness entries
  EXPECT_DOUBLE_EQ(ls.getOperator()(1, 1), 4.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(2, 2), 4.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(1, 2), 1.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(2, 1), 1.0);
}

// Test sparse eliminate with offset
TEST_F(LinearSystemTest, SparseEliminateWithOffset)
{
  using LSType = LinearSystem<SparseMatrix<Real>, Vector<Real>>;

  LSType ls;

  // Build a 4x4 sparse identity
  SparseMatrix<Real> A(4, 4);
  std::vector<Eigen::Triplet<Real>> trips;
  for (int i = 0; i < 4; ++i)
    trips.emplace_back(i, i, 2.0);
  trips.emplace_back(0, 1, 1.0);
  trips.emplace_back(1, 0, 1.0);
  trips.emplace_back(2, 3, 1.0);
  trips.emplace_back(3, 2, 1.0);
  A.setFromTriplets(trips.begin(), trips.end());

  ls.getOperator() = A;
  ls.getVector().resize(4);
  ls.getVector() << 1, 2, 3, 4;
  ls.getSolution().resize(4);

  // Eliminate DOF 0 with offset 2 → affects global DOF 2
  IndexMap<Real> dofs;
  dofs[0] = 10.0;
  ls.eliminate(dofs, 2);

  // DOF at index 2 (0+offset) should be eliminated
  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(2, 2), 1.0);
  EXPECT_DOUBLE_EQ(ls.getVector()(2), 10.0);
}

// Test dense end-to-end: eliminate then solve
TEST_F(LinearSystemTest, DenseEliminateThenSolve)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;

  // Set up system: simple Laplacian-like [2,-1,0; -1,2,-1; 0,-1,2] x = [1,0,1]
  // With BCs x(0) = 0, x(2) = 0
  LSType ls;
  ls.getOperator().resize(3, 3);
  ls.getVector().resize(3);
  ls.getSolution().resize(3);

  ls.getOperator() << 2, -1, 0,
                       -1, 2, -1,
                        0, -1, 2;
  ls.getVector() << 0, 1, 0;

  // Eliminate DOFs 0 and 2 with value 0
  IndexMap<Real> dofs;
  dofs[0] = 0.0;
  dofs[2] = 0.0;
  ls.eliminate(dofs);

  // After elimination:
  // Row 0: [1, 0, 0] x = [0]
  // Row 1: [0, 2, 0] x = [1]  (contributions from eliminated DOFs are 0)
  // Row 2: [0, 0, 1] x = [0]
  // So x = [0, 0.5, 0]
  EXPECT_DOUBLE_EQ(ls.getOperator()(1, 1), 2.0);
  EXPECT_NEAR(ls.getVector()(1), 1.0, 1e-10);
}
