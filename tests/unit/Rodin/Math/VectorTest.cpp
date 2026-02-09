/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

class VectorTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test basic vector type aliases
TEST_F(VectorTest, TypeAliases)
{
  // Test Vector type alias
  static_assert(std::is_same_v<Vector<Real>, Eigen::VectorX<Real>>);
  static_assert(std::is_same_v<Vector<Complex>, ComplexVector>);

  // Test SpatialPoint type
  static_assert(std::is_same_v<SpatialPoint, SpatialVector<Real>>);

  // Test fixed size vectors
  static_assert(std::is_same_v<Vector2<Real>, FixedSizeVector<Real, 2>>);
  static_assert(std::is_same_v<Vector3<Real>, FixedSizeVector<Real, 3>>);
  static_assert(std::is_same_v<Vector4<Real>, FixedSizeVector<Real, 4>>);
  static_assert(std::is_same_v<Vector8<Real>, FixedSizeVector<Real, 8>>);
  static_assert(std::is_same_v<Vector16<Real>, FixedSizeVector<Real, 16>>);
}

// Test vector construction and basic operations
TEST_F(VectorTest, Construction)
{
  // Test default construction
  Vector<Real> v1;
  EXPECT_EQ(v1.size(), 0);
  
  // Test sized construction
  Vector<Real> v2(5);
  EXPECT_EQ(v2.size(), 5);
  
  // Test zero initialization
  Vector<Real> v3 = Vector<Real>::Zero(3);
  EXPECT_EQ(v3.size(), 3);
  for (int i = 0; i < v3.size(); ++i)
    EXPECT_DOUBLE_EQ(v3[i], 0.0);
  
  // Test ones initialization
  Vector<Real> v4 = Vector<Real>::Ones(4);
  EXPECT_EQ(v4.size(), 4);
  for (int i = 0; i < v4.size(); ++i)
    EXPECT_DOUBLE_EQ(v4[i], 1.0);
}

// Test vector element access and assignment
TEST_F(VectorTest, AccessAndAssignment)
{
  Vector<Real> v(3);
  
  // Test element assignment using operator[]
  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;
  
  // Test element access using operator[]
  EXPECT_DOUBLE_EQ(v[0], 1.0);
  EXPECT_DOUBLE_EQ(v[1], 2.0);
  EXPECT_DOUBLE_EQ(v[2], 3.0);
  
  // Test element access using operator()
  EXPECT_DOUBLE_EQ(v(0), 1.0);
  EXPECT_DOUBLE_EQ(v(1), 2.0);
  EXPECT_DOUBLE_EQ(v(2), 3.0);
  
  // Test assignment operator
  v(1) = 5.0;
  EXPECT_DOUBLE_EQ(v(1), 5.0);
}

// Test vector arithmetic operations
TEST_F(VectorTest, ArithmeticOperations)
{
  Vector<Real> v1(3);
  v1 << 1, 2, 3;
  
  Vector<Real> v2(3);
  v2 << 4, 5, 6;
  
  // Test vector addition
  Vector<Real> sum = v1 + v2;
  EXPECT_DOUBLE_EQ(sum[0], 5.0);
  EXPECT_DOUBLE_EQ(sum[1], 7.0);
  EXPECT_DOUBLE_EQ(sum[2], 9.0);
  
  // Test vector subtraction
  Vector<Real> diff = v2 - v1;
  EXPECT_DOUBLE_EQ(diff[0], 3.0);
  EXPECT_DOUBLE_EQ(diff[1], 3.0);
  EXPECT_DOUBLE_EQ(diff[2], 3.0);
  
  // Test scalar multiplication
  Vector<Real> scaled = 2.0 * v1;
  EXPECT_DOUBLE_EQ(scaled[0], 2.0);
  EXPECT_DOUBLE_EQ(scaled[1], 4.0);
  EXPECT_DOUBLE_EQ(scaled[2], 6.0);
  
  // Test scalar division
  Vector<Real> divided = v2 / 2.0;
  EXPECT_DOUBLE_EQ(divided[0], 2.0);
  EXPECT_DOUBLE_EQ(divided[1], 2.5);
  EXPECT_DOUBLE_EQ(divided[2], 3.0);
}

// Test dot product and vector norms
TEST_F(VectorTest, DotProductAndNorms)
{
  Vector<Real> v1(3);
  v1 << 1, 2, 3;
  
  Vector<Real> v2(3);
  v2 << 4, 5, 6;
  
  // Test dot product
  Real dotProduct = v1.dot(v2);
  EXPECT_DOUBLE_EQ(dotProduct, 32.0);  // 1*4 + 2*5 + 3*6 = 32
  
  // Test L2 norm (Euclidean norm)
  Real norm = v1.norm();
  EXPECT_DOUBLE_EQ(norm, std::sqrt(14.0));  // sqrt(1^2 + 2^2 + 3^2) = sqrt(14)
  
  // Test squared norm
  Real squaredNorm = v1.squaredNorm();
  EXPECT_DOUBLE_EQ(squaredNorm, 14.0);  // 1^2 + 2^2 + 3^2 = 14
  
  // Test normalized vector
  Vector<Real> normalized = v1.normalized();
  Real normalizedNorm = normalized.norm();
  EXPECT_NEAR(normalizedNorm, 1.0, 1e-10);
}

// Test FixedSizeVector functionality
TEST_F(VectorTest, FixedSizeVectors)
{
  // Test Vector2
  Vector2<Real> v2;
  v2 << 1.0, 2.0;
  EXPECT_EQ(v2.size(), 2);
  EXPECT_DOUBLE_EQ(v2[0], 1.0);
  EXPECT_DOUBLE_EQ(v2[1], 2.0);
  
  // Test Vector3
  Vector3<Real> v3;
  v3 << 1.0, 2.0, 3.0;
  EXPECT_EQ(v3.size(), 3);
  EXPECT_DOUBLE_EQ(v3[0], 1.0);
  EXPECT_DOUBLE_EQ(v3[1], 2.0);
  EXPECT_DOUBLE_EQ(v3[2], 3.0);
  
  // Test Vector4
  Vector4<Real> v4;
  v4 << 1.0, 2.0, 3.0, 4.0;
  EXPECT_EQ(v4.size(), 4);
  EXPECT_DOUBLE_EQ(v4[3], 4.0);
  
  // Test compile-time size
  static_assert(Vector2<Real>::SizeAtCompileTime == 2);
  static_assert(Vector3<Real>::SizeAtCompileTime == 3);
  static_assert(Vector4<Real>::SizeAtCompileTime == 4);
  static_assert(Vector8<Real>::SizeAtCompileTime == 8);
  static_assert(Vector16<Real>::SizeAtCompileTime == 16);
}

// Test SpatialVector functionality
TEST_F(VectorTest, SpatialVector)
{
  SpatialVector<Real> sv(3);
  sv[0] = 1.0;
  sv[1] = 2.0;
  sv[2] = 3.0;

  EXPECT_EQ(sv.size(), 3);
  EXPECT_DOUBLE_EQ(sv[0], 1.0);
  EXPECT_DOUBLE_EQ(sv[1], 2.0);
  EXPECT_DOUBLE_EQ(sv[2], 3.0);

  // Test SpatialPoint (alias for SpatialVector<Real>)
  SpatialPoint sp(3);
  sp[0] = 4.0;
  sp[1] = 5.0;
  sp[2] = 6.0;

  EXPECT_EQ(sp.size(), 3);
  EXPECT_DOUBLE_EQ(sp[0], 4.0);
  EXPECT_DOUBLE_EQ(sp[1], 5.0);
  EXPECT_DOUBLE_EQ(sp[2], 6.0);
}

// Test complex vector operations
TEST_F(VectorTest, ComplexVector)
{
  ComplexVector cv(3);
  cv[0] = Complex(1.0, 2.0);  // 1 + 2i
  cv[1] = Complex(3.0, -1.0); // 3 - i
  cv[2] = Complex(-2.0, 1.0); // -2 + i

  // Test conjugate
  ComplexVector conjugated = cv.conjugate();
  EXPECT_DOUBLE_EQ(conjugated[0].real(), 1.0);
  EXPECT_DOUBLE_EQ(conjugated[0].imag(), -2.0);
  EXPECT_DOUBLE_EQ(conjugated[1].real(), 3.0);
  EXPECT_DOUBLE_EQ(conjugated[1].imag(), 1.0);
  EXPECT_DOUBLE_EQ(conjugated[2].real(), -2.0);
  EXPECT_DOUBLE_EQ(conjugated[2].imag(), -1.0);
  
  // Test complex norm
  Real complexNorm = cv.norm();
  EXPECT_GT(complexNorm, 0.0);
}

// Test vector resizing
TEST_F(VectorTest, Resizing)
{
  Vector<Real> v(3);
  v << 1, 2, 3;
  
  // Test resize
  v.resize(5);
  EXPECT_EQ(v.size(), 5);
  
  // Test conservativeResize
  Vector<Real> v2(3);
  v2 << 1, 2, 3;
  v2.conservativeResize(5);
  EXPECT_EQ(v2.size(), 5);
  EXPECT_DOUBLE_EQ(v2[0], 1.0);
  EXPECT_DOUBLE_EQ(v2[1], 2.0);
  EXPECT_DOUBLE_EQ(v2[2], 3.0);
}

// Test vector cross product (for 3D vectors)
TEST_F(VectorTest, CrossProduct)
{
  Vector3<Real> v1;
  v1 << 1, 0, 0;  // x-axis unit vector
  
  Vector3<Real> v2;
  v2 << 0, 1, 0;  // y-axis unit vector
  
  // Cross product should give z-axis unit vector
  Vector3<Real> cross = v1.cross(v2);
  EXPECT_DOUBLE_EQ(cross[0], 0.0);
  EXPECT_DOUBLE_EQ(cross[1], 0.0);
  EXPECT_DOUBLE_EQ(cross[2], 1.0);
  
  // Test anti-commutativity: v1 x v2 = -(v2 x v1)
  Vector3<Real> cross_rev = v2.cross(v1);
  EXPECT_DOUBLE_EQ(cross_rev[0], 0.0);
  EXPECT_DOUBLE_EQ(cross_rev[1], 0.0);
  EXPECT_DOUBLE_EQ(cross_rev[2], -1.0);
}

// Test vector max and min operations
TEST_F(VectorTest, MaxMinOperations)
{
  Vector<Real> v(4);
  v << 3.0, -1.0, 5.0, 2.0;
  
  // Test max coefficient
  Real maxVal = v.maxCoeff();
  EXPECT_DOUBLE_EQ(maxVal, 5.0);
  
  // Test min coefficient
  Real minVal = v.minCoeff();
  EXPECT_DOUBLE_EQ(minVal, -1.0);
  
  // Test max coefficient with index
  int maxIndex;
  Real maxValWithIndex = v.maxCoeff(&maxIndex);
  EXPECT_DOUBLE_EQ(maxValWithIndex, 5.0);
  EXPECT_EQ(maxIndex, 2);
  
  // Test min coefficient with index
  int minIndex;
  Real minValWithIndex = v.minCoeff(&minIndex);
  EXPECT_DOUBLE_EQ(minValWithIndex, -1.0);
  EXPECT_EQ(minIndex, 1);
}

// Test vector segment operations
TEST_F(VectorTest, SegmentOperations)
{
  Vector<Real> v(6);
  v << 1, 2, 3, 4, 5, 6;
  
  // Test head (first n elements)
  auto head3 = v.head(3);
  EXPECT_EQ(head3.size(), 3);
  EXPECT_DOUBLE_EQ(head3[0], 1.0);
  EXPECT_DOUBLE_EQ(head3[1], 2.0);
  EXPECT_DOUBLE_EQ(head3[2], 3.0);
  
  // Test tail (last n elements)
  auto tail3 = v.tail(3);
  EXPECT_EQ(tail3.size(), 3);
  EXPECT_DOUBLE_EQ(tail3[0], 4.0);
  EXPECT_DOUBLE_EQ(tail3[1], 5.0);
  EXPECT_DOUBLE_EQ(tail3[2], 6.0);
  
  // Test segment (middle elements)
  auto segment2 = v.segment(2, 2);  // start at index 2, length 2
  EXPECT_EQ(segment2.size(), 2);
  EXPECT_DOUBLE_EQ(segment2[0], 3.0);
  EXPECT_DOUBLE_EQ(segment2[1], 4.0);
}
