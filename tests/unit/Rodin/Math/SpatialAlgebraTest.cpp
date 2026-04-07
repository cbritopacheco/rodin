#include <gtest/gtest.h>

#include <complex>

#include <Eigen/Dense>

#include "Rodin/Math/Common.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

class SpatialAlgebraTest : public ::testing::Test
{};

TEST_F(SpatialAlgebraTest, SpatialVectorOperations)
{
  SpatialVector<Real> v{1.0, 2.0, 3.0};
  SpatialVector<Real> w{2.0, -1.0, 4.0};

  const auto s = v + w;
  EXPECT_DOUBLE_EQ(s(0), 3.0);
  EXPECT_DOUBLE_EQ(s(1), 1.0);
  EXPECT_DOUBLE_EQ(s(2), 7.0);

  const auto d = v - w;
  EXPECT_DOUBLE_EQ(d(0), -1.0);
  EXPECT_DOUBLE_EQ(d(1), 3.0);
  EXPECT_DOUBLE_EQ(d(2), -1.0);

  const auto c = v.cross(w);
  EXPECT_DOUBLE_EQ(c(0), 11.0);
  EXPECT_DOUBLE_EQ(c(1), 2.0);
  EXPECT_DOUBLE_EQ(c(2), -5.0);

  EXPECT_DOUBLE_EQ(v.dot(w), 12.0);
  EXPECT_DOUBLE_EQ(v.squaredNorm(), 14.0);
  EXPECT_NEAR(v.norm(), std::sqrt(14.0), 1e-12);
  EXPECT_NEAR(v.lpNorm<2>(), v.norm(), 1e-12);
  EXPECT_NEAR(v.stableNorm(), v.norm(), 1e-12);
  EXPECT_NEAR(v.blueNorm(), v.norm(), 1e-12);

  const auto vn = v.normalized();
  EXPECT_NEAR(vn.norm(), 1.0, 1e-12);

  auto vc = v;
  vc.normalize();
  EXPECT_NEAR(vc.norm(), 1.0, 1e-12);

  const auto t = v.transpose();
  ASSERT_EQ(t.rows(), 1);
  ASSERT_EQ(t.cols(), 3);
  EXPECT_DOUBLE_EQ(t(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(t(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(t(0, 2), 3.0);
}

TEST_F(SpatialAlgebraTest, SpatialVectorComplexConjugate)
{
  SpatialVector<Complex> v{
    Complex(1.0, 2.0),
    Complex(-3.0, 1.0),
    Complex(0.5, -0.75)
  };

  const auto c = v.conjugate();
  EXPECT_EQ(c(0), Complex(1.0, -2.0));
  EXPECT_EQ(c(1), Complex(-3.0, -1.0));
  EXPECT_EQ(c(2), Complex(0.5, 0.75));
}

TEST_F(SpatialAlgebraTest, SpatialVectorDotWithEigen)
{
  SpatialVector<Real> v{1.0, 2.0, 3.0};
  Eigen::Vector3d e;
  e << 2.0, 3.0, 4.0;

  EXPECT_DOUBLE_EQ(v.dot(e), 20.0);
  EXPECT_DOUBLE_EQ(Math::dot(v, e), 20.0);
  EXPECT_DOUBLE_EQ(Math::dot(e, v), 20.0);
}

namespace
{
  SpatialMatrix<Real> makeSpatial(std::uint8_t rows, std::uint8_t cols)
  {
    SpatialMatrix<Real> m(rows, cols);
    Real v = 1.0;
    for (std::uint8_t i = 0; i < rows; ++i)
    {
      for (std::uint8_t j = 0; j < cols; ++j)
      {
        m(i, j) = v;
        v += 1.0;
      }
    }
    return m;
  }
}

TEST_F(SpatialAlgebraTest, SpatialMatrixMultiplyAllShapeCombos)
{
  for (std::uint8_t r = 0; r <= 3; ++r)
  {
    for (std::uint8_t k = 0; k <= 3; ++k)
    {
      for (std::uint8_t c = 0; c <= 3; ++c)
      {
        SpatialMatrix<Real> A = makeSpatial(r, k);
        SpatialMatrix<Real> B = makeSpatial(k, c);
        SpatialMatrix<Real> C = A;
        C *= B;

        Eigen::MatrixXd EA = Eigen::MatrixXd::Zero(r, k);
        Eigen::MatrixXd EB = Eigen::MatrixXd::Zero(k, c);
        for (std::uint8_t i = 0; i < r; ++i)
          for (std::uint8_t j = 0; j < k; ++j)
            EA(i, j) = A(i, j);
        for (std::uint8_t i = 0; i < k; ++i)
          for (std::uint8_t j = 0; j < c; ++j)
            EB(i, j) = B(i, j);
        const Eigen::MatrixXd EC = EA * EB;

        EXPECT_EQ(C.rows(), r);
        EXPECT_EQ(C.cols(), c);
        for (std::uint8_t i = 0; i < r; ++i)
        {
          for (std::uint8_t j = 0; j < c; ++j)
          {
            EXPECT_NEAR(C(i, j), EC(i, j), 1e-12);
          }
        }
      }
    }
  }
}

TEST_F(SpatialAlgebraTest, SpatialMatrixSolveDeterminantAndInverse)
{
  SpatialMatrix<Real> A(3, 3);
  A(0, 0) = 4.0; A(0, 1) = 1.0; A(0, 2) = 2.0;
  A(1, 0) = 0.0; A(1, 1) = 3.0; A(1, 2) = -1.0;
  A(2, 0) = 2.0; A(2, 1) = 0.0; A(2, 2) = 5.0;

  const auto det = A.determinant();
  EXPECT_NEAR(det, 50.0, 1e-12);

  const auto inv = A.inverse();
  SpatialMatrix<Real> I = A;
  I *= inv;
  for (std::uint8_t i = 0; i < 3; ++i)
  {
    for (std::uint8_t j = 0; j < 3; ++j)
    {
      EXPECT_NEAR(I(i, j), i == j ? 1.0 : 0.0, 1e-10);
    }
  }

  SpatialVector<Real> b{1.0, 2.0, 3.0};
  const auto x = A.solve(b);
  const auto check = A * x;
  EXPECT_NEAR(check(0, 0), b(0), 1e-10);
  EXPECT_NEAR(check(1, 0), b(1), 1e-10);
  EXPECT_NEAR(check(2, 0), b(2), 1e-10);
}

TEST_F(SpatialAlgebraTest, SpatialMatrixPseudoInverseRectangular)
{
  SpatialMatrix<Real> A(3, 2);
  A(0, 0) = 1.0; A(0, 1) = 0.0;
  A(1, 0) = 0.0; A(1, 1) = 1.0;
  A(2, 0) = 1.0; A(2, 1) = 1.0;

  const auto Ap = A.pseudoInverse();
  const auto triple = A * Ap * A;
  for (std::uint8_t i = 0; i < A.rows(); ++i)
  {
    for (std::uint8_t j = 0; j < A.cols(); ++j)
    {
      EXPECT_NEAR(triple(i, j), A(i, j), 1e-8);
    }
  }
}

TEST_F(SpatialAlgebraTest, SpatialMatrixTransposeConjugateAdjointRowColTrace)
{
  SpatialMatrix<Complex> A(2, 2);
  A(0, 0) = Complex(1.0, 2.0);
  A(0, 1) = Complex(2.0, -1.0);
  A(1, 0) = Complex(-3.0, 4.0);
  A(1, 1) = Complex(0.5, 0.25);

  const auto T = A.transpose();
  EXPECT_EQ(T(0, 1), A(1, 0));
  EXPECT_EQ(T(1, 0), A(0, 1));

  const auto C = A.conjugate();
  EXPECT_EQ(C(0, 0), std::conj(A(0, 0)));
  EXPECT_EQ(C(1, 0), std::conj(A(1, 0)));

  const auto H = A.adjoint();
  EXPECT_EQ(H(0, 1), std::conj(A(1, 0)));
  EXPECT_EQ(H(1, 0), std::conj(A(0, 1)));

  SpatialMatrix<Real> R(2, 2);
  R(0, 0) = 3.0; R(0, 1) = 1.0;
  R(1, 0) = 2.0; R(1, 1) = 4.0;
  EXPECT_DOUBLE_EQ(R.trace(), 7.0);

  const auto r0 = R.row(0);
  const auto c1 = R.col(1);
  EXPECT_DOUBLE_EQ(r0(0), 3.0);
  EXPECT_DOUBLE_EQ(r0(1), 1.0);
  EXPECT_DOUBLE_EQ(c1(0), 1.0);
  EXPECT_DOUBLE_EQ(c1(1), 4.0);
}

TEST_F(SpatialAlgebraTest, SpatialMatrixDotOverloads)
{
  SpatialMatrix<Real> A(2, 2);
  A(0, 0) = 1.0; A(0, 1) = 2.0;
  A(1, 0) = 3.0; A(1, 1) = 4.0;

  SpatialMatrix<Real> B(2, 2);
  B(0, 0) = 2.0; B(0, 1) = 0.0;
  B(1, 0) = 1.0; B(1, 1) = 2.0;

  EXPECT_DOUBLE_EQ(Math::dot(A, B), 13.0);

  Eigen::Matrix2d E;
  E << 2.0, 0.0,
       1.0, 2.0;
  EXPECT_DOUBLE_EQ(Math::dot(A, E), 13.0);
  EXPECT_DOUBLE_EQ(Math::dot(E, A), 13.0);
}
