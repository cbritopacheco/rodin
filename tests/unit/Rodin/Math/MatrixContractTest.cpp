/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Dense matrix contractions, ownership and Eigen interoperability.
 */
#include <gtest/gtest.h>
#include <Eigen/IterativeLinearSolvers>
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  /// @brief Exercises dynamic, rectangular, square and vector-shaped matrix owners.
  template <class Matrix>
  class MatrixContractTest : public ::testing::Test
  {};

  using MatrixTypes = ::testing::Types<Math::Matrix<Real>, Math::Matrix<Complex>,
    Math::Matrix<std::complex<float>>, Math::Matrix<std::complex<long double>>,
    Math::FixedSizeMatrix<Real, 2, 3>, Math::FixedSizeMatrix<Complex, 2, 3>,
    Math::FixedSizeMatrix<std::complex<float>, 3, 3>,
    Math::FixedSizeMatrix<std::complex<long double>, 2, 3>,
    Math::FixedSizeMatrix<Complex, 1, 3>, Math::FixedSizeMatrix<Complex, 3, 1>,
    Math::FixedSizeMatrix<Complex, 1, 1>>;
  TYPED_TEST_SUITE(MatrixContractTest, MatrixTypes);

  /// @brief Member and free contractions agree with an independent right-conjugating sum.
  TYPED_TEST(MatrixContractTest, DotConjugatesRight)
  {
    using Matrix = TypeParam;
    using Scalar = typename Matrix::Scalar;
    constexpr int Rows =
      Matrix::RowsAtCompileTime == Eigen::Dynamic ? 2 : Matrix::RowsAtCompileTime;
    constexpr int Cols =
      Matrix::ColsAtCompileTime == Eigen::Dynamic ? 3 : Matrix::ColsAtCompileTime;
    Matrix left(Rows, Cols), right(Rows, Cols);
    Scalar expected = 0;
    for (int i = 0; i < Rows; ++i)
    {
      for (int j = 0; j < Cols; ++j)
      {
        if constexpr (Eigen::NumTraits<Scalar>::IsComplex)
        {
          left(i, j) = Scalar(i + j + 1, i + 2);
          right(i, j) = Scalar(j + 2, -i - j - 3);
          expected += left(i, j) * std::conj(right(i, j));
        }
        else
        {
          left(i, j) = Scalar(i + j + 1);
          right(i, j) = Scalar(j + 2);
          expected += left(i, j) * right(i, j);
        }
      }
    }
    EXPECT_EQ(left.dot(right), expected);
    EXPECT_EQ(left.dot(right.block(0, 0, Rows, Cols)), expected);
    EXPECT_EQ(left.dot(Scalar(2) * right), Scalar(2) * expected);
    EXPECT_EQ(Math::dot(left, right), expected);
    EXPECT_EQ(Math::dot(left + left, right), Scalar(2) * expected);
    if constexpr (Eigen::NumTraits<Scalar>::IsComplex)
    {
      EXPECT_NE(expected, std::conj(expected));
      EXPECT_EQ(right.dot(left), std::conj(expected));
      if constexpr (Matrix::IsVectorAtCompileTime)
      {
        const Eigen::Matrix<Scalar, Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime>
          native = left;
        EXPECT_EQ(native.dot(right), std::conj(expected));
      }
    }
  }

  /// @brief Expressions, maps, references and storage transfers retain Eigen's behavior.
  TYPED_TEST(MatrixContractTest, EigenExpressionsAndOwnership)
  {
    using Matrix = TypeParam;
    using Scalar = typename Matrix::Scalar;
    using Parent =
      Eigen::Matrix<Scalar, Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime>;
    constexpr int Rows =
      Matrix::RowsAtCompileTime == Eigen::Dynamic ? 2 : Matrix::RowsAtCompileTime;
    constexpr int Cols =
      Matrix::ColsAtCompileTime == Eigen::Dynamic ? 3 : Matrix::ColsAtCompileTime;
    static_assert(std::is_same_v<typename Matrix::PlainObject, Parent>);
    static_assert(sizeof(Matrix) == sizeof(Parent));
    static_assert(alignof(Matrix) == alignof(Parent));
    static_assert(FormLanguage::IsEigenObject<const Matrix&>::Value);
    static_assert(FormLanguage::IsEigenObject<Parent>::Value);
    static_assert(std::is_same_v<typename FormLanguage::Traits<Matrix>::ScalarType,
      typename FormLanguage::Traits<Parent>::ScalarType>);
    static_assert(
      std::is_same_v<decltype(std::declval<const Matrix&>().eval()), const Matrix&>);
    Matrix left = Matrix::Ones(Rows, Cols);
    EXPECT_EQ(left.eval().data(), left.data());
    const Matrix right = Matrix::Constant(Rows, Cols, Scalar(2));
    Matrix sum = left + right;
    EXPECT_EQ(sum(0, 0), Scalar(3));
    sum = Scalar(2) * left + right;
    EXPECT_EQ(sum(0, 0), Scalar(4));
    Matrix copy(sum);
    copy(0, 0) = Scalar(7);
    EXPECT_EQ(sum(0, 0), Scalar(4));
    Matrix moved(std::move(copy));
    EXPECT_EQ(moved(0, 0), Scalar(7));
    left = std::move(moved);
    EXPECT_EQ(left(0, 0), Scalar(7));
    left = sum;
    Eigen::Map<Matrix> mapped(left.data(), Rows, Cols);
    EXPECT_EQ(Math::dot(mapped, right), left.dot(right));
    Eigen::Ref<Matrix> reference(left);
    reference(0, 0) = Scalar(5);
    EXPECT_EQ(left(0, 0), Scalar(5));
    const Eigen::Map<const Matrix> readOnly(left.data(), Rows, Cols);
    EXPECT_EQ(readOnly(0, 0), Scalar(5));
    Parent copied = Parent::Constant(Rows, Cols, Scalar(8));
    Matrix owned(copied);
    copied(0, 0) = Scalar(12);
    EXPECT_EQ(owned(0, 0), Scalar(8));
    Parent native = Parent::Constant(Rows, Cols, Scalar(9));
    const auto* data = native.data();
    Matrix converted(std::move(native));
    EXPECT_EQ(converted(0, 0), Scalar(9));
    if constexpr (Matrix::SizeAtCompileTime == Eigen::Dynamic)
    {
      EXPECT_EQ(converted.data(), data);
    }
    Parent assigned = Parent::Constant(Rows, Cols, Scalar(11));
    data = assigned.data();
    converted = std::move(assigned);
    EXPECT_EQ(converted(0, 0), Scalar(11));
    if constexpr (Matrix::SizeAtCompileTime == Eigen::Dynamic)
    {
      EXPECT_EQ(converted.data(), data);
    }
    Math::Matrix<Scalar> identity = Math::Matrix<Scalar>::Identity(Rows, Rows);
    Matrix solution = identity.fullPivLu().solve(right);
    EXPECT_EQ(solution, right);
    solution.noalias() = identity * left;
    EXPECT_EQ(solution, left);
  }

  /// @brief Empty shapes and mixed scalar precisions retain the same contraction contract.
  TEST(MatrixProducts, EmptyShapesAndMixedPrecision)
  {
    for (int rows = 0; rows <= 3; ++rows)
    {
      for (int cols = 0; cols <= 3; ++cols)
      {
        if (rows != 0 && cols != 0)
          continue;
        const Math::Matrix<Complex> empty(rows, cols);
        EXPECT_EQ(empty.dot(empty), Complex(0));
      }
    }
    const Math::FixedSizeMatrix<Complex, 0, 0> zero;
    const Math::FixedSizeMatrix<Complex, 0, 2> noRows;
    const Math::FixedSizeMatrix<Complex, 2, 0> noCols;
    EXPECT_EQ(zero.dot(zero), Complex(0));
    EXPECT_EQ(noRows.dot(noRows), Complex(0));
    EXPECT_EQ(noCols.dot(noCols), Complex(0));
    Math::Matrix<std::complex<float>> left(1, 1);
    Math::FixedSizeMatrix<std::complex<long double>, 1, 1> right;
    left(0, 0) = std::complex<float>(1, 2);
    right(0, 0) = std::complex<long double>(2, -3);
    static_assert(std::is_same_v<decltype(left.dot(right)), std::complex<long double>>);
    EXPECT_EQ(left.dot(right), std::complex<long double>(-4, 7));
    EXPECT_EQ(right.dot(left), std::complex<long double>(-4, -7));
    const Math::Matrix<double> real = Math::Matrix<double>::Constant(1, 1, 2);
    EXPECT_EQ(real.dot(left), Complex(2, -4));
    EXPECT_EQ(left.dot(real), Complex(2, 4));
    Math::FixedSizeVector<Complex, 1> vector;
    vector(0) = Complex(1, 2);
    EXPECT_EQ(vector.dot(right), std::complex<long double>(-4, 7));
    EXPECT_EQ(right.dot(vector), std::complex<long double>(-4, -7));
  }

  /// @brief Eigen decompositions accept Rodin owners as their explicit matrix type.
  TEST(MatrixProducts, NativeEigenDecompositionsAcceptRodinMatrices)
  {
    Math::Matrix<Complex> matrix = Math::Matrix<Complex>::Identity(3, 3);
    matrix(0, 0) = Complex(4);
    matrix(1, 1) = Complex(5);
    matrix(2, 2) = Complex(6);
    matrix(0, 1) = Complex(1, 1);
    matrix(1, 0) = Complex(1, -1);
    Math::Vector<Complex> rhs(3);
    rhs << Complex(1, 2), Complex(3, -1), Complex(-2, 1);
    Eigen::FullPivLU<Math::Matrix<Complex>> lu(matrix);
    const Math::Vector<Complex> direct = lu.solve(rhs);
    EXPECT_LT((matrix * direct - rhs).norm(), 1e-13);
    Eigen::LLT<Math::Matrix<Complex>> llt(matrix);
    ASSERT_EQ(llt.info(), Eigen::Success);
    const Math::Vector<Complex> cholesky = llt.solve(rhs);
    EXPECT_LT((matrix * cholesky - rhs).norm(), 1e-13);
    Eigen::ConjugateGradient<Math::Matrix<Complex>, Eigen::Lower | Eigen::Upper> cg;
    cg.compute(matrix);
    const Math::Vector<Complex> iterative = cg.solve(rhs);
    ASSERT_EQ(cg.info(), Eigen::Success);
    EXPECT_LT((matrix * iterative - rhs).norm(), 1e-12);
    const Math::FixedSizeMatrix<Complex, 3, 3> fixed = matrix;
    Eigen::JacobiSVD<Math::FixedSizeMatrix<Complex, 3, 3>> svd(
      fixed, Eigen::ComputeFullU | Eigen::ComputeFullV);
    EXPECT_LT(
      (svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().adjoint() -
        fixed)
        .norm(),
      1e-12);
  }

  /// @brief Rectangular SVD and QR keep references to the Rodin matrix storage.
  TEST(MatrixProducts, RectangularEigenDecompositions)
  {
    Math::Matrix<Real> matrix(40, 27);
    for (Eigen::Index i = 0; i < matrix.rows(); ++i)
    {
      for (Eigen::Index j = 0; j < matrix.cols(); ++j)
      {
        matrix(i, j) = std::sin(Real(i + 1) * Real(j + 2)) + (i == j ? 5 : 0);
      }
    }
#if EIGEN_VERSION_AT_LEAST(3, 4, 90)
    Eigen::BDCSVD<Math::Matrix<Real>, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(
      matrix);
#else
    Eigen::BDCSVD<Math::Matrix<Real>> svd(
      matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
#endif
    ASSERT_EQ(svd.info(), Eigen::Success);
    EXPECT_LT(
      (svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().adjoint() -
        matrix)
        .norm(),
      1e-12);
    const Math::Matrix<Real> expected = Math::Matrix<Real>::Ones(27, 2);
    const Math::Matrix<Real> rhs = matrix * expected;
    const Math::Matrix<Real> solved = svd.solve(rhs);
    EXPECT_LT((solved - expected).norm(), 1e-12);
    Eigen::ColPivHouseholderQR<Math::Matrix<Real>> qr(matrix);
    const Math::Matrix<Real> qrSolved = qr.solve(rhs);
    EXPECT_LT((qrSolved - expected).norm(), 1e-12);
  }

  /// @brief Coefficient constructors, permutations and native Eigen references interoperate.
  TEST(MatrixProducts, ConstructorsAndNativeReferences)
  {
    const Math::Matrix<Real> dynamic{{1, 2, 3}, {4, 5, 6}};
    const Math::FixedSizeMatrix<Real, 2, 3> fixed{{1, 2, 3}, {4, 5, 6}};
    EXPECT_EQ(dynamic, fixed);
    const Real data[]{1, 4, 2, 5, 3, 6};
    const Math::FixedSizeMatrix<Real, 2, 3> copied(data);
    EXPECT_EQ(copied, fixed);
    const Math::FixedSizeMatrix<Real, 1, 2> two(1, 2);
    const Math::FixedSizeMatrix<Real, 1, 3> three(1, 2, 3);
    const Math::FixedSizeMatrix<Real, 4, 1> four(1, 2, 3, 4);
    const Math::FixedSizeMatrix<Real, 1, 1> one(2);
    EXPECT_EQ(two(1), 2);
    EXPECT_EQ(three(2), 3);
    EXPECT_EQ(four(3), 4);
    const Math::FixedSizeMatrix<Real, 5, 1> five(1, 2, 3, 4, 5);
    EXPECT_EQ(five(4), 5);
    const Eigen::Rotation2D<Real> rotation(Real(0.3));
    Math::FixedSizeMatrix<Real, 2, 2> rotated(rotation);
    EXPECT_EQ(rotated, rotation.toRotationMatrix());
    rotated = rotation;
    EXPECT_EQ(rotated, rotation.toRotationMatrix());
    EXPECT_EQ(one(0), 2);
    Math::Matrix<Real> matrix = dynamic;
    Eigen::Ref<Eigen::MatrixXd> reference(matrix);
    Eigen::Ref<const Eigen::MatrixXd> readOnly(matrix);
    reference(0, 0) = 9;
    EXPECT_EQ(matrix(0, 0), 9);
    EXPECT_EQ(readOnly.data(), matrix.data());
    Eigen::PermutationMatrix<2> permutation;
    permutation.indices() << 1, 0;
    const Math::Matrix<Real> permuted = permutation * matrix;
    EXPECT_EQ(permuted.row(0), matrix.row(1));
    Math::Matrix<Real> identity = permutation;
    EXPECT_EQ(identity, permutation.toDenseMatrix().cast<Real>());
    static_assert(std::is_same_v<decltype(identity = permutation), Math::Matrix<Real>&>);
    identity = permutation;
    EXPECT_EQ(identity, permutation.toDenseMatrix().cast<Real>());
  }

  /// @brief Conservative resizing matches native storage preservation and append behavior.
  TEST(MatrixProducts, ConservativeResizeMatchesEigen)
  {
    using Matrix = Math::Matrix<Complex>;
    using Native = Eigen::MatrixXcd;
    for (Eigen::Index rows = 0; rows <= 4; ++rows)
    {
      for (Eigen::Index cols = 0; cols <= 4; ++cols)
      {
        for (Eigen::Index newRows = 0; newRows <= 4; ++newRows)
        {
          for (Eigen::Index newCols = 0; newCols <= 4; ++newCols)
          {
            Matrix matrix = Matrix::Constant(rows, cols, Complex(3, 2));
            Native native = matrix;
            matrix.conservativeResize(newRows, newCols);
            native.conservativeResize(newRows, newCols);
            const auto commonRows = std::min(rows, newRows);
            const auto commonCols = std::min(cols, newCols);
            EXPECT_EQ(matrix.rows(), native.rows());
            EXPECT_EQ(matrix.cols(), native.cols());
            EXPECT_EQ(matrix.topLeftCorner(commonRows, commonCols),
              native.topLeftCorner(commonRows, commonCols));
            matrix.setConstant(Complex(3, 2));
            native.setConstant(Complex(3, 2));
            const Matrix other = Matrix::Constant(rows, cols, Complex(5, -7));
            matrix.conservativeResizeLike(other);
            native.conservativeResizeLike(other);
            EXPECT_EQ(matrix, native);
          }
        }
      }
    }
    Matrix matrix = Matrix::Constant(2, 3, Complex(1, 2));
    matrix.conservativeResize(4, Eigen::NoChange);
    EXPECT_EQ(matrix.topRows(2), Matrix::Constant(2, 3, Complex(1, 2)));
    matrix.conservativeResize(Eigen::NoChange, 1);
    EXPECT_EQ(matrix.topRows(2), Matrix::Constant(2, 1, Complex(1, 2)));
  }

  /// @brief Blocked Hermitian decompositions preserve complex arithmetic and triangular views.
  TEST(MatrixProducts, BlockedComplexEigenDecompositions)
  {
    using Matrix = Math::Matrix<Complex>;
    Matrix seed(40, 40);
    for (Eigen::Index i = 0; i < seed.rows(); ++i)
    {
      for (Eigen::Index j = 0; j < seed.cols(); ++j)
      {
        seed(i, j) = Complex(std::sin(Real(i + j)), std::cos(Real(i - j)));
      }
    }
    const Matrix matrix = seed.adjoint() * seed + Matrix::Identity(40, 40);
    const Matrix rhs = Matrix::Ones(40, 2);
    Eigen::LLT<Matrix> lower(matrix);
    Eigen::LLT<Matrix, Eigen::Upper> upper(matrix);
    Eigen::LDLT<Matrix> ldlt(matrix);
    ASSERT_EQ(lower.info(), Eigen::Success);
    ASSERT_EQ(upper.info(), Eigen::Success);
    ASSERT_EQ(ldlt.info(), Eigen::Success);
    EXPECT_LT((matrix * lower.solve(rhs) - rhs).norm(), 1e-11);
    EXPECT_LT((matrix * upper.solve(rhs) - rhs).norm(), 1e-11);
    EXPECT_LT((matrix * ldlt.solve(rhs) - rhs).norm(), 1e-11);
    Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(matrix);
    ASSERT_EQ(eigensolver.info(), Eigen::Success);
    EXPECT_LT((matrix * eigensolver.eigenvectors() -
                eigensolver.eigenvectors() * eigensolver.eigenvalues().asDiagonal())
                .norm(),
      1e-10);
  }

}
