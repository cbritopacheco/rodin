/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Sparse elimination must enforce the complete constrained system.
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <random>
#include "Rodin/Math/LinearSystem.h"

using namespace Rodin;

/// @brief Omitting a zero transpose entry must not change Dirichlet elimination.
TEST(SparseElimination, IndependentOfStoredTransposeZero)
{
  for (const bool storeZero : {false, true})
  {
    SCOPED_TRACE(storeZero);
    Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>> system;
    std::vector<Eigen::Triplet<Real>> entries{{0, 0, 2}, {0, 1, 3}, {1, 1, 4}};
    if (storeZero)
      entries.emplace_back(1, 0, 0);
    system.getOperator().resize(2, 2);
    system.getOperator().setFromTriplets(entries.begin(), entries.end());
    system.getVector().resize(2);
    system.getVector() << 5, 6;
    IndexMap<Real> fixed;
    fixed.emplace(0, 1);
    system.eliminate(fixed);
    Math::Matrix<Real> expected(2, 2);
    expected << 1, 0, 0, 4;
    const Math::Matrix<Real> actual = system.getOperator();
    EXPECT_EQ((actual - expected).norm(), 0);
    Math::Vector<Real> load(2), solution(2);
    load << 1, 6;
    solution << 1, Real(1.5);
    EXPECT_EQ((system.getVector() - load).norm(), 0);
    EXPECT_NEAR(
      (actual.fullPivLu().solve(system.getVector()) - solution).norm(), 0, 1e-12);
  }
}

/// @brief A missing constrained diagonal must become an identity constraint row.
TEST(SparseElimination, InsertsMissingConstrainedDiagonal)
{
  Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>> system;
  std::vector<Eigen::Triplet<Real>> entries{{0, 1, 1}, {1, 0, 1}, {1, 1, 2}};
  system.getOperator().resize(2, 2);
  system.getOperator().setFromTriplets(entries.begin(), entries.end());
  system.getVector().resize(2);
  system.getVector() << 5, 6;
  IndexMap<Real> fixed;
  fixed.emplace(0, 1);
  system.eliminate(fixed);
  Math::Matrix<Real> expected(2, 2);
  expected << 1, 0, 0, 2;
  const Math::Matrix<Real> actual = system.getOperator();
  EXPECT_EQ((actual - expected).norm(), 0);
  Math::Vector<Real> load(2), solution(2);
  load << 1, 5;
  solution << 1, Real(2.5);
  EXPECT_EQ((system.getVector() - load).norm(), 0);
  EXPECT_NEAR((actual.fullPivLu().solve(system.getVector()) - solution).norm(), 0, 1e-12);
}

namespace
{
  /**
   * @brief Compares elimination with an independent dense row/column oracle.
   * @tparam MatrixScalar Operator coefficient type.
   * @tparam VectorScalar Load coefficient type.
   * @tparam DOFScalar Prescribed coefficient type.
   * @param matrix Original sparse matrix, possibly uncompressed.
   * @param load Original right-hand side, independent of the constraints.
   * @param fixed Local constraint map.
   * @param offset Shift applied to the map's indices.
   */
  template <class MatrixScalar, class VectorScalar, class DOFScalar>
  void checkElimination(const Math::SparseMatrix<MatrixScalar>& matrix,
    const Math::Vector<VectorScalar>& load, const IndexMap<DOFScalar>& fixed,
    size_t offset = 0)
  {
    Math::LinearSystem<Math::SparseMatrix<MatrixScalar>, Math::Vector<VectorScalar>>
      system;
    system.getOperator() = matrix;
    system.getVector() = load;
    system.getSolution() = Math::Vector<VectorScalar>::Ones(load.size());
    Math::Matrix<MatrixScalar> expected = matrix;
    Math::Vector<VectorScalar> prescribed = Math::Vector<VectorScalar>::Zero(load.size());
    size_t absent = 0;
    for (const auto& [index, value] : fixed)
    {
      const Index column = index + offset;
      prescribed[column] = value;
      bool found = false;
      for (typename Math::SparseMatrix<MatrixScalar>::InnerIterator it(matrix, column);
           it; ++it)
        found |= it.row() == static_cast<Eigen::Index>(column);
      absent += !found;
    }
    Math::Vector<VectorScalar> expectedLoad = load - matrix * prescribed;
    for (const auto& [index, value] : fixed)
    {
      const Index column = index + offset;
      expected.row(column).setZero();
      expected.col(column).setZero();
      expected(column, column) = MatrixScalar(1);
      expectedLoad[column] = value;
    }
    EXPECT_EQ(&system.eliminate(fixed, offset), &system);
    const Math::Matrix<MatrixScalar> actual = system.getOperator();
    EXPECT_LT((actual - expected).norm(), 1e-11 * std::max(Real(1), expected.norm()));
    EXPECT_LT((system.getVector() - expectedLoad).norm(),
      1e-11 * std::max(Real(1), expectedLoad.norm()));
    EXPECT_TRUE(system.getOperator().isCompressed());
    EXPECT_EQ(system.getOperator().nonZeros(), matrix.nonZeros() + absent);
    EXPECT_TRUE(system.getSolution().isOnes());
    for (Eigen::Index column = 0; column < matrix.cols(); ++column)
    {
      for (typename Math::SparseMatrix<MatrixScalar>::InnerIterator original(
             matrix, column);
           original; ++original)
      {
        const auto* rows = system.getOperator().innerIndexPtr();
        const auto* columns = system.getOperator().outerIndexPtr();
        EXPECT_TRUE(std::binary_search(
          rows + columns[column], rows + columns[column + 1], original.row()));
      }
      Optional<Eigen::Index> previous;
      for (typename Math::SparseMatrix<MatrixScalar>::InnerIterator it(
             system.getOperator(), column);
           it; ++it)
      {
        if (previous)
        {
          EXPECT_LT(*previous, it.row());
        }
        previous = it.row();
      }
    }
    const auto once = system.getVector();
    system.eliminate(fixed, offset);
    EXPECT_LT((Math::Matrix<MatrixScalar>(system.getOperator()) - expected).norm(),
      1e-11 * std::max(Real(1), expected.norm()));
    EXPECT_LT((system.getVector() - once).norm(), 1e-11 * std::max(Real(1), once.norm()));
  }
}

/// @brief Independent nonsymmetric patterns, loads, offsets and scalar domains.
TEST(SparseElimination, RandomizedRealComplexAndPromotedLoad)
{
  std::mt19937 random(1942);
  for (size_t sample = 0; sample < 160; ++sample)
  {
    const size_t count = 1 + random() % 71;
    const size_t offset = random() % count;
    std::vector<Math::SparseTriplet<Real>> entries;
    std::vector<Math::SparseTriplet<Complex>> complexEntries;
    Math::Vector<Real> load(count);
    Math::Vector<Complex> complexLoad(count);
    for (size_t column = 0; column < count; ++column)
    {
      load[column] = Real(int(random() % 29) - 14) / 7;
      complexLoad[column] = Complex(load[column], Real(int(random() % 23) - 11) / 5);
      for (size_t row = 0; row < count; ++row)
      {
        if (random() % 7 == 0)
        {
          const Real value = Real(int(random() % 17) - 8) / 3;
          entries.emplace_back(row, column, value);
          complexEntries.emplace_back(row, column, Complex(value, Real(row % 5) / 7));
        }
      }
    }
    Math::SparseMatrix<Real> matrix(count, count);
    Math::SparseMatrix<Complex> complexMatrix(count, count);
    matrix.setFromTriplets(entries.begin(), entries.end());
    complexMatrix.setFromTriplets(complexEntries.begin(), complexEntries.end());
    IndexMap<Real> fixed;
    IndexMap<Complex> complexFixed;
    for (size_t index = 0; index < count - offset; ++index)
    {
      if (random() % 3 == 0)
      {
        fixed.emplace(index, Real(int(random() % 13) - 6) / 2);
        complexFixed.emplace(index, Complex(fixed.at(index), Real(index % 3) / 2));
      }
    }
    SCOPED_TRACE(sample);
    checkElimination(matrix, load, fixed, offset);
    checkElimination(complexMatrix, complexLoad, complexFixed, offset);
    checkElimination(matrix, complexLoad, complexFixed, offset);
  }
}

/// @brief Empty systems and completely constrained storage need no special pattern.
TEST(SparseElimination, EmptyUncompressedAndAllConstraints)
{
  for (size_t count : {size_t(0), size_t(1), size_t(5), size_t(33)})
  {
    Math::SparseMatrix<Real> matrix(count, count);
    Math::Vector<Real> load = Math::Vector<Real>::Ones(count);
    IndexMap<Real> fixed;
    checkElimination(matrix, load, fixed, std::numeric_limits<size_t>::max());
    for (size_t index = 0; index < count; ++index)
      fixed.emplace(index, Real(index) - 2);
    checkElimination(matrix, load, fixed);
    if (count > 1)
    {
      matrix.insert(0, count - 1) = 3;
      ASSERT_FALSE(matrix.isCompressed());
      checkElimination(matrix, load, IndexMap<Real>{});
      checkElimination(matrix, load, fixed);
      IndexMap<Real> one;
      one.emplace(count - 1, 2);
      checkElimination(matrix, load, one);
    }
  }
}

/// @brief Consecutive and fragmented constraints exercise pruning and its fallback.
TEST(SparseElimination, LocalAndNonlocalStorageRanges)
{
  constexpr size_t Count = 192;
  for (bool local : {true, false})
  {
    std::vector<Math::SparseTriplet<Real>> entries;
    for (size_t column = 0; column < Count; ++column)
    {
      if (column % 31 == 0)
        continue;
      for (size_t step = 0; step < 17; ++step)
      {
        const size_t row =
          local ? (column + step + Count - 8) % Count : (column * 19 + step * 11) % Count;
        entries.emplace_back(row, column, Real(step + 1) / 19);
      }
    }
    Math::SparseMatrix<Real> matrix(Count, Count);
    matrix.setFromTriplets(entries.begin(), entries.end());
    const Math::Vector<Real> load = Math::Vector<Real>::LinSpaced(Count, -2, 3);
    for (size_t start : {size_t(0), size_t(60), size_t(160)})
    {
      IndexMap<Real> fixed;
      for (size_t index = start; index < start + 24; ++index)
        fixed.emplace(index, Real(index % 7) / 3);
      checkElimination(matrix, load, fixed);
    }
    IndexMap<Real> all;
    for (size_t index = 0; index < Count; ++index)
      all.emplace(index, Real(index % 7) / 3);
    checkElimination(matrix, load, all);
  }
}

/// @brief Validation precedes all coefficient and load mutations, including overflow.
TEST(SparseElimination, RejectsInvalidDimensionsAndIndicesBeforeMutation)
{
  using System = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  System system;
  system.getOperator().resize(3, 2);
  system.getVector() = Math::Vector<Real>::Ones(3);
  IndexMap<Real> fixed{{0, 2}};
  EXPECT_THROW(system.eliminate(fixed), Alert::Exception);
  system.getOperator().resize(3, 3);
  system.getOperator().setIdentity();
  system.getVector().resize(2);
  EXPECT_THROW(system.eliminate(fixed), Alert::Exception);
  system.getVector() = Math::Vector<Real>::Ones(3);
  for (size_t offset : {size_t(3), std::numeric_limits<size_t>::max()})
    EXPECT_THROW(system.eliminate(fixed, offset), Alert::Exception);
  fixed.emplace(3, 1);
  EXPECT_THROW(system.eliminate(fixed), Alert::Exception);
  EXPECT_TRUE(system.getVector().isOnes());
  EXPECT_EQ(system.getOperator().diagonal().sum(), 3);
}

/// @brief Zero assignment clears nonfinite constrained entries and preserves free bits.
TEST(SparseElimination, NonfiniteAndSignedZeroCoefficientClearing)
{
  constexpr size_t Count = 24;
  for (bool complex : {false, true})
  {
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const Real infinity = std::numeric_limits<Real>::infinity();
    std::vector<Math::SparseTriplet<Real>> entries;
    for (size_t row = 0; row < Count; ++row)
      entries.emplace_back(row, 0, row % 2 ? nan : infinity);
    entries.emplace_back(0, 1, -0.0);
    entries.emplace_back(1, 1, infinity);
    entries.emplace_back(2, 1, nan);
    entries.emplace_back(2, 2, 2);
    entries.emplace_back(0, 2, -0.0);
    Math::SparseMatrix<Real> original(Count, Count);
    original.setFromTriplets(entries.begin(), entries.end());
    IndexMap<Real> fixed{{1, 0}};
    Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>> system;
    system.getOperator() = original;
    system.getVector() = Math::Vector<Real>::Zero(Count);
    system.eliminate(fixed);
    EXPECT_EQ(system.getOperator().coeff(1, 0), 0);
    EXPECT_EQ(system.getOperator().coeff(2, 1), 0);
    EXPECT_FALSE(std::signbit(system.getOperator().coeff(0, 1)));
    EXPECT_TRUE(std::signbit(system.getOperator().coeff(0, 2)));
    EXPECT_TRUE(std::isinf(system.getOperator().coeff(0, 0)));
    EXPECT_TRUE(std::isnan(system.getOperator().coeff(3, 0)));
    if (complex)
    {
      Math::LinearSystem<Math::SparseMatrix<Complex>, Math::Vector<Complex>> other;
      other.getOperator() = original.cast<Complex>();
      other.getVector() = Math::Vector<Complex>::Zero(Count);
      other.eliminate(fixed);
      EXPECT_EQ(other.getOperator().coeff(1, 0), Complex(0));
      EXPECT_EQ(other.getOperator().coeff(2, 1), Complex(0));
    }
  }
}
