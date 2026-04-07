#include <gtest/gtest.h>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

class LinearSystemAdvancedTest : public ::testing::Test
{};

TEST_F(LinearSystemAdvancedTest, SparseEliminateWithOffset)
{
  using LSType = LinearSystem<SparseMatrix<Real>, Vector<Real>>;
  LSType ls;
  ls.getOperator().resize(4, 4);
  ls.getVector().resize(4);
  ls.getVector() << 1.0, 2.0, 3.0, 4.0;

  std::vector<Eigen::Triplet<Real>> t;
  t.emplace_back(0, 0, 4.0);
  t.emplace_back(1, 1, 3.0);
  t.emplace_back(2, 2, 2.0);
  t.emplace_back(3, 3, 5.0);
  t.emplace_back(1, 2, 1.0);
  t.emplace_back(2, 1, 1.0);
  ls.getOperator().setFromTriplets(t.begin(), t.end());

  IndexMap<Real> dofs;
  dofs[1] = 7.0; // with offset=1 this constrains global index 2
  ls.eliminate(dofs, 1);

  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(2, 2), 1.0);
  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(2, 1), 0.0);
  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(ls.getVector().coeff(2), 7.0);
}

TEST_F(LinearSystemAdvancedTest, DenseMergeImposesParentConstraint)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;
  LSType ls;
  ls.getOperator().setIdentity(3, 3);
  ls.getVector().setOnes(3);

  IndexMap<std::pair<IndexArray, Vector<Real>>> dofs;
  IndexArray children(2);
  children << 0, 1;
  Vector<Real> coeffs(2);
  coeffs << 0.5, 0.5;
  dofs[2] = {children, coeffs};

  ls.merge(dofs);

  EXPECT_DOUBLE_EQ(ls.getOperator()(2, 2), 1.0);
  EXPECT_DOUBLE_EQ(ls.getOperator()(2, 0), -0.5);
  EXPECT_DOUBLE_EQ(ls.getOperator()(2, 1), -0.5);
}

TEST_F(LinearSystemAdvancedTest, SparseMergeImposesParentConstraint)
{
  using LSType = LinearSystem<SparseMatrix<Real>, Vector<Real>>;
  LSType ls;
  ls.getOperator().resize(3, 3);
  ls.getVector().setOnes(3);

  std::vector<Eigen::Triplet<Real>> t;
  t.emplace_back(0, 0, 1.0);
  t.emplace_back(1, 1, 1.0);
  t.emplace_back(2, 2, 1.0);
  ls.getOperator().setFromTriplets(t.begin(), t.end());

  IndexMap<std::pair<IndexArray, Vector<Real>>> dofs;
  IndexArray children(2);
  children << 0, 1;
  Vector<Real> coeffs(2);
  coeffs << 0.25, 0.75;
  dofs[2] = {children, coeffs};

  ls.merge(dofs);

  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(2, 2), 1.0);
  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(2, 0), -0.25);
  EXPECT_DOUBLE_EQ(ls.getOperator().coeff(2, 1), -0.75);
}

TEST_F(LinearSystemAdvancedTest, TupleGetSupportsRvalueAndConstRvalue)
{
  using LSType = LinearSystem<Matrix<Real>, Vector<Real>>;
  LSType ls;
  ls.getOperator().setIdentity(2, 2);
  ls.getSolution().setOnes(2);
  ls.getVector().setZero(2);

  auto&& op = Math::get<0>(std::move(ls));
  EXPECT_EQ(op.rows(), 2);

  const LSType cls{};
  auto&& rhs = Math::get<2>(std::move(cls));
  EXPECT_EQ(rhs.size(), 0);
}

TEST_F(LinearSystemAdvancedTest, CRTPBaseDispatchDenseAndSparse)
{
  using DenseLS = LinearSystem<Matrix<Real>, Vector<Real>>;
  using DenseBase = LinearSystemBase<Matrix<Real>, Vector<Real>, DenseLS>;
  DenseLS dls;
  DenseBase& db = dls;
  db.getOperator().resize(2, 2);
  EXPECT_EQ(db.getOperator().rows(), 2);

  using SparseLS = LinearSystem<SparseMatrix<Real>, Vector<Real>>;
  using SparseBase = LinearSystemBase<SparseMatrix<Real>, Vector<Real>, SparseLS>;
  SparseLS sls;
  SparseBase& sb = sls;
  sb.getOperator().resize(2, 2);
  EXPECT_EQ(sb.getOperator().rows(), 2);
}
