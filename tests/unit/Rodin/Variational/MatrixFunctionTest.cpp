/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/*
 * @file MatrixFunctionTest.cpp
 * @brief Unit tests for MatrixFunction class.
 */
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static Geometry::Point makePoint(Mesh<Context::Local>& mesh)
{
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
  return Geometry::Point(*it, rc);
}

/// @brief Verifies constant matrix for variational matrix function by checking tolerance-based numerical results.
TEST(Rodin_Variational_MatrixFunction, ConstantMatrix)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto p = makePoint(mesh);

  Math::Matrix<Real> m(2, 2);
  m << 1.0, 2.0, 3.0, 4.0;
  MatrixFunction mf(m);

  const auto& val = mf.getValue(p);
  EXPECT_DOUBLE_EQ(val(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(val(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(val(1, 0), 3.0);
  EXPECT_DOUBLE_EQ(val(1, 1), 4.0);
}

/// @brief Verifies rows and columns for variational matrix function by checking exact expected values.
TEST(Rodin_Variational_MatrixFunction, RowsAndColumns)
{
  Math::Matrix<Real> m(3, 2);
  m.setZero();
  MatrixFunction mf(m);

  EXPECT_EQ(mf.getRows(), 3u);
  EXPECT_EQ(mf.getColumns(), 2u);
}

/// @brief Verifies copy for variational matrix function by checking tolerance-based numerical results, copy semantics.
TEST(Rodin_Variational_MatrixFunction, Copy)
{
  Math::Matrix<Real> m(2, 2);
  m << 5.0, 6.0, 7.0, 8.0;
  MatrixFunction mf(m);

  auto* mfCopy = mf.copy();
  ASSERT_NE(mfCopy, nullptr);

  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto p = makePoint(mesh);

  const auto& val = mfCopy->getValue(p);
  EXPECT_DOUBLE_EQ(val(0, 0), 5.0);
  EXPECT_DOUBLE_EQ(val(1, 1), 8.0);

  delete mfCopy;
}

/// @brief Verifies get order for variational matrix function by checking exact expected values, true predicates.
TEST(Rodin_Variational_MatrixFunction, GetOrder)
{
  Math::Matrix<Real> m(2, 2);
  m.setIdentity();
  MatrixFunction mf(m);

  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getCell();
  auto polytope = *it;

  auto order = mf.getOrder(polytope);
  ASSERT_TRUE(order.has_value());
  EXPECT_EQ(order.value(), 0u);
}

/// @brief Verifies identity matrix 2 x 2 for variational matrix function by checking tolerance-based numerical results.
TEST(Rodin_Variational_MatrixFunction, IdentityMatrix2x2)
{
  Math::Matrix<Real> m(2, 2);
  m.setIdentity();
  MatrixFunction mf(m);

  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto p = makePoint(mesh);

  const auto& val = mf.getValue(p);
  EXPECT_DOUBLE_EQ(val(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(val(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(val(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(val(1, 1), 1.0);
}

/// @brief Invalid declared extents are rejected before evaluating either scalar callback type.
TEST(Rodin_Variational_MatrixFunction, CallableRejectsInvalidDeclaredDimensions)
{
  const auto check = []<class Scalar>() {
    size_t calls = 0;
    const auto callback = [&](const Point&) {
      ++calls;
      return Math::SpatialMatrix<Scalar>::Identity(2, 2);
    };
    const size_t tooLarge = RODIN_MAXIMAL_SPACE_DIMENSION + 1;
    for (size_t invalid : {size_t{0}, tooLarge, std::numeric_limits<size_t>::max()})
    {
      for (bool invalidRows : {false, true})
      {
        try
        {
          MatrixFunction coefficient(invalidRows ? invalid : size_t{2},
            invalidRows ? size_t{2} : invalid, callback);
          ADD_FAILURE() << "Invalid declared matrix extents were accepted";
        }
        catch (const Alert::Exception& error)
        {
          const std::string message = error.what();
          EXPECT_NE(message.find("member function MatrixFunction"), std::string::npos);
          EXPECT_NE(message.find("Invalid spatial matrix coefficient dimensions"),
            std::string::npos);
        }
        EXPECT_EQ(calls, 0u);
      }
    }
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto point = makePoint(mesh);
    MatrixFunction valid(size_t{2}, size_t{2}, callback);
    const auto value = valid.getValue(point);
    EXPECT_EQ(value(0, 0), Scalar(1));
    EXPECT_EQ(value(1, 1), Scalar(1));
    EXPECT_EQ(value(0, 1), Scalar(0));
    EXPECT_EQ(value(1, 0), Scalar(0));
    EXPECT_EQ(calls, 1u);
  };
  check.template operator()<Real>();
  check.template operator()<Complex>();
}

/// @brief A shape-changing callback raises member context and permits correction and copying.
TEST(Rodin_Variational_MatrixFunction, CallableShapeErrorsPermitRetry)
{
  const auto check = []<class Scalar>() {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto point = makePoint(mesh);
    const IntegrationPoint ip(point);
    size_t rows = 2, columns = 2, calls = 0;
    MatrixFunction coefficient(size_t{2}, size_t{2}, [&](const Point&) {
      ++calls;
      Math::SpatialMatrix<Scalar> value(
        static_cast<std::uint8_t>(rows), static_cast<std::uint8_t>(columns));
      value.setConstant(Scalar(5));
      return value;
    });
    coefficient.setOrder(4);
    for (bool wrongRows : {false, true})
    {
      rows = wrongRows ? 1 : 2;
      columns = wrongRows ? 2 : 1;
      const auto previousCalls = calls;
      try
      {
        (void)coefficient.getValue(ip);
        ADD_FAILURE() << "Callback returned incompatible matrix extents";
      }
      catch (const Alert::Exception& error)
      {
        const std::string message = error.what();
        EXPECT_NE(message.find("member function getValue"), std::string::npos);
        EXPECT_NE(message.find("Matrix coefficient returned incompatible dimensions"),
          std::string::npos);
      }
      EXPECT_EQ(calls, previousCalls + 1);
      EXPECT_EQ(coefficient.getRows(), 2u);
      EXPECT_EQ(coefficient.getColumns(), 2u);
      EXPECT_EQ(coefficient.getOrder(point.getPolytope()), Optional<size_t>(4));
      rows = columns = 2;
      const auto value = coefficient.getValue(point);
      EXPECT_EQ(value.rows(), 2);
      EXPECT_EQ(value.cols(), 2);
      std::unique_ptr<decltype(coefficient)> copied(coefficient.copy());
      EXPECT_EQ(copied->getRows(), 2u);
      EXPECT_EQ(copied->getColumns(), 2u);
      EXPECT_EQ(copied->getOrder(point.getPolytope()), Optional<size_t>(4));
      const auto copiedValue = copied->getValue(ip);
      for (size_t row = 0; row < 2; ++row)
      {
        for (size_t column = 0; column < 2; ++column)
        {
          EXPECT_EQ(value(row, column), Scalar(5));
          EXPECT_EQ(copiedValue(row, column), Scalar(5));
        }
      }
    }
  };
  check.template operator()<Real>();
  check.template operator()<Complex>();
}

namespace Rodin::Tests::Unit
{
  /// @brief Native Eigen constants keep their copied value and right-conjugating contractions.
  TEST(MatrixFunctionContract, ComplexEigenStorageAndTemporaries)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = mesh.getPolytope(mesh.getDimension(), 0);
    Point point(*cell, Math::SpatialPoint{0.25, 0.25});
    Eigen::MatrixX<Complex> native(2, 2);
    native << Complex(1, 2), Complex(3, -1), Complex(-2, 1), Complex(2, -3);
    const auto saved = native;
    MatrixFunction owned(native);
    MatrixFunction temporary((native + native).eval());
    Math::Matrix<Complex> right(2, 2);
    right << Complex(2, -3), Complex(-1, 4), Complex(3, 2), Complex(1, -2);
    MatrixFunction rhs(right);
    Complex expected = 0;
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
        expected += saved(i, j) * std::conj(right(i, j));
    }
    native.setZero();
    right.setZero();
    EXPECT_EQ(owned.getValue(point), saved);
    EXPECT_EQ(Dot(owned, rhs).getValue(point), expected);
    EXPECT_EQ(Dot(temporary, rhs).getValue(point), Complex(2) * expected);
  }
}
