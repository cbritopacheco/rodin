/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <type_traits>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief Arithmetic results own dense coefficients after local operands expire.
  TEST(ExpressionLifetime, DenseArithmetic)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = mesh.getCell();
    const Point point(*cell, Math::SpatialPoint{0.2, 0.3});
    const IntegrationPoint ip(point);
    auto function = VectorFunction(
      size_t(4), [](const Point&) { return Math::Vector<Real>{{2, 4, 6, 8}}; });
    const auto check = [&]<class Location>(const Location& location) {
      const auto sum = (function + function).getValue(location);
      const auto quotient = (function / Real(2)).getValue(location);
      constexpr bool ownsSum =
        std::is_same_v<std::decay_t<decltype(sum)>, typename decltype(sum)::PlainObject>;
      constexpr bool ownsQuotient = std::is_same_v<std::decay_t<decltype(quotient)>,
        typename decltype(quotient)::PlainObject>;
      ASSERT_TRUE(ownsSum);
      ASSERT_TRUE(ownsQuotient);
      for (size_t i = 0; i < 4; ++i)
      {
        EXPECT_EQ(sum[i], Real(4 * (i + 1)));
        EXPECT_EQ(quotient[i], Real(i + 1));
      }
    };
    check(point);
    check(ip);
  }

  /// @brief Both interface evaluation overloads materialize dense trace results.
  TEST(ExpressionLifetime, DenseInterfaceTraces)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 2);
    auto face = mesh.getPolytope(1);
    for (; face; ++face)
    {
      if (mesh.getConnectivity().getIncidence({1, 2}, face->getIndex()).size() == 2)
        break;
    }
    ASSERT_TRUE(face);
    const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
    const IntegrationPoint ip(point);
    auto function = VectorFunction(size_t(4), [](const Point& location) {
      return Math::Vector<Real>{{2 + location.x(), 4 + location.y(), 6, 8}};
    });
    const auto check = [&]<class Location>(const Location& location) {
      const auto average = Average(function).getValue(location);
      const auto jump = Jump(function).getValue(location);
      constexpr bool ownsAverage = std::is_same_v<std::decay_t<decltype(average)>,
        typename decltype(average)::PlainObject>;
      constexpr bool ownsJump = std::is_same_v<std::decay_t<decltype(jump)>,
        typename decltype(jump)::PlainObject>;
      ASSERT_TRUE(ownsAverage);
      ASSERT_TRUE(ownsJump);
      EXPECT_NEAR(average[0], 2 + point.x(), 1e-14);
      EXPECT_NEAR(average[1], 4 + point.y(), 1e-14);
      EXPECT_EQ(average[2], 6);
      EXPECT_EQ(average[3], 8);
      EXPECT_NEAR(jump.norm(), 0, 1e-14);
    };
    check(point);
    check(ip);
  }
}

/// @brief Pointwise wrappers forward side selection to their owned derivative operands.
TEST(ExpressionLifetime, WrappedDerivativeTrace)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(1, 2);
  for (auto cell = mesh.getCell(); cell; ++cell)
  {
    const Point center(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
    mesh.setAttribute({2, cell->getIndex()}, center.x() + center.y() > 1 ? 2 : 1);
  }
  P1<Real> scalarSpace(mesh);
  P1<Math::SpatialVector<Real>> vectorSpace(mesh, 2);
  GridFunction scalar(scalarSpace);
  GridFunction vector(vectorSpace);
  scalar = RealFunction([](const Point& p) {
    const Real sum = p.x() + p.y();
    return sum > 1 ? 3 * sum - 2 : sum;
  });
  vector = VectorFunction(size_t(2),
    [&scalar](const Point& p) { return Math::SpatialVector<Real>{scalar(p), 0}; });
  auto face = mesh.getPolytope(1);
  for (; face; ++face)
  {
    if (mesh.getConnectivity().getIncidence({1, 2}, face->getIndex()).size() == 2)
      break;
  }
  ASSERT_TRUE(face);
  const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
  const IntegrationPoint ip(point);
  for (const size_t attribute : {1u, 2u})
  {
    const Real expected = attribute == 1 ? 1 : 3;
    auto component = Component(Grad(scalar), size_t(0));
    auto transpose = Transpose(Jacobian(vector));
    auto conjugate = Conjugate(Div(vector));
    component.traceOf(attribute);
    transpose.traceOf(attribute);
    conjugate.traceOf(attribute);
    const auto check = [&](const auto& location) {
      EXPECT_NEAR(component.getValue(location), expected, 1e-13);
      EXPECT_NEAR(transpose.getValue(location)(0, 0), expected, 1e-13);
      EXPECT_NEAR(conjugate.getValue(location), expected, 1e-13);
      const auto copied = component;
      EXPECT_NEAR(copied.getValue(location), expected, 1e-13);
    };
    check(point);
    check(ip);
  }
}

/// @brief Spatial conjugation returns an owning value for real and complex ranges.
TEST(ExpressionLifetime, SpatialConjugation)
{
  const auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  const Point point(*mesh.getCell(), Math::SpatialPoint{0.5});
  const auto vector = VectorFunction(size_t(2), [](const Point&) {
    return Math::SpatialVector<Complex>{Complex(2, 3), Complex(5, -7)};
  });
  const auto matrix = MatrixFunction(size_t(2), size_t(1), [](const Point&) {
    Math::SpatialMatrix<Complex> result(2, 1);
    result(0, 0) = Complex(11, 13);
    result(1, 0) = Complex(17, -19);
    return result;
  });
  const auto conjugateVector = Conjugate(vector).getValue(point);
  const auto conjugateMatrix = Conjugate(matrix).getValue(point);
  EXPECT_EQ(conjugateVector[0], Complex(2, -3));
  EXPECT_EQ(conjugateVector[1], Complex(5, 7));
  EXPECT_EQ(conjugateMatrix(0, 0), Complex(11, -13));
  EXPECT_EQ(conjugateMatrix(1, 0), Complex(17, 19));
}
