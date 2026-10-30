/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <array>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include <Rodin/Variational/P1/P1Element.h>

namespace Rodin::Tests::Unit
{
  template <class Scalar>
  class P1HigherDerivative : public ::testing::Test
  {};

  using P1DerivativeScalars = ::testing::Types<Real, Complex>;
  TYPED_TEST_SUITE(P1HigherDerivative, P1DerivativeScalars);

  /// @brief Pure reference derivatives agree with analytic and finite-difference oracles.
  TYPED_TEST(P1HigherDerivative, PureCoordinateDerivatives)
  {
    using Scalar = TypeParam;
    // Interior samples avoid the pyramid apex, where rational derivatives need
    // not have a direction-independent limit. Coordinates are dimensionless.
    constexpr Real DifferenceStep = 1e-6;
    constexpr Real DifferenceTolerance = 1e-7;
    constexpr Real AnalyticTolerance = 1e-11;
    constexpr std::array<Real, 5> Factorials{1, 1, 2, 6, 24};
    for (const auto geometry : {Geometry::Polytope::Type::Segment,
           Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Quadrilateral,
           Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Pyramid,
           Geometry::Polytope::Type::Hexahedron, Geometry::Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const Geometry::Polytope::Traits traits(geometry);
      const Variational::P1Element<Scalar> element(geometry);
      std::vector<Math::SpatialPoint> samples{traits.getCentroid()};
      if (geometry == Geometry::Polytope::Type::Pyramid)
      {
        samples.emplace_back(Math::SpatialPoint{0.1, 0.2, 0.35});
        samples.emplace_back(Math::SpatialPoint{0.1, 0.15, 0.6});
      }
      for (const auto& point : samples)
      {
        for (size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
        {
          SCOPED_TRACE(vertex);
          const auto basis = element.getBasis(vertex);
          for (size_t coordinate = 0; coordinate < traits.getDimension(); ++coordinate)
          {
            SCOPED_TRACE(coordinate);
            const auto check = [&]<size_t Order>() {
              SCOPED_TRACE(Order);
              const auto derivative = basis.template getDerivative<Order>(coordinate);
              const Scalar actual = derivative(point);
              if constexpr (Order == 0)
              {
                EXPECT_EQ(actual, basis(point));
              }
              else
              {
                auto plus = point, minus = point;
                plus(coordinate) += DifferenceStep;
                minus(coordinate) -= DifferenceStep;
                const auto previous = basis.template getDerivative<Order - 1>(coordinate);
                const Scalar difference =
                  (previous(plus) - previous(minus)) / (2 * DifferenceStep);
                EXPECT_LE(std::abs(actual - difference), DifferenceTolerance);
                if constexpr (Order >= 2)
                {
                  // All nonpyramid bases are linear in each coordinate. The
                  // pyramid's four base functions contain +/- xy/(1-z), whose
                  // Order-th z derivative is +/- Order! xy/(1-z)^(Order+1).
                  Real expected = 0;
                  if (geometry == Geometry::Polytope::Type::Pyramid && coordinate == 2 &&
                    vertex < 4)
                  {
                    const Real sign = vertex % 2 == 0 ? 1 : -1;
                    expected = sign * Factorials[Order] * point.x() * point.y() /
                      std::pow(1 - point.z(), Order + 1);
                  }
                  EXPECT_LE(std::abs(actual - Scalar(expected)), AnalyticTolerance);
                }
              }
            };
            check.template operator()<0>();
            check.template operator()<1>();
            check.template operator()<2>();
            check.template operator()<3>();
            check.template operator()<4>();
          }
        }
      }
    }
    const Variational::P1Element<Scalar> pyramid(Geometry::Polytope::Type::Pyramid);
    const Math::SpatialPoint apex{0, 0, 1};
    for (size_t vertex = 0; vertex < 5; ++vertex)
    {
      for (size_t coordinate = 0; coordinate < 3; ++coordinate)
      {
        const auto basis = pyramid.getBasis(vertex);
        EXPECT_EQ(basis.template getDerivative<2>(coordinate)(apex), Scalar(0));
        EXPECT_EQ(basis.template getDerivative<3>(coordinate)(apex), Scalar(0));
        EXPECT_EQ(basis.template getDerivative<4>(coordinate)(apex), Scalar(0));
      }
    }
  }

  /// @brief Every vector component retains the nonzero rational pyramid derivatives.
  TYPED_TEST(P1HigherDerivative, VectorPyramidDerivatives)
  {
    using Scalar = TypeParam;
    constexpr size_t Components = 3;
    constexpr Real Tolerance = 1e-11;
    constexpr std::array<Real, 5> Factorials{1, 1, 2, 6, 24};
    const Math::SpatialPoint point{0.1, 0.2, 0.35};
    const Variational::P1Element<Math::SpatialVector<Scalar>> element(
      Geometry::Polytope::Type::Pyramid, Components);
    const auto check = [&]<size_t Order>() {
      SCOPED_TRACE(Order);
      for (size_t vertex = 0; vertex < 5; ++vertex)
      {
        SCOPED_TRACE(vertex);
        for (size_t component = 0; component < Components; ++component)
        {
          SCOPED_TRACE(component);
          const auto& basis = element.getBasis(Components * vertex + component);
          for (size_t direction = 0; direction < 3; ++direction)
          {
            for (size_t row = 0; row < Components; ++row)
            {
              SCOPED_TRACE(row);
              Real expected = 0;
              if (direction == 2 && vertex < 4 && row == component)
              {
                const Real sign = vertex % 2 == 0 ? 1 : -1;
                expected = sign * Factorials[Order] * point.x() * point.y() /
                  std::pow(1 - point.z(), Order + 1);
              }
              const auto derivative = basis.template getDerivative<Order>(row, direction);
              ASSERT_LE(std::abs(derivative(point) - Scalar(expected)), Tolerance);
              Scalar output(101);
              derivative(output, point);
              ASSERT_LE(std::abs(output - Scalar(expected)), Tolerance);
            }
          }
        }
      }
    };
    check.template operator()<2>();
    check.template operator()<3>();
    check.template operator()<4>();
  }

  /// @brief Rectangular matrix bases preserve active entries and exact inactive zeros.
  TYPED_TEST(P1HigherDerivative, MatrixPyramidDerivatives)
  {
    using Scalar = TypeParam;
    constexpr size_t Rows = 2, Columns = 3;
    constexpr Real Tolerance = 1e-11;
    constexpr std::array<Real, 5> Factorials{1, 1, 2, 6, 24};
    const Math::SpatialPoint point{0.1, 0.2, 0.35};
    const Variational::P1Element<Math::SpatialMatrix<Scalar>> element(
      Geometry::Polytope::Type::Pyramid, Rows, Columns);
    const auto check = [&]<size_t Order>() {
      SCOPED_TRACE(Order);
      for (size_t vertex = 0; vertex < 5; ++vertex)
      {
        SCOPED_TRACE(vertex);
        for (size_t component = 0; component < Rows * Columns; ++component)
        {
          SCOPED_TRACE(component);
          const auto basis = element.getBasis(Rows * Columns * vertex + component);
          for (size_t direction = 0; direction < 3; ++direction)
          {
            const auto value = basis.template getDerivative<Order>(direction)(point);
            ASSERT_EQ(value.rows(), Rows);
            ASSERT_EQ(value.cols(), Columns);
            for (size_t row = 0; row < Rows; ++row)
            {
              for (size_t column = 0; column < Columns; ++column)
              {
                Real expected = 0;
                if (direction == 2 && vertex < 4 && row * Columns + column == component)
                {
                  const Real sign = vertex % 2 == 0 ? 1 : -1;
                  expected = sign * Factorials[Order] * point.x() * point.y() /
                    std::pow(1 - point.z(), Order + 1);
                }
                ASSERT_LE(std::abs(value(row, column) - Scalar(expected)), Tolerance);
              }
            }
          }
        }
      }
    };
    check.template operator()<2>();
    check.template operator()<3>();
    check.template operator()<4>();
  }
}
