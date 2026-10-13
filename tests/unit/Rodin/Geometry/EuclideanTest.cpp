/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Analytic contracts for the Euclidean primitive templates. */
#include <gtest/gtest.h>
#include <string>
#include <type_traits>
#include <utility>
#include "Rodin/Alert/Exception.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Geometry/Euclidean/Common.h"
#include "Rodin/Geometry/Euclidean/Point2D.h"
#include "Rodin/Geometry/Euclidean/Rectangle.h"

using namespace Rodin;
using namespace Rodin::Geometry::Euclidean;

TEST(Rodin_Geometry_Euclidean, LinesAndSegments)
{
  const Point2D<Real> origin(0, 0), endpoint(3, 4);
  LineSegment2D<Real> segment(origin, endpoint);
  EXPECT_DOUBLE_EQ(segment.length(), 5);
  EXPECT_NEAR(segment.direction().norm(), 1, 1e-12);
  EXPECT_LT((segment(0.5) - Point2D<Real>(1.5, 2)).norm(), 1e-12);
  EXPECT_EQ(segment.reverse().start(), endpoint);
  segment.setStart({1, 0}).setEnd({1, 4});
  EXPECT_DOUBLE_EQ(segment.length(), 4);
  Line2D<Real> vertical(segment);
  Line2D<Real> horizontal({0, 1, 2});
  EXPECT_TRUE(vertical.isVertical());
  EXPECT_FALSE(vertical.slope());
  EXPECT_FALSE(vertical.yIntercept());
  EXPECT_FALSE(vertical(1));
  EXPECT_TRUE(horizontal.isHorizontal());
  EXPECT_FALSE(horizontal.xIntercept());
  EXPECT_DOUBLE_EQ(*vertical.xIntercept(), 1);
  EXPECT_DOUBLE_EQ(*horizontal.yIntercept(), 2);
  EXPECT_DOUBLE_EQ(*horizontal.slope(), 0);
  EXPECT_DOUBLE_EQ(*horizontal(7), 2);
  ASSERT_TRUE(vertical.intersect(horizontal));
  EXPECT_LT((vertical.intersect(horizontal).value() - Point2D<Real>(1, 2)).norm(), 1e-12);
  EXPECT_FALSE(vertical.intersect(Line2D<Real>(1, 0, 3)));
  EXPECT_DOUBLE_EQ(vertical.distance(Line2D<Real>(1, 0, 3)), 2);
  EXPECT_DOUBLE_EQ(horizontal.distance(Line2D<Real>(0, 1, 5)), 3);
  EXPECT_DOUBLE_EQ(vertical.distance(horizontal), 0);
  EXPECT_DOUBLE_EQ(Line2D<Real>(origin, endpoint)(endpoint), 0);
  LineSegment2D<Real> crossing({0, 2}, {2, 2});
  ASSERT_TRUE(segment.intersect(crossing));
  EXPECT_LT((segment.intersect(crossing).value() - Point2D<Real>(1, 2)).norm(), 1e-12);
  EXPECT_FALSE(segment.intersect(LineSegment2D<Real>({2, 0}, {2, 4})));
  EXPECT_FALSE(segment.intersect(LineSegment2D<Real>({3, 5}, {4, 5})));
  EXPECT_FALSE(origin.connect(origin));
  ASSERT_TRUE(origin.connect(endpoint));
  EXPECT_DOUBLE_EQ(origin.connect(endpoint)->length(), 5);
  ASSERT_TRUE(origin.connect(horizontal));
  EXPECT_LT((origin.connect(horizontal)->end() - Point2D<Real>(0, 2)).norm(), 1e-12);
  EXPECT_FALSE(Point2D<Real>(0, 2).connect(horizontal));
}

TEST(Rodin_Geometry_Euclidean, CircleIntersectionsAndClosestPoints)
{
  Circle<Real> circle({2, 3}, 2);
  EXPECT_DOUBLE_EQ(circle({2, 3}), -4);
  EXPECT_DOUBLE_EQ(circle({4, 3}), 0);
  EXPECT_GT(circle({5, 3}), 0);
  const Math::Rad angle(0);
  EXPECT_LT((circle(angle) - Point2D<Real>(4, 3)).norm(), 1e-12);
  auto tangent = circle.tangent(angle);
  EXPECT_NEAR(tangent({4, 3}), 0, 1e-12);
  // Scaling the coefficients must not move the tangent intersection.
  auto touch = Line2D<Real>(2, 0, 8).intersect(circle);
  ASSERT_TRUE(std::holds_alternative<Point2D<Real>>(touch));
  EXPECT_LT((std::get<Point2D<Real>>(touch) - Point2D<Real>(4, 3)).norm(), 1e-12);
  auto two = Line2D<Real>(0, 1, 3).intersect(circle);
  using Pair = std::pair<Point2D<Real>, Point2D<Real>>;
  ASSERT_TRUE(std::holds_alternative<Pair>(two));
  for (const auto& point : {std::get<Pair>(two).first, std::get<Pair>(two).second})
  {
    EXPECT_NEAR(circle(point), 0, 1e-12);
    EXPECT_NEAR(point.y(), 3, 1e-12);
  }
  EXPECT_TRUE(
    std::holds_alternative<std::nullopt_t>(Line2D<Real>(1, 0, 5).intersect(circle)));
  for (Real x : {3.0, 5.0})
  {
    auto connection = Point2D<Real>(x, 3).connect(circle);
    ASSERT_TRUE(connection);
    EXPECT_LT((connection->end() - Point2D<Real>(4, 3)).norm(), 1e-12);
  }
  EXPECT_FALSE(Point2D<Real>(2, 3).connect(circle));
  EXPECT_FALSE(Point2D<Real>(4, 3).connect(circle));
  auto connection = Line2D<Real>(1, 0, 5).connect(circle);
  ASSERT_TRUE(connection);
  EXPECT_LT((connection->start() - Point2D<Real>(5, 3)).norm(), 1e-12);
  EXPECT_LT((connection->end() - Point2D<Real>(4, 3)).norm(), 1e-12);
  EXPECT_FALSE(Line2D<Real>(1, 0, 2).connect(circle));
  EXPECT_FALSE(Line2D<Real>(1, 0, 3).connect(circle));
  circle.setCenter({0, 0}).setRadius(1);
  EXPECT_DOUBLE_EQ(circle.radius(), 1);
  EXPECT_EQ(circle.center(), Point2D<Real>(0, 0));
}

TEST(Rodin_Geometry_Euclidean, RectangleCornersAndInclusiveContainment)
{
  using Rect = Rectangle<Real>;
  Rect rectangle({1, 2}, {4, 6});
  EXPECT_DOUBLE_EQ(rectangle.width(), 3);
  EXPECT_DOUBLE_EQ(rectangle.height(), 4);
  EXPECT_DOUBLE_EQ(rectangle.area(), 12);
  EXPECT_EQ(rectangle.getCorner<Rect::BottomLeft>(), Point2D<Real>(1, 2));
  EXPECT_EQ(rectangle.getCorner<Rect::BottomRight>(), Point2D<Real>(4, 2));
  EXPECT_EQ(rectangle.getCorner<Rect::TopLeft>(), Point2D<Real>(1, 6));
  EXPECT_EQ(rectangle.getCorner<Rect::TopRight>(), Point2D<Real>(4, 6));
  EXPECT_TRUE(rectangle.contains({1, 2}));
  EXPECT_TRUE(rectangle.contains({4, 6}));
  EXPECT_TRUE(rectangle.contains({2, 3}));
  EXPECT_FALSE(rectangle.contains({0, 3}));
  EXPECT_FALSE(rectangle.contains({2, 7}));
}

namespace Rodin::Tests::Unit
{
  template <class Scalar>
  class EuclideanSpatialStorage : public ::testing::Test
  {};

  using EuclideanScalars = ::testing::Types<Float, Real>;
  TYPED_TEST_SUITE(EuclideanSpatialStorage, EuclideanScalars);

  /// @brief Coordinates and directions use fixed-capacity spatial storage with value semantics.
  TYPED_TEST(EuclideanSpatialStorage, CoordinatesAndValueSemantics)
  {
    using Scalar = TypeParam;
    using Point = Point2D<Scalar>;
    using Vector = Math::SpatialVector<Scalar>;
    static_assert(
      std::is_same_v<std::remove_cvref_t<decltype(std::declval<Point>().asVector())>,
        Vector>);
    static_assert(
      std::is_same_v<decltype(std::declval<LineSegment2D<Scalar>>().direction()),
        Vector>);
    static_assert(Vector::Data::RowsAtCompileTime == 3);
    Point origin;
    EXPECT_EQ(origin.asVector().size(), 2);
    EXPECT_EQ(origin.x(), Scalar(0));
    EXPECT_EQ(origin.y(), Scalar(0));
    Point point(Scalar(3), Scalar(4));
    EXPECT_EQ(point.asVector().size(), 2);
    EXPECT_EQ(point.norm(), Scalar(5));
    Point copied(point);
    copied.asVector()[0] = Scalar(8);
    EXPECT_EQ(point.x(), Scalar(3));
    Point moved(std::move(copied));
    EXPECT_EQ(moved.x(), Scalar(8));
    Point assigned;
    assigned = point;
    assigned.y() = Scalar(9);
    EXPECT_EQ(point.y(), Scalar(4));
    assigned = Vector{Scalar(-2), Scalar(7)};
    EXPECT_EQ(assigned, Point(Scalar(-2), Scalar(7)));
    try
    {
      assigned = Vector{Scalar(1)};
      FAIL() << "Invalid coordinate dimension was accepted";
    }
    catch (const Alert::Exception& error)
    {
      const std::string message = error.what();
      EXPECT_NE(message.find("In member function requireTwoCoordinates of class "),
        std::string::npos);
      EXPECT_NE(message.find("Point2D"), std::string::npos);
      EXPECT_NE(
        message.find("Point2D requires two coordinates; received 1"), std::string::npos);
    }
    EXPECT_EQ(assigned, Point(Scalar(-2), Scalar(7)));
    EXPECT_THROW((Point(Vector{})), Alert::Exception);
    EXPECT_THROW((Point(Vector{Scalar(1), Scalar(2), Scalar(3)})), Alert::Exception);
    EXPECT_THROW((Point({Scalar(1)})), Alert::Exception);
  }

  /// @brief Spatial arithmetic preserves affine combinations and segment parameterization.
  TYPED_TEST(EuclideanSpatialStorage, AffineArithmeticAndInterpolation)
  {
    using Scalar = TypeParam;
    using Point = Point2D<Scalar>;
    using Vector = Math::SpatialVector<Scalar>;
    constexpr Real tolerance = std::is_same_v<Scalar, Float> ? 1e-5 : 1e-12;
    const Point first(Scalar(1), Scalar(2)), last(Scalar(4), Scalar(6));
    static_assert(std::is_same_v<decltype(last - first), Vector>);
    const Vector displacement = last - first;
    EXPECT_EQ(displacement.size(), 2);
    EXPECT_EQ(displacement[0], Scalar(3));
    EXPECT_EQ(displacement[1], Scalar(4));
    EXPECT_EQ(first + displacement, last);
    EXPECT_EQ(displacement + first, last);
    EXPECT_EQ(last - displacement, first);
    EXPECT_EQ((Scalar(2) * first) / Scalar(2), first);
    const LineSegment2D<Scalar> segment(first, last);
    EXPECT_NEAR(segment.direction()[0], Real(0.6), tolerance);
    EXPECT_NEAR(segment.direction()[1], Real(0.8), tolerance);
    EXPECT_NEAR(segment.direction().norm(), 1, tolerance);
    EXPECT_NEAR(segment(Scalar(0.25)).x(), Real(1.75), tolerance);
    EXPECT_NEAR(segment(Scalar(0.25)).y(), 3, tolerance);
    EXPECT_EQ(segment(Scalar(0)), first);
    EXPECT_EQ(segment(Scalar(1)), last);
    EXPECT_LT((segment.reverse().direction() + segment.direction()).norm(), tolerance);
    const Point center = barycenter(Point(Scalar(0), Scalar(0)),
      Point(Scalar(2), Scalar(0)), Point(Scalar(1), Scalar(2)));
    EXPECT_NEAR(center.x(), 1, tolerance);
    EXPECT_NEAR(center.y(), Real(2) / Real(3), tolerance);
  }

  /// @brief Projection and intersection results satisfy independent geometric equations.
  TYPED_TEST(EuclideanSpatialStorage, ProjectionIntersectionsAndCorners)
  {
    using Scalar = TypeParam;
    using Point = Point2D<Scalar>;
    constexpr Real tolerance = std::is_same_v<Scalar, Float> ? 1e-4 : 1e-12;
    const Line2D<Scalar> line(Scalar(3), Scalar(4), Scalar(5));
    const Point source(Scalar(2), Scalar(3));
    const auto projection = source.connect(line);
    ASSERT_TRUE(projection);
    EXPECT_NEAR(projection->end().x(), Real(11) / Real(25), tolerance);
    EXPECT_NEAR(projection->end().y(), Real(23) / Real(25), tolerance);
    EXPECT_NEAR(line(projection->end()), 0, tolerance);
    EXPECT_NEAR(projection->length(), Real(13) / Real(5), tolerance);
    Circle<Scalar> circle(Point(Scalar(1), Scalar(-2)), Scalar(5));
    const auto closest = Point(Scalar(7), Scalar(6)).connect(circle);
    ASSERT_TRUE(closest);
    EXPECT_LT((closest->end() - Point(Scalar(4), Scalar(2))).norm(), tolerance);
    EXPECT_NEAR(closest->length(), 5, tolerance);
    using Pair = std::pair<Point, Point>;
    const auto secant =
      Line2D<Scalar>(Scalar(3), Scalar(4), Scalar(-5)).intersect(circle);
    ASSERT_TRUE(std::holds_alternative<Pair>(secant));
    const auto& endpoints = std::get<Pair>(secant);
    EXPECT_NEAR((endpoints.first - endpoints.second).norm(), 10, tolerance);
    for (const auto& point : {endpoints.first, endpoints.second})
    {
      EXPECT_NEAR(circle(point), 0, tolerance);
      EXPECT_NEAR(Scalar(3) * point.x() + Scalar(4) * point.y(), -5, tolerance);
    }
    const auto tangent =
      Line2D<Scalar>(Scalar(-6), Scalar(-8), Scalar(-40)).intersect(circle);
    ASSERT_TRUE(std::holds_alternative<Point>(tangent));
    EXPECT_LT((std::get<Point>(tangent) - Point(Scalar(4), Scalar(2))).norm(), tolerance);
    using Rect = Rectangle<Scalar>;
    const Rect rectangle(Point(Scalar(-2), Scalar(1)), Point(Scalar(3), Scalar(5)));
    EXPECT_EQ(
      rectangle.template getCorner<Rect::BottomRight>(), Point(Scalar(3), Scalar(1)));
    EXPECT_EQ(
      rectangle.template getCorner<Rect::TopLeft>(), Point(Scalar(-2), Scalar(5)));
    EXPECT_EQ(rectangle.area(), Scalar(20));
    EXPECT_TRUE(rectangle.contains(Point(Scalar(-2), Scalar(5))));
    EXPECT_FALSE(rectangle.contains(Point(Scalar(4), Scalar(5))));
  }
}
