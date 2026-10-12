/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Analytic contracts for the Euclidean primitive templates. */
#include <gtest/gtest.h>
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
