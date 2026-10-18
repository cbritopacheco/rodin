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
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Deg.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Geometry/Euclidean.h"

using namespace Rodin;
using namespace Rodin::Geometry::Euclidean;

static_assert(HasIntersectMethod<Line<3>, Line<3>>::value);
static_assert(!HasIntersectMethod<Line<3>, int>::value);
static_assert(HasDistanceMethod<Line<4>, Point<4>>::value);
static_assert(!HasDistanceMethod<Line<4>, int>::value);
static_assert(HasConnectMethod<Point<3>, Line<3>>::value);
static_assert(!HasConnectMethod<Point<3>, int>::value);

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
  class EuclideanFixedStorage : public ::testing::Test
  {};

  using EuclideanScalars = ::testing::Types<Float, Real>;
  TYPED_TEST_SUITE(EuclideanFixedStorage, EuclideanScalars);

  /// @brief Coordinates and directions use fixed-capacity spatial storage with value semantics.
  TYPED_TEST(EuclideanFixedStorage, CoordinatesAndValueSemantics)
  {
    using Scalar = TypeParam;
    using Point = Point2D<Scalar>;
    using Vector = Math::FixedSizeVector<Scalar, 2>;
    static_assert(std::is_same_v<Point2D<Scalar>, Geometry::Euclidean::Point<2, Scalar>>);
    static_assert(std::is_same_v<LineSegment2D<Scalar>, LineSegment<2, Scalar>>);
    static_assert(std::is_same_v<Line2D<Scalar>, Line<2, Scalar>>);
    static_assert(std::is_same_v<Circle<Scalar>, Sphere<2, Scalar>>);
    static_assert(std::is_same_v<Rectangle<Scalar>, Box<2, Scalar>>);
    static_assert(
      std::is_same_v<std::remove_cvref_t<decltype(std::declval<Point>().asVector())>,
        Vector>);
    static_assert(
      std::is_same_v<decltype(std::declval<LineSegment2D<Scalar>>().direction()),
        Vector>);
    static_assert(Vector::RowsAtCompileTime == 2);
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
      assigned = Math::Vector<Scalar>::Ones(1);
      FAIL() << "Invalid coordinate dimension was accepted";
    }
    catch (const Alert::Exception& error)
    {
      const std::string message = error.what();
      EXPECT_NE(message.find("In member function require of class "), std::string::npos);
      EXPECT_NE(message.find("Point<"), std::string::npos);
      EXPECT_NE(
        message.find("Point requires 2 coordinates; received 1"), std::string::npos);
    }
    EXPECT_EQ(assigned, Point(Scalar(-2), Scalar(7)));
    EXPECT_THROW((Point(Math::Vector<Scalar>::Zero(0))), Alert::Exception);
    EXPECT_THROW((Point(Math::Vector<Scalar>::Ones(3))), Alert::Exception);
    EXPECT_THROW((Point({Scalar(1)})), Alert::Exception);
  }

  /// @brief Spatial arithmetic preserves affine combinations and segment parameterization.
  TYPED_TEST(EuclideanFixedStorage, AffineArithmeticAndInterpolation)
  {
    using Scalar = TypeParam;
    using Point = Point2D<Scalar>;
    using Vector = Math::FixedSizeVector<Scalar, 2>;
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
  TYPED_TEST(EuclideanFixedStorage, ProjectionIntersectionsAndCorners)
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

namespace Rodin::Tests::Unit
{
  template <size_t N, class S>
  struct EuclideanConfiguration
  {
    static constexpr size_t Dimension = N;
    using Scalar = S;
  };
  template <class Configuration>
  class EuclideanDimensions : public ::testing::Test
  {};
  using EuclideanConfigurations =
    ::testing::Types<EuclideanConfiguration<1, Float>, EuclideanConfiguration<2, Float>,
      EuclideanConfiguration<3, Float>, EuclideanConfiguration<4, Float>,
      EuclideanConfiguration<6, Float>, EuclideanConfiguration<1, Real>,
      EuclideanConfiguration<2, Real>, EuclideanConfiguration<3, Real>,
      EuclideanConfiguration<4, Real>, EuclideanConfiguration<6, Real>>;
  TYPED_TEST_SUITE(EuclideanDimensions, EuclideanConfigurations);

  TYPED_TEST(EuclideanDimensions, CoordinatesAndAffineOperations)
  {
    using S = typename TypeParam::Scalar;
    constexpr size_t N = TypeParam::Dimension;
    using P = Geometry::Euclidean::Point<N, S>;
    using V = typename P::Vector;
    static_assert(P::getDimension() == N);
    EXPECT_EQ(P::getDimension(), N);
    static_assert(V::RowsAtCompileTime == N);
    static_assert(V::SizeAtCompileTime == N);
    V coordinates = V::Zero(), offset = V::Zero();
    S squaredNorm = 0;
    for (size_t i = 0; i < N; ++i)
    {
      coordinates[i] = S(i + 1);
      offset[i] = S(2 * (i + 1));
      squaredNorm += S((i + 1) * (i + 1));
    }
    const P point(coordinates);
    const P variadic = []<size_t... I>(std::index_sequence<I...>) {
      return P(S(I + 1)...);
    }(std::make_index_sequence<N>{});
    const P listed = []<size_t... I>(std::index_sequence<I...>) {
      return P({S(I + 1)...});
    }(std::make_index_sequence<N>{});
    EXPECT_EQ(variadic, point);
    EXPECT_EQ(listed, point);
    const P copy(point);
    P moving(point);
    const P moved(std::move(moving));
    P named;
    named.x() = S(1);
    EXPECT_EQ(std::as_const(named).x(), S(1));
    if constexpr (N >= 2)
    {
      named.y() = S(2);
      EXPECT_EQ(std::as_const(named).y(), S(2));
    }
    if constexpr (N >= 3)
    {
      named.z() = S(3);
      EXPECT_EQ(std::as_const(named).z(), S(3));
    }
    static_assert((requires(P p) { p.y(); }) == (N >= 2));
    static_assert((requires(P p) { p.z(); }) == (N >= 3));
    if constexpr (N == 1)
    {
      EXPECT_EQ(P(S(1)), P({S(1)}));
    }
    if constexpr (N == 2)
    {
      EXPECT_EQ(P(S(1), S(2)), P({S(1), S(2)}));
    }
    if constexpr (N == 3)
    {
      EXPECT_EQ(P(S(1), S(2), S(3)), P({S(1), S(2), S(3)}));
    }
    EXPECT_EQ(point, copy);
    EXPECT_EQ(point, moved);
    const auto sum = point + offset;
    EXPECT_EQ(offset + point, sum);
    const auto middle = barycenter(point, P(offset));
    for (size_t i = 0; i < N; ++i)
    {
      EXPECT_EQ(point[i], S(i + 1));
      EXPECT_EQ(sum[i], S(3 * (i + 1)));
      EXPECT_EQ((sum - offset)[i], point[i]);
      EXPECT_EQ((point * S(2))[i], offset[i]);
      EXPECT_EQ((S(2) * point / S(2))[i], point[i]);
      EXPECT_EQ(middle[i], S(1.5) * S(i + 1));
      EXPECT_EQ((P(offset) - point)[i], point[i]);
    }
    EXPECT_NEAR(point.norm(), std::sqrt(squaredNorm),
      S(8) * std::numeric_limits<S>::epsilon() * squaredNorm);
    P alias(point);
    alias = alias.asVector().reverse();
    for (size_t i = 0; i < N; ++i)
      EXPECT_EQ(alias[i], S(N - i));
    P changed(point);
    changed[N - 1] += S(1);
    EXPECT_NE(point, changed);
    Math::Vector<S> wrong = Math::Vector<S>::Zero(N - 1);
    EXPECT_THROW(changed = wrong, Alert::Exception);
    EXPECT_EQ(changed[N - 1], S(N + 1));
    EXPECT_THROW((void)P(wrong), Alert::Exception);
    EXPECT_THROW(point + wrong, Alert::Exception);
    if constexpr (N >= 2)
    {
      const Math::Matrix<S> row = Math::Matrix<S>::Ones(1, N);
      EXPECT_THROW((void)P(row), Alert::Exception);
      EXPECT_THROW(changed = row, Alert::Exception);
      EXPECT_EQ(changed[N - 1], S(N + 1));
      EXPECT_THROW((Line<N, S>(point, row)), Alert::Exception);
      EXPECT_THROW((Hyperplane<N, S>(point, row)), Alert::Exception);
      EXPECT_THROW((Sphere<N, S>().pointAt(row)), Alert::Exception);
    }
    EXPECT_THROW(P({}), Alert::Exception);
    EXPECT_THROW(P(std::initializer_list<S>{S(1), S(2), S(3), S(4), S(5), S(6), S(7)}),
      Alert::Exception);
  }

  TYPED_TEST(EuclideanDimensions, SegmentInterpolationAndIntersections)
  {
    using S = typename TypeParam::Scalar;
    constexpr size_t N = TypeParam::Dimension;
    using P = Geometry::Euclidean::Point<N, S>;
    using Segment = LineSegment<N, S>;
    P origin, end;
    end[N - 1] = S(8);
    Segment segment(origin, end);
    EXPECT_EQ(segment.length(), S(8));
    EXPECT_EQ(segment.direction()[N - 1], S(1));
    const P middle = segment(S(0.25));
    EXPECT_EQ(middle[N - 1], S(2));
    EXPECT_EQ(segment.reverse()(S(0.75)), middle);
    EXPECT_EQ(segment.project(origin - S(2) * segment.direction()), origin);
    EXPECT_EQ(segment.project(end + S(2) * segment.direction()), end);
    EXPECT_EQ(segment.distance(end + S(2) * segment.direction()), S(2));
    P shortEnd, adjacentEnd;
    shortEnd[N - 1] = S(0.1);
    adjacentEnd[N - 1] = S(0.4);
    const Segment shortSegment(origin, shortEnd), adjacent(shortEnd, adjacentEnd);
    ASSERT_TRUE(shortSegment.intersect(adjacent));
    EXPECT_EQ(*shortSegment.intersect(adjacent), shortEnd);
    EXPECT_FALSE(segment.intersect(segment));
    P further = end;
    further[N - 1] = S(12);
    const auto touch = segment.intersect(Segment(end, further));
    ASSERT_TRUE(touch);
    EXPECT_EQ(*touch, end);
    if constexpr (N >= 2)
    {
      P low = middle, high = middle;
      low[0] = S(-2);
      high[0] = S(2);
      ASSERT_TRUE(segment.intersect(Segment(low, high)));
      EXPECT_EQ(*segment.intersect(Segment(low, high)), middle);
      low[N - 1] = high[N - 1] = S(10);
      EXPECT_FALSE(segment.intersect(Segment(low, high)));
      if constexpr (N >= 3)
      {
        low[1] = high[1] = S(1);
        low[N - 1] = high[N - 1] = S(2);
        EXPECT_FALSE(segment.intersect(Segment(low, high)));
      }
    }
    EXPECT_THROW(Segment(origin, origin), Alert::Exception);
    EXPECT_THROW(segment.setStart(end), Alert::Exception);
    EXPECT_THROW(segment.setEnd(origin), Alert::Exception);
    EXPECT_EQ(segment.start(), origin);
    EXPECT_EQ(segment.end(), end);
    segment.setStart(origin + segment.direction()).setEnd(end + segment.direction());
    EXPECT_EQ(segment.start()[N - 1], S(1));
    EXPECT_EQ(segment.end()[N - 1], S(9));
    EXPECT_FALSE(origin.connect(origin));
    ASSERT_TRUE(origin.connect(end));
    EXPECT_EQ(origin.connect(end)->length(), S(8));
  }

  TYPED_TEST(EuclideanDimensions, LineProjectionAndSkewDistance)
  {
    using S = typename TypeParam::Scalar;
    constexpr size_t N = TypeParam::Dimension;
    using P = Geometry::Euclidean::Point<N, S>;
    using L = Line<N, S>;
    using V = typename P::Vector;
    P origin, end;
    end[N - 1] = S(4);
    L line(origin, end);
    EXPECT_EQ(line.getDimension(), N);
    EXPECT_EQ(line.origin(), origin);
    EXPECT_EQ(line.direction()[N - 1], S(1));
    EXPECT_EQ(line.pointAt(S(0.5))[N - 1], S(2));
    EXPECT_EQ(line.project(end), end);
    EXPECT_EQ(line.distance(end), S(0));
    EXPECT_FALSE(end.connect(line));
    const L scaled(origin, V::Constant(S(0.1))), sameLine(origin, V::Constant(S(0.3)));
    EXPECT_FALSE(scaled.intersect(sameLine));
    for (const S scale :
      {std::numeric_limits<S>::min(), std::numeric_limits<S>::max() / S(4)})
    {
      V scaledDirection = V::Zero();
      scaledDirection[N - 1] = scale;
      const L rescaled(origin, scaledDirection);
      EXPECT_EQ(rescaled.direction()[N - 1], S(1));
      EXPECT_EQ(rescaled.project(end), end);
      EXPECT_EQ(rescaled.distance(end), S(0));
      EXPECT_FALSE(rescaled.intersect(line));
      const Hyperplane<N, S> plane(origin, scaledDirection);
      EXPECT_EQ(plane.project(end), origin);
      const Sphere<N, S> sphere(origin, S(2));
      EXPECT_EQ(sphere.pointAt(scaledDirection)[N - 1], S(2));
    }
    EXPECT_FALSE(line.intersect(line));
    EXPECT_EQ(line.closestPoint(line), origin);
    EXPECT_EQ(line.distance(line), S(0));
    EXPECT_THROW(L(origin, V::Zero()), Alert::Exception);
    EXPECT_THROW(L(origin, Math::Vector<S>::Zero(N - 1)), Alert::Exception);
    if constexpr (N == 2)
    {
      const P query(S(1), S(1));
      const L nearlyHorizontal(std::numeric_limits<S>::min(), S(1), S(1));
      const L nearlyVertical(S(1), std::numeric_limits<S>::min(), S(1));
      EXPECT_EQ(nearlyHorizontal.project(query), query);
      EXPECT_EQ(nearlyVertical.project(query), query);
      EXPECT_FALSE(query.connect(nearlyHorizontal));
      EXPECT_FALSE(query.connect(nearlyVertical));
      const L coefficients(std::initializer_list<S>{S(3), S(4), S(5)});
      EXPECT_EQ(coefficients.a(), S(3));
      EXPECT_EQ(coefficients.b(), S(4));
      EXPECT_EQ(coefficients.c(), S(5));
      EXPECT_EQ(*coefficients.slope(), S(-0.75));
      EXPECT_EQ(*coefficients.xIntercept(), S(5) / S(3));
      EXPECT_EQ(*coefficients.yIntercept(), S(1.25));
      EXPECT_EQ(*coefficients(S(1)), S(0.5));
      EXPECT_EQ(coefficients(P(S(3), S(-1))), S(0));
      EXPECT_FALSE(coefficients.isHorizontal());
      EXPECT_FALSE(coefficients.isVertical());
      const L tangent = Sphere<N, S>().tangent(Math::Rad(0));
      EXPECT_EQ(tangent(P(S(1), S(0))), S(0));
    }
    if constexpr (N >= 2)
    {
      P p = end;
      p[0] = S(3);
      EXPECT_EQ(line.project(p), end);
      EXPECT_EQ(line.distance(p), S(3));
      const auto connection = p.connect(line);
      ASSERT_TRUE(connection);
      EXPECT_EQ(connection->end(), end);
      P q = p;
      q[N - 1] = S(8);
      const L parallel(p, q);
      EXPECT_EQ(line.distance(parallel), S(3));
      EXPECT_FALSE(line.intersect(parallel));
      const L crossing(end, p);
      ASSERT_TRUE(line.intersect(crossing));
      EXPECT_EQ(*line.intersect(crossing), end);
      EXPECT_EQ(line.distance(crossing), S(0));
      if constexpr (N >= 3)
      {
        P a = end, b = p;
        a[1] = b[1] = S(2);
        const L skew(a, b);
        EXPECT_EQ(line.distance(skew), S(2));
        EXPECT_FALSE(line.intersect(skew));
      }
    }
  }

  TYPED_TEST(EuclideanDimensions, SphereSurfaceTangentsAndLineIntersections)
  {
    using S = typename TypeParam::Scalar;
    constexpr size_t N = TypeParam::Dimension;
    using P = Geometry::Euclidean::Point<N, S>;
    using V = typename P::Vector;
    using Ball = Sphere<N, S>;
    using Pair = std::pair<P, P>;
    P center;
    for (size_t i = 0; i < N; ++i)
      center[i] = S(i + 1);
    Ball sphere(center, S(2));
    EXPECT_EQ(sphere.getDimension(), N);
    V direction = V::Zero();
    direction[N - 1] = S(4);
    const P surface = sphere.pointAt(direction);
    EXPECT_EQ(surface[N - 1], center[N - 1] + S(2));
    EXPECT_EQ(sphere(surface), S(0));
    EXPECT_EQ(sphere(center), S(-4));
    EXPECT_EQ(sphere.distance(center), S(2));
    const auto tangent = sphere.tangent(direction);
    EXPECT_EQ(tangent.getDimension(), N);
    EXPECT_EQ(tangent(surface), S(0));
    EXPECT_EQ(tangent.distance(center), S(2));
    EXPECT_EQ(tangent.project(center), surface);
    EXPECT_EQ(tangent.normal()[N - 1], S(4));
    EXPECT_EQ(tangent.origin(), surface);
    const Line<N, S> diameter(center, direction);
    EXPECT_FALSE(diameter.connect(sphere));
    const P inside = center + S(0.25) * direction;
    ASSERT_TRUE(inside.connect(sphere));
    EXPECT_EQ(inside.connect(sphere)->end(), surface);
    const auto result = diameter.intersect(sphere);
    ASSERT_TRUE(std::holds_alternative<Pair>(result));
    const auto points = std::get<Pair>(result);
    EXPECT_EQ(sphere(points.first), S(0));
    EXPECT_EQ(sphere(points.second), S(0));
    EXPECT_EQ((points.first - points.second).norm(), S(4));
    EXPECT_FALSE(center.connect(sphere));
    EXPECT_FALSE(surface.connect(sphere));
    const P outside = center + direction;
    ASSERT_TRUE(outside.connect(sphere));
    EXPECT_EQ(outside.connect(sphere)->end(), surface);
    EXPECT_THROW(sphere.setRadius(S(0)), Alert::Exception);
    EXPECT_EQ(sphere.radius(), S(2));
    EXPECT_THROW(sphere.pointAt(V::Zero()), Alert::Exception);
    EXPECT_THROW(sphere.pointAt(Math::Vector<S>::Zero(N - 1)), Alert::Exception);
    EXPECT_THROW((Hyperplane<N, S>(center, V::Zero())), Alert::Exception);
    if constexpr (N >= 2)
    {
      P tangentOrigin = center, separatedOrigin = center;
      tangentOrigin[0] += S(2);
      separatedOrigin[0] += S(3);
      const Line<N, S> tangentLine(tangentOrigin, direction),
        separated(separatedOrigin, direction);
      const auto touch = tangentLine.intersect(sphere);
      ASSERT_TRUE(std::holds_alternative<P>(touch));
      EXPECT_EQ(std::get<P>(touch), tangentOrigin);
      EXPECT_TRUE(std::holds_alternative<std::nullopt_t>(separated.intersect(sphere)));
      ASSERT_TRUE(separated.connect(sphere));
      EXPECT_EQ(separated.connect(sphere)->length(), S(1));
      EXPECT_FALSE(tangentLine.connect(sphere));
    }
    sphere.setCenter(P()).setRadius(S(1));
    EXPECT_EQ(sphere.center(), P());
    EXPECT_EQ(sphere.radius(), S(1));
  }

  TYPED_TEST(EuclideanDimensions, BoxContainmentProjectionAndIntersection)
  {
    using S = typename TypeParam::Scalar;
    constexpr size_t N = TypeParam::Dimension;
    using P = Geometry::Euclidean::Point<N, S>;
    using B = Box<N, S>;
    P lower, upper;
    S measure = S(1);
    for (size_t i = 0; i < N; ++i)
    {
      lower[i] = S(i + 1);
      upper[i] = lower[i] + S(2);
      measure *= S(2);
    }
    const B box(lower, upper);
    EXPECT_EQ(box.getDimension(), N);
    if constexpr (N == 2)
    {
      EXPECT_EQ(box.width(), S(2));
      EXPECT_EQ(box.height(), S(2));
      EXPECT_EQ(box.area(), S(4));
    }
    EXPECT_EQ(box.measure(), measure);
    EXPECT_EQ(box.lower(), lower);
    EXPECT_EQ(box.upper(), upper);
    EXPECT_TRUE(box.contains(lower));
    EXPECT_TRUE(box.contains(upper));
    EXPECT_FALSE(box.contains(P()));
    std::array<bool, N> mask{};
    for (size_t i = 0; i < N; ++i)
    {
      mask[i] = true;
      EXPECT_EQ(box.getCorner(mask)[i], upper[i]);
      EXPECT_EQ(box.extent(i), S(2));
    }
    EXPECT_EQ(box.getCorner(mask), upper);
    P outside = upper;
    outside[N - 1] += S(3);
    EXPECT_EQ(box.project(outside), upper);
    EXPECT_EQ(box.distance(outside), S(3));
    EXPECT_EQ(box.distance(lower), S(0));
    const auto overlap = box.intersect(B(lower, lower + upper.asVector()));
    ASSERT_TRUE(overlap);
    EXPECT_EQ(overlap->measure(), measure);
    const auto touch = box.intersect(B(upper, outside));
    ASSERT_TRUE(touch);
    EXPECT_EQ(touch->measure(), S(0));
    P separated = upper;
    separated[N - 1] += S(1);
    EXPECT_FALSE(box.intersect(B(separated, outside)));
    EXPECT_THROW(B(upper, lower), Alert::Exception);
    EXPECT_THROW(box.extent(N), Alert::Exception);
  }
}

namespace Rodin::Tests::Unit
{
  TYPED_TEST(EuclideanDimensions, ManufacturedObliqueGeometry)
  {
    using S = typename TypeParam::Scalar;
    constexpr size_t N = TypeParam::Dimension;
    using P = Geometry::Euclidean::Point<N, S>;
    using V = typename P::Vector;
    using L = Line<N, S>;
    P origin;
    V u = V::Zero(), normal = V::Zero();
    S squared = 0;
    for (size_t i = 0; i < N; ++i)
    {
      origin[i] = S(3 * (i + 1));
      u[i] = S(i + 1);
      squared += S((i + 1) * (i + 1));
    }
    const S epsilon = std::numeric_limits<S>::epsilon();
    const S tolerance = S(64) * epsilon;
    const L line(origin, u);
    const auto along = origin + S(2) * u;
    if constexpr (N >= 2)
    {
      normal[0] = S(2);
      normal[1] = S(-1);
      const auto sample = along + S(3) * normal;
      const auto projected = line.project(sample);
      for (size_t i = 0; i < N; ++i)
        EXPECT_NEAR(projected[i], origin[i] + S(2 * (i + 1)), tolerance);
      EXPECT_NEAR(line.distance(sample), S(3) * std::sqrt(S(5)), tolerance);
      const L crossing(along, normal);
      ASSERT_TRUE(line.intersect(crossing));
      EXPECT_LT((*line.intersect(crossing) - along).norm(), tolerance);
      if constexpr (N == 3)
      {
        const V transverse(S(3), S(6), S(-5));
        const L skew(along + transverse, normal);
        EXPECT_NEAR(line.distance(skew), std::sqrt(S(70)), tolerance);
        EXPECT_FALSE(line.intersect(skew));
      }
    }
    const S radius = S(2) * std::sqrt(squared);
    const Sphere<N, S> sphere(origin, radius);
    const auto surface = sphere.pointAt(u);
    for (size_t i = 0; i < N; ++i)
      EXPECT_NEAR(surface[i], origin[i] + S(2 * (i + 1)), tolerance);
    EXPECT_NEAR(sphere(surface), S(0), tolerance * squared);
    const auto tangent = sphere.tangent(u);
    EXPECT_NEAR(tangent.distance(origin), radius, tolerance);
    for (size_t i = 0; i < N; ++i)
      EXPECT_NEAR(tangent.project(origin)[i], surface[i], tolerance);
    const auto result = line.intersect(sphere);
    using Pair = std::pair<P, P>;
    ASSERT_TRUE(std::holds_alternative<Pair>(result));
    const auto pair = std::get<Pair>(result);
    for (size_t i = 0; i < N; ++i)
    {
      EXPECT_NEAR(pair.first[i], origin[i] + S(2 * (i + 1)), tolerance);
      EXPECT_NEAR(pair.second[i], origin[i] - S(2 * (i + 1)), tolerance);
    }
    const Box<N, S> box(origin, origin + S(2) * u);
    for (size_t bits = 0; bits < (size_t(1) << N); ++bits)
    {
      std::array<bool, N> mask{};
      for (size_t i = 0; i < N; ++i)
        mask[i] = (bits & (size_t(1) << i)) != 0;
      const auto corner = box.getCorner(mask);
      EXPECT_TRUE(box.contains(corner));
      for (size_t i = 0; i < N; ++i)
        EXPECT_EQ(corner[i], origin[i] + (mask[i] ? S(2 * (i + 1)) : S(0)));
    }
  }
}

TEST(Rodin_Geometry_Euclidean, InvalidInputsCarryMemberContext)
{
  const auto requireMember = [](auto operation, const std::string& member) {
    try
    {
      operation();
      FAIL() << "Expected an Alert exception";
    }
    catch (const Alert::Exception& error)
    {
      const std::string message = error.what();
      EXPECT_NE(
        message.find("In member function " + member + " of class "), std::string::npos)
        << message;
    }
  };
  requireMember([] { (void)Line<2>(Real(0), Real(0), Real(1)); }, "Line");
  requireMember(
    [] { (void)Line<2>(std::initializer_list<Real>{Real(1), Real(2)}); }, "Line");
  requireMember(
    [] { (void)Line<3>(Point<3>(), Math::FixedSizeVector<Real, 3>::Zero()); }, "Line");
  requireMember(
    [] {
      Sphere<3> sphere;
      sphere.setRadius(-1);
    },
    "setRadius");
  requireMember(
    [] {
      Sphere<3> sphere;
      (void)sphere.pointAt(Math::FixedSizeVector<Real, 3>::Zero());
    },
    "pointAt");
  requireMember(
    [] { (void)Hyperplane<3>(Point<3>(), Math::Vector<Real>::Zero(2)); }, "Hyperplane");
  requireMember([] { (void)Box<3>({1, 2, 3}, {1, 1, 4}); }, "Box");
  requireMember([] { (void)Box<3>({0, 0, 0}, {1, 1, 1}).extent(3); }, "extent");
}

// Instantiate every supported public member, including variants hidden by shared line coverage.
template class Rodin::Geometry::Euclidean::Point<1, Rodin::Float>;
template class Rodin::Geometry::Euclidean::LineSegment<1, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Line<1, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Sphere<1, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Hyperplane<1, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Box<1, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Point<1, Rodin::Real>;
template class Rodin::Geometry::Euclidean::LineSegment<1, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Line<1, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Sphere<1, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Hyperplane<1, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Box<1, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Point<2, Rodin::Float>;
template class Rodin::Geometry::Euclidean::LineSegment<2, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Line<2, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Sphere<2, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Hyperplane<2, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Box<2, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Point<2, Rodin::Real>;
template class Rodin::Geometry::Euclidean::LineSegment<2, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Line<2, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Sphere<2, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Hyperplane<2, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Box<2, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Point<3, Rodin::Float>;
template class Rodin::Geometry::Euclidean::LineSegment<3, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Line<3, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Sphere<3, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Hyperplane<3, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Box<3, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Point<3, Rodin::Real>;
template class Rodin::Geometry::Euclidean::LineSegment<3, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Line<3, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Sphere<3, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Hyperplane<3, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Box<3, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Point<4, Rodin::Float>;
template class Rodin::Geometry::Euclidean::LineSegment<4, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Line<4, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Sphere<4, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Hyperplane<4, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Box<4, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Point<4, Rodin::Real>;
template class Rodin::Geometry::Euclidean::LineSegment<4, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Line<4, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Sphere<4, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Hyperplane<4, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Box<4, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Point<6, Rodin::Float>;
template class Rodin::Geometry::Euclidean::LineSegment<6, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Line<6, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Sphere<6, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Hyperplane<6, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Box<6, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Point<6, Rodin::Real>;
template class Rodin::Geometry::Euclidean::LineSegment<6, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Line<6, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Sphere<6, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Hyperplane<6, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Box<6, Rodin::Real>;

namespace Rodin::Tests::Unit
{
  template <class Scalar>
  class EuclideanZeroDimensions : public ::testing::Test
  {};
  TYPED_TEST_SUITE(EuclideanZeroDimensions, EuclideanScalars);

  TYPED_TEST(EuclideanZeroDimensions, PointAndBoxCountingMeasure)
  {
    using S = TypeParam;
    using P = Geometry::Euclidean::Point<0, S>;
    using V = Math::FixedSizeVector<S, 0>;
    using B = Box<0, S>;
    static_assert(V::SizeAtCompileTime == 0);
    static_assert(!(requires(P p) { p.x(); }));
    const P origin;
    EXPECT_EQ(P::getDimension(), 0);
    EXPECT_EQ(origin.asVector().size(), 0);
    EXPECT_EQ(origin.norm(), S(0));
    EXPECT_EQ(origin, P({}));
    EXPECT_EQ(origin, P(V::Zero()));
    EXPECT_EQ((origin - origin).size(), 0);
    EXPECT_EQ(origin + origin, origin);
    EXPECT_EQ(origin + V::Zero(), origin);
    EXPECT_EQ(V::Zero() + origin, origin);
    EXPECT_EQ(origin - V::Zero(), origin);
    EXPECT_EQ(origin * S(2), origin);
    EXPECT_EQ(barycenter(origin, origin, origin), origin);
    EXPECT_FALSE(origin.connect(origin));
    const B box(origin, origin);
    EXPECT_EQ(box.getDimension(), 0);
    EXPECT_TRUE(box.contains(origin));
    EXPECT_EQ(box.measure(), S(1));
    EXPECT_EQ(box.getCorner(std::array<bool, 0>{}), origin);
    EXPECT_EQ(box.project(origin), origin);
    EXPECT_EQ(box.distance(origin), S(0));
    const auto intersection = box.intersect(box);
    ASSERT_TRUE(intersection);
    EXPECT_EQ(intersection->measure(), S(1));
    EXPECT_THROW(box.extent(0), Alert::Exception);
    EXPECT_THROW((void)P(Math::Vector<S>::Ones(1)), Alert::Exception);
  }

  /// @brief The unique zero-dimensional point and box still own distinct coordinate objects.
  TYPED_TEST(EuclideanZeroDimensions, StorageAndBoxCornersRemainOwned)
  {
    using S = TypeParam;
    using P = Geometry::Euclidean::Point<0, S>;
    using B = Box<0, S>;
    P point;
    auto& coordinates = point.asVector();
    EXPECT_EQ(coordinates.rows(), 0);
    EXPECT_EQ(coordinates.cols(), 1);
    EXPECT_EQ(coordinates.size(), 0);
    using V = Math::FixedSizeVector<S, 0>;
    const V base = V::Zero(), offset = V::Zero();
    const auto affine = base + S(2) * offset / S(3);
    const P constructed(affine);
    EXPECT_EQ(constructed, point);
    EXPECT_EQ(&(point = affine), &point);
    const P& view = point;
    EXPECT_EQ(&coordinates, &view.asVector());
    const B box(point, point);
    EXPECT_EQ(box.lower(), point);
    EXPECT_EQ(box.upper(), point);
    EXPECT_NE(&box.lower(), &point);
    EXPECT_NE(&box.upper(), &point);
    EXPECT_NE(&box.lower(), &box.upper());
    const B copied(box);
    EXPECT_EQ(copied.lower(), point);
    EXPECT_EQ(copied.upper(), point);
    EXPECT_NE(&copied.lower(), &box.lower());
    EXPECT_NE(&copied.upper(), &box.upper());
    EXPECT_NE(&copied.lower(), &copied.upper());
  }

  TYPED_TEST(EuclideanZeroDimensions, RejectsNonemptyPositiveDimensionalPrimitives)
  {
    using S = TypeParam;
    using P = Geometry::Euclidean::Point<0, S>;
    using V = Math::FixedSizeVector<S, 0>;
    using L = Line<0, S>;
    using Segment = LineSegment<0, S>;
    using Ball = Sphere<0, S>;
    using Plane = Hyperplane<0, S>;
    const P origin;
    EXPECT_THROW((void)L(origin, origin), Alert::Exception);
    EXPECT_THROW((void)Segment(origin, origin), Alert::Exception);
    EXPECT_THROW((void)Ball(), Alert::Exception);
    EXPECT_THROW((void)Plane(origin, V::Zero()), Alert::Exception);
  }
}

template class Rodin::Geometry::Euclidean::Point<0, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Point<0, Rodin::Real>;
template class Rodin::Geometry::Euclidean::Box<0, Rodin::Float>;
template class Rodin::Geometry::Euclidean::Box<0, Rodin::Real>;

/// @brief Concrete degree arithmetic remains usable by the planar sphere angle API.
TEST(Rodin_Geometry_Euclidean, DegreeArithmeticFeedsSphereParameterization)
{
  const Sphere<2> circle({2, 3}, 2);
  const Math::Deg start(30), increment(60);
  const auto quarterTurn = (start + increment).toRad();
  const auto point = circle(quarterTurn);
  EXPECT_NEAR(point[0], 2, 1e-12);
  EXPECT_NEAR(point[1], 5, 1e-12);
  const auto tangent = circle.tangent(quarterTurn);
  EXPECT_NEAR(tangent(point), 0, 1e-12);
  EXPECT_NEAR(tangent.distance(Point<2>(2, 7)), 2, 1e-12);
  const auto zero = ((start + increment) - Math::Deg(90)).toRad();
  EXPECT_NEAR(circle(zero)[0], 4, 1e-12);
  EXPECT_NEAR(circle(zero)[1], 3, 1e-12);
  const auto negative = (-(start + increment)).toRad();
  EXPECT_NEAR(circle(negative)[0], 2, 1e-12);
  EXPECT_NEAR(circle(negative)[1], 1, 1e-12);
}
