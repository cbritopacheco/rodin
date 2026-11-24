/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <limits>
#include <gtest/gtest.h>
#include <Rodin/Geometry/Ordering.h>

using namespace Rodin;
using namespace Rodin::Geometry;

TEST(Geometry_OrderingMorton, KnownKeysAndStableQuantizationTies)
{
  LocalMesh::Builder builder;
  builder.initialize(2).nodes(6);
  builder.vertex({1, 0}).vertex({0, 1}).vertex({0, 0}).vertex({1, 1});
  builder.vertex({0.1, 0.1}).vertex({0, 0});
  const auto mesh = builder.finalize();
  Ordering::Morton ordering;
  EXPECT_FALSE(ordering.getBitsPerCoordinate());
  ordering.setBitsPerCoordinate(1);
  EXPECT_EQ(ordering.getBitsPerCoordinate(), Optional<size_t>(1));
  EXPECT_EQ(ordering(mesh), (IndexVector{2, 4, 5, 1, 0, 3}));
  ordering.setBitsPerCoordinate(std::nullopt);
  EXPECT_FALSE(ordering.getBitsPerCoordinate());
  EXPECT_EQ(ordering(mesh), (IndexVector{2, 5, 4, 1, 0, 3}));
}

TEST(Geometry_OrderingMorton, PrecisionValidationAndFiniteCoordinates)
{
  Ordering::Morton ordering;
  EXPECT_THROW(ordering.setBitsPerCoordinate(0), Alert::Exception);
  EXPECT_THROW(ordering.setBitsPerCoordinate(53), Alert::Exception);
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  ordering.setBitsPerCoordinate(32);
  EXPECT_THROW(ordering(mesh), Alert::Exception);
  ordering.setBitsPerCoordinate(31);
  EXPECT_EQ(ordering(mesh).size(), mesh.getVertexCount());
  LocalMesh::Builder builder;
  builder.initialize(1).nodes(1).vertex({std::numeric_limits<Real>::quiet_NaN()});
  EXPECT_THROW(ordering(builder.finalize()), Alert::Exception);
}

TEST(Geometry_OrderingMorton, EmptyConstantAndExtremeCoordinates)
{
  LocalMesh::Builder empty;
  empty.initialize(2).nodes(0);
  EXPECT_TRUE(Ordering::Morton{}(empty.finalize()).empty());
  LocalMesh::Builder builder;
  const Real largest = std::numeric_limits<Real>::max();
  builder.initialize(2).nodes(4);
  builder.vertex({largest, 7}).vertex({-largest, 7}).vertex({0, 7}).vertex({0, 7});
  EXPECT_EQ(Ordering::Morton{}(builder.finalize()), (IndexVector{1, 2, 3, 0}));
}

TEST(Geometry_OrderingMorton, MaximumOneDimensionalPrecision)
{
  LocalMesh::Builder builder;
  builder.initialize(1).nodes(3).vertex({1}).vertex({0}).vertex({0.5});
  Ordering::Morton ordering;
  ordering.setBitsPerCoordinate(52);
  EXPECT_EQ(ordering(builder.finalize()), (IndexVector{1, 2, 0}));
}

TEST(Geometry_OrderingHilbert, KnownTwoDimensionalCurve)
{
  LocalMesh::Builder builder;
  builder.initialize(2).nodes(16);
  for (size_t x = 0; x < 4; ++x)
  {
    for (size_t y = 0; y < 4; ++y)
      builder.vertex({Real(x), Real(y)});
  }
  Ordering::Hilbert ordering;
  ordering.setBitsPerCoordinate(2);
  EXPECT_EQ(ordering(builder.finalize()),
    (IndexVector{0, 4, 5, 1, 2, 3, 7, 6, 10, 11, 15, 14, 13, 9, 8, 12}));
}

TEST(Geometry_OrderingHilbert, CompleteGridsHaveFaceAdjacentSuccessors)
{
  for (size_t dim = 1; dim <= 3; ++dim)
  {
    const size_t count = size_t(1) << (2 * dim);
    LocalMesh::Builder builder;
    builder.initialize(dim).nodes(count);
    for (Index v = 0; v < count; ++v)
    {
      Math::SpatialPoint point(dim);
      for (size_t d = 0; d < dim; ++d)
        point(d) = Real((v >> (2 * d)) & 3);
      builder.vertex(point);
    }
    const auto mesh = builder.finalize();
    Ordering::Hilbert ordering;
    ordering.setBitsPerCoordinate(2);
    const auto order = ordering(mesh);
    auto sorted = order;
    std::sort(sorted.begin(), sorted.end());
    ASSERT_EQ(sorted.size(), count);
    for (Index v = 0; v < count; ++v)
      EXPECT_EQ(sorted[v], v);
    for (size_t i = 1; i < order.size(); ++i)
    {
      const auto delta =
        mesh.getVertexCoordinates(order[i]) - mesh.getVertexCoordinates(order[i - 1]);
      Real distance = 0;
      for (size_t d = 0; d < dim; ++d)
        distance += std::abs(delta(d));
      EXPECT_EQ(distance, 1) << "dimension " << dim << ", step " << i;
    }
  }
}

TEST(Geometry_OrderingHilbert, PrecisionEmptyAndCoincidentVertices)
{
  Ordering::Hilbert ordering;
  EXPECT_FALSE(ordering.getBitsPerCoordinate());
  EXPECT_THROW(ordering.setBitsPerCoordinate(0), Alert::Exception);
  EXPECT_THROW(ordering.setBitsPerCoordinate(53), Alert::Exception);
  LocalMesh::Builder empty;
  empty.initialize(3).nodes(0);
  const auto emptyMesh = empty.finalize();
  EXPECT_TRUE(ordering(emptyMesh).empty());
  ordering.setBitsPerCoordinate(22);
  EXPECT_THROW(ordering(emptyMesh), Alert::Exception);
  ordering.setBitsPerCoordinate(21);
  LocalMesh::Builder builder;
  builder.initialize(3).nodes(3).vertex({2, 2, 2}).vertex({2, 2, 2}).vertex({2, 2, 2});
  EXPECT_EQ(ordering(builder.finalize()), (IndexVector{0, 1, 2}));
  ordering.setBitsPerCoordinate(std::nullopt);
  EXPECT_FALSE(ordering.getBitsPerCoordinate());
}

TEST(Geometry_OrderingRCM, ScrambledPathAndIsolatedVertex)
{
  LocalMesh::Builder builder;
  builder.initialize(1).nodes(6);
  for (size_t i = 0; i < 6; ++i)
    builder.vertex({Real(i)});
  builder.polytope(Polytope::Type::Segment, {0, 3});
  builder.polytope(Polytope::Type::Segment, {3, 1});
  builder.polytope(Polytope::Type::Segment, {1, 4});
  builder.polytope(Polytope::Type::Segment, {4, 2});
  const auto mesh = builder.finalize();
  const auto order = Ordering::RCM{}(mesh);
  EXPECT_EQ(order, (IndexVector{5, 0, 3, 1, 4, 2}));
  EXPECT_EQ(order, Ordering::RCM{}(mesh));
  IndexVector inverse(order.size());
  for (size_t i = 0; i < order.size(); ++i)
    inverse[order[i]] = i;
  for (Index c = 0; c < mesh.getCellCount(); ++c)
  {
    const auto& vertices = mesh.getConnectivity().getPolytope(1, c);
    EXPECT_EQ(std::max(inverse[vertices[0]], inverse[vertices[1]]) -
        std::min(inverse[vertices[0]], inverse[vertices[1]]),
      1);
  }
}

TEST(Geometry_OrderingRCM, NonSimplexCellsInduceCliques)
{
  LocalMesh::Builder builder;
  builder.initialize(2).nodes(4);
  builder.vertex({0, 0}).vertex({1, 0}).vertex({1, 1}).vertex({0, 1});
  builder.polytope(Polytope::Type::Quadrilateral, {0, 1, 2, 3});
  const auto mesh = builder.finalize();
  EXPECT_EQ(mesh.getConnectivity().getCount(1), 0);
  // All degrees are three. An edge-only graph would produce a different order.
  EXPECT_EQ(Ordering::RCM{}(mesh), (IndexVector{3, 2, 0, 1}));
  EXPECT_EQ(mesh.getConnectivity().getCount(1), 0);
}

TEST(Geometry_OrderingRCM, DisconnectedComponentsAndEmptyMesh)
{
  LocalMesh::Builder empty;
  empty.initialize(1).nodes(0);
  EXPECT_TRUE(Ordering::RCM{}(empty.finalize()).empty());
  LocalMesh::Builder builder;
  builder.initialize(1).nodes(4);
  for (size_t i = 0; i < 4; ++i)
    builder.vertex({Real(i)});
  builder.polytope(Polytope::Type::Segment, {0, 2});
  builder.polytope(Polytope::Type::Segment, {1, 3});
  EXPECT_EQ(Ordering::RCM{}(builder.finalize()), (IndexVector{0, 2, 1, 3}));
}
