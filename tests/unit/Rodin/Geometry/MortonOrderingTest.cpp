/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <limits>
#include <gtest/gtest.h>
#include <Rodin/Geometry/MortonOrdering.h>

using namespace Rodin;
using namespace Rodin::Geometry;

TEST(Geometry_MortonOrdering, KnownKeysAndStableQuantizationTies)
{
  LocalMesh::Builder builder;
  builder.initialize(2).nodes(6);
  builder.vertex({1, 0}).vertex({0, 1}).vertex({0, 0}).vertex({1, 1});
  builder.vertex({0.1, 0.1}).vertex({0, 0});
  const auto mesh = builder.finalize();
  MortonOrdering ordering;
  EXPECT_FALSE(ordering.getBitsPerCoordinate());
  ordering.setBitsPerCoordinate(1);
  EXPECT_EQ(ordering.getBitsPerCoordinate(), Optional<size_t>(1));
  EXPECT_EQ(ordering(mesh), (IndexVector{2, 4, 5, 1, 0, 3}));
  ordering.setBitsPerCoordinate(std::nullopt);
  EXPECT_FALSE(ordering.getBitsPerCoordinate());
  EXPECT_EQ(ordering(mesh), (IndexVector{2, 5, 4, 1, 0, 3}));
}

TEST(Geometry_MortonOrdering, PrecisionValidationAndFiniteCoordinates)
{
  MortonOrdering ordering;
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

TEST(Geometry_MortonOrdering, EmptyConstantAndExtremeCoordinates)
{
  LocalMesh::Builder empty;
  empty.initialize(2).nodes(0);
  EXPECT_TRUE(MortonOrdering{}(empty.finalize()).empty());
  LocalMesh::Builder builder;
  const Real largest = std::numeric_limits<Real>::max();
  builder.initialize(2).nodes(4);
  builder.vertex({largest, 7}).vertex({-largest, 7}).vertex({0, 7}).vertex({0, 7});
  EXPECT_EQ(MortonOrdering{}(builder.finalize()), (IndexVector{1, 2, 3, 0}));
}

TEST(Geometry_MortonOrdering, MaximumOneDimensionalPrecision)
{
  LocalMesh::Builder builder;
  builder.initialize(1).nodes(3).vertex({1}).vertex({0}).vertex({0.5});
  MortonOrdering ordering;
  ordering.setBitsPerCoordinate(52);
  EXPECT_EQ(ordering(builder.finalize()), (IndexVector{1, 2, 0}));
}
