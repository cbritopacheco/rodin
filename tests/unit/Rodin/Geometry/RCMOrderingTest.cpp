/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <Rodin/Geometry/RCMOrdering.h>

using namespace Rodin;
using namespace Rodin::Geometry;

TEST(Geometry_RCMOrdering, ScrambledPathAndIsolatedVertex)
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
  const auto order = RCMOrdering{}(mesh);
  EXPECT_EQ(order, (IndexVector{5, 0, 3, 1, 4, 2}));
  EXPECT_EQ(order, RCMOrdering{}(mesh));
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

TEST(Geometry_RCMOrdering, NonSimplexCellsInduceCliques)
{
  LocalMesh::Builder builder;
  builder.initialize(2).nodes(4);
  builder.vertex({0, 0}).vertex({1, 0}).vertex({1, 1}).vertex({0, 1});
  builder.polytope(Polytope::Type::Quadrilateral, {0, 1, 2, 3});
  const auto mesh = builder.finalize();
  EXPECT_EQ(mesh.getConnectivity().getCount(1), 0);
  // All degrees are three. An edge-only graph would produce a different order.
  EXPECT_EQ(RCMOrdering{}(mesh), (IndexVector{3, 2, 0, 1}));
  EXPECT_EQ(mesh.getConnectivity().getCount(1), 0);
}

TEST(Geometry_RCMOrdering, DisconnectedComponentsAndEmptyMesh)
{
  LocalMesh::Builder empty;
  empty.initialize(1).nodes(0);
  EXPECT_TRUE(RCMOrdering{}(empty.finalize()).empty());
  LocalMesh::Builder builder;
  builder.initialize(1).nodes(4);
  for (size_t i = 0; i < 4; ++i)
    builder.vertex({Real(i)});
  builder.polytope(Polytope::Type::Segment, {0, 2});
  builder.polytope(Polytope::Type::Segment, {1, 3});
  EXPECT_EQ(RCMOrdering{}(builder.finalize()), (IndexVector{0, 2, 1, 3}));
}
