/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or copy at https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <Rodin/Location.h>
#include <Rodin/Variational.h>
#include <Rodin/Geometry/Shard.h>

using namespace Rodin;
using namespace Rodin::Geometry;

TEST(Location_AABBCandidates, PreservesSparseIndicesAndFiltersFallback)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
  Location::AABB<LocalMesh>::Candidates candidates(2);
  candidates[1] = {3, 1};
  Location::AABB locator(mesh, candidates);
  locator.setExhaustiveFallback(true).setProjectionPruning(true);
  for (Index i = 0; i < 4; ++i)
  {
    Point sample(*mesh.getCell(i), Polytope::Traits(Polytope::Type::Segment).getCentroid());
    auto hit = locator.locate(sample.getPhysicalCoordinates());
    ASSERT_EQ(hit.has_value(), i == 1 || i == 3);
    if (hit)
      EXPECT_EQ(hit->getPolytope().getIndex(), i);
  }
  locator.setTolerance(locator.getTolerance());
  EXPECT_TRUE(locator.locate(Point(*mesh.getCell(3),
    Polytope::Traits(Polytope::Type::Segment).getCentroid()).getPhysicalCoordinates()));
}

TEST(Location_AABBCandidates, EmptySelectionAndInvalidIndices)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  Location::AABB<LocalMesh>::Candidates candidates(2);
  Location::AABB locator(mesh, candidates);
  const auto x = Point(*mesh.getCell(),
    Polytope::Traits(Polytope::Type::Segment).getCentroid()).getPhysicalCoordinates();
  EXPECT_FALSE(locator.setExhaustiveFallback(true).locate(x));
  candidates[1] = {2};
  EXPECT_THROW((Location::AABB(mesh, candidates)), std::invalid_argument);
  candidates[1] = {0, 0};
  EXPECT_THROW((Location::AABB(mesh, candidates)), std::invalid_argument);
  candidates.resize(3);
  candidates[1].clear();
  candidates[2] = {0};
  EXPECT_THROW((Location::AABB(mesh, candidates)), std::invalid_argument);
}

TEST(Location_AABBCandidates, ParentShardPreservesCurvedTransformation)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {1});
  Variational::RealH1Element<3> element(Polytope::Type::Segment);
  PointCloud nodes(1, element.getCount());
  // Monotone cubic, preserving endpoints and visibly differing from the affine map.
  constexpr Real Curvature = Real(0.2);
  for (size_t a = 0; a < element.getCount(); ++a)
  {
    const Real t = element.getNode(a)[0];
    nodes(0, a) = t + Curvature * t * (1 - t) * (1 - 2 * t);
  }
  mesh.setPolytopeTransformation({1, 0},
    new ParametricTransformation<Variational::RealH1Element<3>>(std::move(nodes), element));
  Shard::Builder builder;
  builder.initialize(mesh);
  for (Index v = 0; v < mesh.getVertexCount(); ++v)
    builder.include({0, v}, Shard::State::Owned);
  builder.include({1, 0}, Shard::State::Owned);
  auto shard = builder.finalize();
  Math::SpatialPoint rc(1);
  rc[0] = Real(0.25);
  Point parent(*mesh.getCell(), rc), child(*shard.getCell(), rc);
  EXPECT_EQ(child.getPolytope().getTransformation().getOrder(), 3);
  EXPECT_LT((parent.getPhysicalCoordinates() - child.getPhysicalCoordinates()).norm(), 1e-12);
  EXPECT_LT((parent.getJacobian() - child.getJacobian()).norm(), 1e-12);
}
