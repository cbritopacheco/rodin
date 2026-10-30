/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <utility>
#include <gtest/gtest.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/Shard.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Builder operations reject uninitialized and incompatible construction modes.
  TEST(Rodin_Geometry_ShardBuilderContract, RejectsIncompatibleModes)
  {
    const Math::SpatialPoint first{0.0}, second{2.0};
    const IndexArray vertices{{0, 1}};
    Shard::Builder uninitialized;
    EXPECT_THROW(uninitialized.vertex(11, first, Shard::State::Owned), std::exception);
    EXPECT_THROW(uninitialized.polytope(
                   1, 41, Polytope::Type::Segment, vertices, Shard::State::Owned),
      std::exception);
    EXPECT_THROW(uninitialized.include({0, 0}, Shard::State::Owned), std::exception);

    Shard::Builder direct;
    direct.initialize(1, 1);
    EXPECT_THROW(direct.include({0, 0}, Shard::State::Owned), std::exception);
    direct.vertex(11, first, Shard::State::Owned);
    direct.vertex(19, second, Shard::State::Shared);
    direct.polytope(1, 41, Polytope::Type::Segment, vertices, Shard::State::Owned);
    const auto directShard = direct.finalize();
    EXPECT_EQ(directShard.getVertexCount(), 2u);
    EXPECT_EQ(directShard.getCellCount(), 1u);

    auto parent = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
    Shard::Builder fromParent;
    fromParent.initialize(parent);
    EXPECT_THROW(fromParent.vertex(11, first, Shard::State::Owned), std::exception);
    EXPECT_THROW(
      fromParent.polytope(1, 41, Polytope::Type::Segment, vertices, Shard::State::Owned),
      std::exception);
    fromParent.include({0, 0}, Shard::State::Owned);
    fromParent.include({0, 1}, Shard::State::Owned);
    fromParent.include({1, 0}, Shard::State::Owned);
    const auto parentShard = fromParent.finalize();
    EXPECT_EQ(parentShard.getVertexCount(), 2u);
    EXPECT_EQ(parentShard.getCellCount(), 1u);

    Shard::Builder incomplete;
    incomplete.initialize(parent);
    incomplete.include({0, 0}, Shard::State::Owned);
    EXPECT_THROW(incomplete.include({1, 0}, Shard::State::Owned), std::exception);
  }

  /// @brief Rejected connectivity preserves pending direct entities and insertion indices.
  TEST(Rodin_Geometry_ShardBuilderContract, InvalidDirectConnectivityPreservesEntities)
  {
    Shard::Builder builder;
    builder.initialize(1, 1);
    ASSERT_EQ(builder.vertex(11, Math::SpatialPoint{0.0}, Shard::State::Owned), 0u);
    ASSERT_EQ(builder.vertex(19, Math::SpatialPoint{2.0}, Shard::State::Shared), 1u);
    const IndexArray invalid{{0, 5}}, valid{{0, 1}}, reversed{{1, 0}};
    EXPECT_THROW(
      builder.polytope(0, 11, Polytope::Type::Point, valid, Shard::State::Owned),
      std::exception);
    EXPECT_THROW(
      builder.polytope(1, 41, Polytope::Type::Segment, invalid, Shard::State::Owned),
      std::exception);
    EXPECT_EQ(builder.getPolytopeCount(0), 2u);
    EXPECT_EQ(builder.getPolytopeCount(1), 0u);
    ASSERT_EQ(
      builder.polytope(1, 41, Polytope::Type::Segment, valid, Shard::State::Owned), 0u);
    builder.attribute({1, 0}, Attribute(9));
    ASSERT_EQ(
      builder.polytope(1, 41, Polytope::Type::Segment, reversed, Shard::State::Ghost),
      0u);
    EXPECT_EQ(builder.getPolytopeCount(1), 1u);
    auto shard = builder.finalize();
    EXPECT_EQ(shard.getVertexCount(), 2u);
    EXPECT_EQ(shard.getCellCount(), 1u);
    EXPECT_TRUE(shard.isOwned(1, 0));
    EXPECT_TRUE(shard.isShared(0, 1));
    EXPECT_EQ(shard.getAttribute(1, 0), Optional<Attribute>(9));
    EXPECT_EQ(shard.getPolytopeMap(0).left[0], 11u);
    EXPECT_EQ(shard.getPolytopeMap(0).left[1], 19u);
    EXPECT_EQ(shard.getPolytopeMap(1).left[0], 41u);
    EXPECT_DOUBLE_EQ(shard.getVertexCoordinates(0)[0], 0);
    EXPECT_DOUBLE_EQ(shard.getVertexCoordinates(1)[0], 2);
    const auto& incidence = shard.getConnectivity().getIncidence(1, 0).at(0);
    ASSERT_EQ(incidence.size(), 2u);
    EXPECT_EQ(incidence[0], 0u);
    EXPECT_EQ(incidence[1], 1u);
  }

  /// @brief Parent overrides and owner/halo views survive finalization without parent mutation.
  TEST(Rodin_Geometry_ShardBuilderContract, ParentOverridesAndMetadataViews)
  {
    auto parent = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
    parent.setAttribute({1, 0}, 3);
    parent.setAttribute({0, 0}, 12);
    Shard::Builder builder;
    builder.initialize(parent);
    builder.include({0, 0}, Shard::State::Owned);
    builder.include({0, 1}, Shard::State::Shared);
    builder.include({1, 0}, Shard::State::Owned);
    builder.attribute({1, 0}, Attribute(17));
    builder.attribute({0, 0}, Optional<Attribute>{});
    builder.setOwner(0, 1, 5);
    builder.halo(0, 0, 5).halo(0, 0, 7).halo(0, 0, 7);
    EXPECT_EQ(builder.getOwner(0).at(1), 5u);
    EXPECT_EQ(builder.getHalo(0).at(0), (IndexSet{5, 7}));
    const auto& view = std::as_const(builder);
    EXPECT_EQ(view.getOwner(0).at(1), 5u);
    EXPECT_EQ(view.getHalo(0).at(0), (IndexSet{5, 7}));
    const auto shard = builder.finalize();
    EXPECT_EQ(shard.getOwner(0).at(1), 5u);
    EXPECT_EQ(shard.getHalo(0).at(0), (IndexSet{5, 7}));
    EXPECT_EQ(shard.getAttribute(1, 0), Optional<Attribute>(17));
    EXPECT_FALSE(shard.getAttribute(0, 0).has_value());
    EXPECT_EQ(parent.getAttribute(1, 0), Optional<Attribute>(3));
    EXPECT_EQ(parent.getAttribute(0, 0), Optional<Attribute>(12));
  }
}
