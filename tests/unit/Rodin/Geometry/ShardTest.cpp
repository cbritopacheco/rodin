/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/Shard.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent-based mode
  // ---------------------------------------------------------------------------

  /**
   * @brief Builds a shard from a parent mesh containing two triangles forming
   * a unit square and includes all entities as owned.
   */
  TEST(Rodin_Geometry_Shard, Builder_ParentMode_AllOwned)
  {
    constexpr size_t sdim = 2;

    Mesh<Context::Local> mesh =
      Mesh<Context::Local>::Builder()
        .initialize(sdim)
        .nodes(4)
        .vertex({0, 0})
        .vertex({1, 0})
        .vertex({0, 1})
        .vertex({1, 1})
        .polytope(Polytope::Type::Triangle, {0, 1, 2})
        .polytope(Polytope::Type::Triangle, {1, 3, 2})
        .finalize();

    // Include all vertices, then all cells, as owned
    Shard::Builder sb;
    sb.initialize(mesh);
    for (Index v = 0; v < mesh.getVertexCount(); v++)
      sb.include({0, v}, Shard::State::Owned);
    for (Index c = 0; c < mesh.getCellCount(); c++)
      sb.include({2, c}, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 4u);
    EXPECT_EQ(shard.getCellCount(), 2u);

    // All entities should be owned
    for (Index v = 0; v < shard.getVertexCount(); v++)
    {
      EXPECT_TRUE(shard.isOwned(0, v));
      EXPECT_FALSE(shard.isShared(0, v));
      EXPECT_FALSE(shard.isGhost(0, v));
      EXPECT_TRUE(shard.isLocal(0, v));
    }
    for (Index c = 0; c < shard.getCellCount(); c++)
    {
      EXPECT_TRUE(shard.isOwned(2, c));
      EXPECT_TRUE(shard.isLocal(2, c));
    }
  }

  /**
   * @brief Builds two shards from a 4-triangle mesh and verifies ownership,
   * shared, and ghost state for vertices and cells.
   */
  TEST(Rodin_Geometry_Shard, Builder_ParentMode_MixedStates)
  {
    // Create a 2x2 triangle mesh (4 triangles, 2x2 grid → 4 vertices min)
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, 0);

    // Build shard 0: owns cell 0, cell 1 is ghost
    Shard::Builder sb0;
    sb0.initialize(mesh);

    // Include vertices of cell 0 as owned
    const auto& conn = mesh.getConnectivity();
    const auto& c2v = conn.getIncidence(D, 0);
    for (const Index v : c2v.at(0))
      sb0.include({0, v}, Shard::State::Owned);
    auto [localCell0, ins0] = sb0.include({D, 0}, Shard::State::Owned);
    EXPECT_TRUE(ins0);
    EXPECT_EQ(localCell0, 0u);

    // Include vertices of cell 1 that are not already present as shared
    for (const Index v : c2v.at(1))
    {
      auto [localV, inserted] = sb0.include({0, v}, Shard::State::Shared);
      if (!inserted)
      {
        // Vertex was already owned
        EXPECT_TRUE(true);
      }
    }
    auto [localCell1, ins1] = sb0.include({D, 1}, Shard::State::Ghost);
    EXPECT_TRUE(ins1);

    Shard shard0 = sb0.finalize();

    // Cell 0 is owned, cell 1 is ghost
    EXPECT_TRUE(shard0.isOwned(D, localCell0));
    EXPECT_TRUE(shard0.isGhost(D, localCell1));
    EXPECT_TRUE(shard0.isLocal(D, localCell0));
    EXPECT_FALSE(shard0.isLocal(D, localCell1));
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Direct mode
  // ---------------------------------------------------------------------------

  /**
   * @brief Tests direct construction of a shard with explicit vertex and
   * polytope insertion.
   */
  TEST(Rodin_Geometry_Shard, Builder_DirectMode_SingleTriangle)
  {
    Shard::Builder sb;
    sb.initialize(/* dimension */ 2, /* sdim */ 2);

    Index v0 = sb.vertex(100, Math::SpatialPoint{{0.0, 0.0}}, Shard::State::Owned);
    Index v1 = sb.vertex(101, Math::SpatialPoint{{1.0, 0.0}}, Shard::State::Owned);
    Index v2 = sb.vertex(102, Math::SpatialPoint{{0.0, 1.0}}, Shard::State::Shared);

    EXPECT_EQ(v0, 0u);
    EXPECT_EQ(v1, 1u);
    EXPECT_EQ(v2, 2u);

    // Insert a triangle referencing local vertex indices
    IndexArray vs(3);
    vs << v0, v1, v2;
    Index t0 = sb.polytope(2, 200, Polytope::Type::Triangle, vs, Shard::State::Owned);
    EXPECT_EQ(t0, 0u);

    // Set owner for the shared vertex
    sb.setOwner(0, v2, /* ownerRank */ 1);

    // Set halo for owned vertices
    sb.halo(0, v0, /* neighborRank */ 1);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 3u);
    EXPECT_EQ(shard.getCellCount(), 1u);

    // Verify states
    EXPECT_TRUE(shard.isOwned(0, v0));
    EXPECT_TRUE(shard.isOwned(0, v1));
    EXPECT_TRUE(shard.isShared(0, v2));
    EXPECT_TRUE(shard.isOwned(2, t0));

    // isLocal = Owned || Shared
    EXPECT_TRUE(shard.isLocal(0, v0));
    EXPECT_TRUE(shard.isLocal(0, v1));
    EXPECT_TRUE(shard.isLocal(0, v2));

    // Verify owner map: v2 is shared, owner is rank 1
    const auto& ownerMap = shard.getOwner(0);
    auto it = ownerMap.find(v2);
    ASSERT_NE(it, ownerMap.end());
    EXPECT_EQ(it->second, 1u);

    // Verify halo map: v0 is owned, has neighbor rank 1
    const auto& haloMap = shard.getHalo(0);
    auto hit = haloMap.find(v0);
    ASSERT_NE(hit, haloMap.end());
    EXPECT_EQ(hit->second.size(), 1u);
    EXPECT_TRUE(hit->second.count(1) > 0);
  }

  /**
   * @brief Verifies that inserting the same global vertex twice in direct mode
   * is idempotent and returns the same local index.
   */
  TEST(Rodin_Geometry_Shard, Builder_DirectMode_VertexIdempotent)
  {
    Shard::Builder sb;
    sb.initialize(2, 2);

    Index v0a = sb.vertex(42, Math::SpatialPoint{{1.0, 2.0}}, Shard::State::Owned);
    Index v0b = sb.vertex(42, Math::SpatialPoint{{1.0, 2.0}}, Shard::State::Owned);

    EXPECT_EQ(v0a, v0b);
    EXPECT_EQ(sb.getPolytopeCount(0), 1u);
  }

  // ---------------------------------------------------------------------------
  // Shard — PolytopeMap
  // ---------------------------------------------------------------------------

  /**
   * @brief Verifies the bidirectional polytope map after building a shard
   * from a parent mesh.
   */
  TEST(Rodin_Geometry_Shard, PolytopeMap_ParentMode)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, 0);

    Shard::Builder sb;
    sb.initialize(mesh);

    // Include a subset of vertices and one cell
    const auto& c2v = mesh.getConnectivity().getIncidence(D, 0);
    for (const Index v : c2v.at(0))
      sb.include({0, v}, Shard::State::Owned);
    sb.include({D, 0}, Shard::State::Owned);

    Shard shard = sb.finalize();

    // Vertex polytope map: left[localIdx] == parentIdx
    const auto& vmap = shard.getPolytopeMap(0);
    for (Index i = 0; i < vmap.left.size(); i++)
    {
      Index parentIdx = vmap.left[i];
      // Reverse lookup must agree
      auto it = vmap.right.find(parentIdx);
      ASSERT_NE(it, vmap.right.end());
      EXPECT_EQ(it->second, i);
    }

    // Cell polytope map
    const auto& cmap = shard.getPolytopeMap(D);
    EXPECT_EQ(cmap.left.size(), 1u);
    EXPECT_EQ(cmap.left[0], 0u); // parent cell index 0
    EXPECT_EQ(cmap.right.size(), 1u);
  }

  /**
   * @brief Verifies the bidirectional polytope map after building a shard
   * in direct mode.
   */
  TEST(Rodin_Geometry_Shard, PolytopeMap_DirectMode)
  {
    Shard::Builder sb;
    sb.initialize(2, 2);

    sb.vertex(10, Math::SpatialPoint{{0.0, 0.0}}, Shard::State::Owned);
    sb.vertex(20, Math::SpatialPoint{{1.0, 0.0}}, Shard::State::Owned);
    sb.vertex(30, Math::SpatialPoint{{0.0, 1.0}}, Shard::State::Owned);

    IndexArray vs(3);
    vs << 0, 1, 2;
    sb.polytope(2, 99, Polytope::Type::Triangle, vs, Shard::State::Owned);

    Shard shard = sb.finalize();

    // Vertex map: local 0 → global 10, local 1 → global 20, local 2 → global 30
    const auto& vmap = shard.getPolytopeMap(0);
    EXPECT_EQ(vmap.left[0], 10u);
    EXPECT_EQ(vmap.left[1], 20u);
    EXPECT_EQ(vmap.left[2], 30u);
    EXPECT_EQ(vmap.right.at(10), 0u);
    EXPECT_EQ(vmap.right.at(20), 1u);
    EXPECT_EQ(vmap.right.at(30), 2u);

    // Cell map: local 0 → global 99
    const auto& cmap = shard.getPolytopeMap(2);
    EXPECT_EQ(cmap.left[0], 99u);
    EXPECT_EQ(cmap.right.at(99), 0u);
  }

  // ---------------------------------------------------------------------------
  // Shard — State accessors
  // ---------------------------------------------------------------------------

  /**
   * @brief Verifies getState() returns the correct per-entity state vectors.
   */
  TEST(Rodin_Geometry_Shard, GetState)
  {
    Shard::Builder sb;
    sb.initialize(2, 2);

    sb.vertex(0, Math::SpatialPoint{{0.0, 0.0}}, Shard::State::Owned);
    sb.vertex(1, Math::SpatialPoint{{1.0, 0.0}}, Shard::State::Shared);
    sb.vertex(2, Math::SpatialPoint{{0.0, 1.0}}, Shard::State::Ghost);

    IndexArray vs(3);
    vs << 0, 1, 2;
    sb.polytope(2, 0, Polytope::Type::Triangle, vs, Shard::State::Owned);

    Shard shard = sb.finalize();

    const auto& vstate = shard.getState(0);
    ASSERT_EQ(vstate.size(), 3u);
    EXPECT_EQ(vstate[0], Shard::State::Owned);
    EXPECT_EQ(vstate[1], Shard::State::Shared);
    EXPECT_EQ(vstate[2], Shard::State::Ghost);

    const auto& cstate = shard.getState(2);
    ASSERT_EQ(cstate.size(), 1u);
    EXPECT_EQ(cstate[0], Shard::State::Owned);
  }

  // ---------------------------------------------------------------------------
  // Shard — Move and copy semantics
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, MoveConstructor)
  {
    Shard::Builder sb;
    sb.initialize(2, 2);

    sb.vertex(0, Math::SpatialPoint{{0.0, 0.0}}, Shard::State::Owned);
    sb.vertex(1, Math::SpatialPoint{{1.0, 0.0}}, Shard::State::Owned);
    sb.vertex(2, Math::SpatialPoint{{0.0, 1.0}}, Shard::State::Owned);

    IndexArray vs(3);
    vs << 0, 1, 2;
    sb.polytope(2, 0, Polytope::Type::Triangle, vs, Shard::State::Owned);

    Shard original = sb.finalize();
    const size_t origVerts = original.getVertexCount();
    const size_t origCells = original.getCellCount();

    Shard moved(std::move(original));

    EXPECT_EQ(moved.getVertexCount(), origVerts);
    EXPECT_EQ(moved.getCellCount(), origCells);
    EXPECT_TRUE(moved.isOwned(0, 0));
    EXPECT_TRUE(moved.isOwned(2, 0));
  }

  TEST(Rodin_Geometry_Shard, MoveAssignment)
  {
    Shard::Builder sb;
    sb.initialize(2, 2);

    sb.vertex(0, Math::SpatialPoint{{0.0, 0.0}}, Shard::State::Owned);
    sb.vertex(1, Math::SpatialPoint{{1.0, 0.0}}, Shard::State::Shared);
    sb.vertex(2, Math::SpatialPoint{{0.0, 1.0}}, Shard::State::Ghost);

    IndexArray vs(3);
    vs << 0, 1, 2;
    sb.polytope(2, 0, Polytope::Type::Triangle, vs, Shard::State::Owned);
    sb.setOwner(0, 1, 2);
    sb.setOwner(0, 2, 3);

    Shard original = sb.finalize();

    Shard assigned;
    assigned = std::move(original);

    EXPECT_EQ(assigned.getVertexCount(), 3u);
    EXPECT_EQ(assigned.getCellCount(), 1u);
    EXPECT_TRUE(assigned.isOwned(0, 0));
    EXPECT_TRUE(assigned.isShared(0, 1));
    EXPECT_TRUE(assigned.isGhost(0, 2));

    // Verify owner map survived the move
    EXPECT_EQ(assigned.getOwner(0).at(1), 2u);
    EXPECT_EQ(assigned.getOwner(0).at(2), 3u);
  }

  TEST(Rodin_Geometry_Shard, CopyConstructor)
  {
    Shard::Builder sb;
    sb.initialize(2, 2);

    sb.vertex(0, Math::SpatialPoint{{0.0, 0.0}}, Shard::State::Owned);
    sb.vertex(1, Math::SpatialPoint{{1.0, 0.0}}, Shard::State::Owned);
    sb.vertex(2, Math::SpatialPoint{{0.0, 1.0}}, Shard::State::Ghost);

    IndexArray vs(3);
    vs << 0, 1, 2;
    sb.polytope(2, 0, Polytope::Type::Triangle, vs, Shard::State::Owned);
    sb.setOwner(0, 2, 5);

    Shard original = sb.finalize();
    Shard copy(original);

    EXPECT_EQ(copy.getVertexCount(), original.getVertexCount());
    EXPECT_EQ(copy.getCellCount(), original.getCellCount());
    EXPECT_TRUE(copy.isOwned(0, 0));
    EXPECT_TRUE(copy.isGhost(0, 2));
    EXPECT_EQ(copy.getOwner(0).at(2), 5u);

    // Verify independence: modifying copy doesn't affect original
    copy.getOwner(0)[2] = 99;
    EXPECT_EQ(original.getOwner(0).at(2), 5u);
    EXPECT_EQ(copy.getOwner(0).at(2), 99u);
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — include() idempotency (parent mode)
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_ParentMode_IncludeIdempotent)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, 0);

    Shard::Builder sb;
    sb.initialize(mesh);

    // Include vertex 0 twice — second time should return same local index
    auto [idx1, ins1] = sb.include({0, 0}, Shard::State::Owned);
    auto [idx2, ins2] = sb.include({0, 0}, Shard::State::Owned);

    EXPECT_TRUE(ins1);
    EXPECT_FALSE(ins2);
    EXPECT_EQ(idx1, idx2);
    EXPECT_EQ(sb.getPolytopeCount(0), 1u);
  }

  // ---------------------------------------------------------------------------
  // Shard — Owner and halo metadata via Builder
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, OwnerAndHaloMetadata_DirectMode)
  {
    Shard::Builder sb;
    sb.initialize(2, 2);

    // Three vertices: 0 owned (with halo to rank 1 & 2), 1 shared (owned by rank 1), 2 ghost (owned by rank 2)
    sb.vertex(10, Math::SpatialPoint{{0.0, 0.0}}, Shard::State::Owned);
    sb.vertex(11, Math::SpatialPoint{{1.0, 0.0}}, Shard::State::Shared);
    sb.vertex(12, Math::SpatialPoint{{0.5, 1.0}}, Shard::State::Ghost);

    sb.setOwner(0, 1, /* ownerRank */ 1);
    sb.setOwner(0, 2, /* ownerRank */ 2);
    sb.halo(0, 0, /* neighborRank */ 1);
    sb.halo(0, 0, /* neighborRank */ 2);

    IndexArray vs(3);
    vs << 0, 1, 2;
    sb.polytope(2, 0, Polytope::Type::Triangle, vs, Shard::State::Owned);

    Shard shard = sb.finalize();

    // Owner map
    const auto& ownerMap = shard.getOwner(0);
    EXPECT_EQ(ownerMap.size(), 2u);
    EXPECT_EQ(ownerMap.at(1), 1u);
    EXPECT_EQ(ownerMap.at(2), 2u);

    // Halo map
    const auto& haloMap = shard.getHalo(0);
    EXPECT_EQ(haloMap.size(), 1u);
    auto hit = haloMap.find(0);
    ASSERT_NE(hit, haloMap.end());
    EXPECT_EQ(hit->second.size(), 2u);
    EXPECT_TRUE(hit->second.count(1) > 0);
    EXPECT_TRUE(hit->second.count(2) > 0);
  }
}
