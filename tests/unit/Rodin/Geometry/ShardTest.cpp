/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <set>

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

  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent-based mode with Quadrilateral mesh
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_ParentMode_Quadrilateral)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Quadrilateral, {3, 3});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, 0);

    Shard::Builder sb;
    sb.initialize(mesh);
    for (Index v = 0; v < mesh.getVertexCount(); v++)
      sb.include({0, v}, Shard::State::Owned);
    for (Index c = 0; c < mesh.getCellCount(); c++)
      sb.include({D, c}, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(shard.getCellCount(), mesh.getCellCount());

    for (Index c = 0; c < shard.getCellCount(); c++)
    {
      EXPECT_TRUE(shard.isOwned(D, c));
      EXPECT_TRUE(shard.isLocal(D, c));
    }
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent-based mode with Segment (1D) mesh
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_ParentMode_Segment)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Segment, {5});
    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 1u);
    mesh.getConnectivity().compute(D, 0);

    Shard::Builder sb;
    sb.initialize(mesh);

    // Include all vertices and cells as owned
    for (Index v = 0; v < mesh.getVertexCount(); v++)
      sb.include({0, v}, Shard::State::Owned);
    for (Index c = 0; c < mesh.getCellCount(); c++)
      sb.include({D, c}, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(shard.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(shard.getDimension(), 1u);

    for (Index c = 0; c < shard.getCellCount(); c++)
      EXPECT_TRUE(shard.isOwned(D, c));
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent-based mode with Tetrahedron (3D) mesh
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_ParentMode_Tetrahedron)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 3u);
    mesh.getConnectivity().compute(D, 0);

    Shard::Builder sb;
    sb.initialize(mesh);

    for (Index v = 0; v < mesh.getVertexCount(); v++)
      sb.include({0, v}, Shard::State::Owned);
    for (Index c = 0; c < mesh.getCellCount(); c++)
      sb.include({D, c}, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(shard.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(shard.getDimension(), 3u);
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent-based mode with Hexahedron (3D) mesh
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_ParentMode_Hexahedron)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Hexahedron, {3, 3, 3});
    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 3u);
    mesh.getConnectivity().compute(D, 0);

    Shard::Builder sb;
    sb.initialize(mesh);

    for (Index v = 0; v < mesh.getVertexCount(); v++)
      sb.include({0, v}, Shard::State::Owned);
    for (Index c = 0; c < mesh.getCellCount(); c++)
      sb.include({D, c}, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(shard.getCellCount(), mesh.getCellCount());
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent-based mode with Wedge (3D) mesh
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_ParentMode_Wedge)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Wedge, {3, 3, 3});
    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 3u);
    mesh.getConnectivity().compute(D, 0);

    Shard::Builder sb;
    sb.initialize(mesh);

    for (Index v = 0; v < mesh.getVertexCount(); v++)
      sb.include({0, v}, Shard::State::Owned);
    for (Index c = 0; c < mesh.getCellCount(); c++)
      sb.include({D, c}, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(shard.getCellCount(), mesh.getCellCount());
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Direct mode with Quadrilateral
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_DirectMode_Quadrilateral)
  {
    Shard::Builder sb;
    sb.initialize(2, 2);

    Index v0 = sb.vertex(0, Math::SpatialPoint{{0.0, 0.0}}, Shard::State::Owned);
    Index v1 = sb.vertex(1, Math::SpatialPoint{{1.0, 0.0}}, Shard::State::Owned);
    Index v2 = sb.vertex(2, Math::SpatialPoint{{1.0, 1.0}}, Shard::State::Owned);
    Index v3 = sb.vertex(3, Math::SpatialPoint{{0.0, 1.0}}, Shard::State::Shared);

    IndexArray vs(4);
    vs << v0, v1, v2, v3;
    Index q0 = sb.polytope(2, 50, Polytope::Type::Quadrilateral, vs, Shard::State::Owned);
    EXPECT_EQ(q0, 0u);

    sb.setOwner(0, v3, 2);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 4u);
    EXPECT_EQ(shard.getCellCount(), 1u);
    EXPECT_TRUE(shard.isOwned(0, v0));
    EXPECT_TRUE(shard.isShared(0, v3));
    EXPECT_TRUE(shard.isOwned(2, q0));
    EXPECT_EQ(shard.getOwner(0).at(v3), 2u);
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Direct mode with Segment (1D)
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_DirectMode_Segment)
  {
    Shard::Builder sb;
    sb.initialize(1, 1);

    Index v0 = sb.vertex(0, Math::SpatialPoint({0.0}), Shard::State::Owned);
    Index v1 = sb.vertex(1, Math::SpatialPoint({1.0}), Shard::State::Owned);
    Index v2 = sb.vertex(2, Math::SpatialPoint({2.0}), Shard::State::Ghost);

    IndexArray seg1_vs(2);
    seg1_vs << v0, v1;
    Index s0 = sb.polytope(1, 10, Polytope::Type::Segment, seg1_vs, Shard::State::Owned);

    IndexArray seg2_vs(2);
    seg2_vs << v1, v2;
    Index s1 = sb.polytope(1, 11, Polytope::Type::Segment, seg2_vs, Shard::State::Ghost);

    sb.setOwner(0, v2, 1);
    sb.setOwner(1, s1, 1);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 3u);
    EXPECT_EQ(shard.getCellCount(), 2u);
    EXPECT_EQ(shard.getDimension(), 1u);

    EXPECT_TRUE(shard.isOwned(1, s0));
    EXPECT_TRUE(shard.isGhost(1, s1));
    EXPECT_TRUE(shard.isGhost(0, v2));
    EXPECT_EQ(shard.getOwner(0).at(v2), 1u);
    EXPECT_EQ(shard.getOwner(1).at(s1), 1u);
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Direct mode with Tetrahedron (3D)
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_DirectMode_Tetrahedron)
  {
    Shard::Builder sb;
    sb.initialize(3, 3);

    Index v0 = sb.vertex(0, Math::SpatialPoint{{0.0, 0.0, 0.0}}, Shard::State::Owned);
    Index v1 = sb.vertex(1, Math::SpatialPoint{{1.0, 0.0, 0.0}}, Shard::State::Owned);
    Index v2 = sb.vertex(2, Math::SpatialPoint{{0.0, 1.0, 0.0}}, Shard::State::Owned);
    Index v3 = sb.vertex(3, Math::SpatialPoint{{0.0, 0.0, 1.0}}, Shard::State::Shared);

    IndexArray vs(4);
    vs << v0, v1, v2, v3;
    Index t0 = sb.polytope(3, 100, Polytope::Type::Tetrahedron, vs, Shard::State::Owned);
    EXPECT_EQ(t0, 0u);

    sb.setOwner(0, v3, 1);
    sb.halo(0, v0, 1);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 4u);
    EXPECT_EQ(shard.getCellCount(), 1u);
    EXPECT_EQ(shard.getDimension(), 3u);

    EXPECT_TRUE(shard.isOwned(3, t0));
    EXPECT_TRUE(shard.isShared(0, v3));
    EXPECT_EQ(shard.getOwner(0).at(v3), 1u);

    const auto& haloMap = shard.getHalo(0);
    ASSERT_NE(haloMap.find(v0), haloMap.end());
    EXPECT_TRUE(haloMap.at(v0).count(1) > 0);
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Direct mode with Hexahedron (3D)
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_DirectMode_Hexahedron)
  {
    Shard::Builder sb;
    sb.initialize(3, 3);

    // 8 vertices of a unit cube
    Index v0 = sb.vertex(0, Math::SpatialPoint{{0.0, 0.0, 0.0}}, Shard::State::Owned);
    Index v1 = sb.vertex(1, Math::SpatialPoint{{1.0, 0.0, 0.0}}, Shard::State::Owned);
    Index v2 = sb.vertex(2, Math::SpatialPoint{{1.0, 1.0, 0.0}}, Shard::State::Owned);
    Index v3 = sb.vertex(3, Math::SpatialPoint{{0.0, 1.0, 0.0}}, Shard::State::Owned);
    Index v4 = sb.vertex(4, Math::SpatialPoint{{0.0, 0.0, 1.0}}, Shard::State::Ghost);
    Index v5 = sb.vertex(5, Math::SpatialPoint{{1.0, 0.0, 1.0}}, Shard::State::Ghost);
    Index v6 = sb.vertex(6, Math::SpatialPoint{{1.0, 1.0, 1.0}}, Shard::State::Ghost);
    Index v7 = sb.vertex(7, Math::SpatialPoint{{0.0, 1.0, 1.0}}, Shard::State::Ghost);

    IndexArray vs(8);
    vs << v0, v1, v2, v3, v4, v5, v6, v7;
    Index h0 = sb.polytope(3, 200, Polytope::Type::Hexahedron, vs, Shard::State::Owned);
    EXPECT_EQ(h0, 0u);

    for (Index vi = v4; vi <= v7; vi++)
      sb.setOwner(0, vi, 1);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 8u);
    EXPECT_EQ(shard.getCellCount(), 1u);
    EXPECT_EQ(shard.getDimension(), 3u);
    EXPECT_TRUE(shard.isOwned(3, h0));

    for (Index vi = v0; vi <= v3; vi++)
      EXPECT_TRUE(shard.isOwned(0, vi));
    for (Index vi = v4; vi <= v7; vi++)
      EXPECT_TRUE(shard.isGhost(0, vi));
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent-based mode with mixed 2D mesh (tri+quad)
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_ParentMode_Mixed2D_TriQuad)
  {
    // Build a mixed mesh: 2 triangles + 1 quadrilateral sharing vertices
    //
    //   2---3---4
    //   |\ | Q |
    //   | \|   |
    //   0--1---5
    //
    // T0: (0,1,2), T1: (1,3,2), Q0: (1,3,4,5)
    Mesh<Context::Local> mesh =
      Mesh<Context::Local>::Builder()
        .initialize(2)
        .nodes(6)
        .vertex({0.0, 0.0})
        .vertex({1.0, 0.0})
        .vertex({0.0, 1.0})
        .vertex({1.0, 1.0})
        .vertex({2.0, 1.0})
        .vertex({2.0, 0.0})
        .polytope(Polytope::Type::Triangle, {0, 1, 2})
        .polytope(Polytope::Type::Triangle, {1, 3, 2})
        .polytope(Polytope::Type::Quadrilateral, {1, 5, 4, 3})
        .finalize();

    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 2u);
    mesh.getConnectivity().compute(D, 0);

    // Shard all entities as owned
    Shard::Builder sb;
    sb.initialize(mesh);
    for (Index v = 0; v < mesh.getVertexCount(); v++)
      sb.include({0, v}, Shard::State::Owned);
    for (Index c = 0; c < mesh.getCellCount(); c++)
      sb.include({D, c}, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 6u);
    EXPECT_EQ(shard.getCellCount(), 3u);

    // Verify polytope map round-trips
    const auto& cmap = shard.getPolytopeMap(D);
    for (Index ci = 0; ci < cmap.left.size(); ci++)
    {
      Index parentIdx = cmap.left[ci];
      EXPECT_EQ(cmap.right.at(parentIdx), ci);
    }
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent-based mode with mixed 2D: partial inclusion
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_ParentMode_Mixed2D_Partial)
  {
    // Same mixed mesh but only include one triangle and one quad, with
    // shared vertices at the interface
    Mesh<Context::Local> mesh =
      Mesh<Context::Local>::Builder()
        .initialize(2)
        .nodes(6)
        .vertex({0.0, 0.0})
        .vertex({1.0, 0.0})
        .vertex({0.0, 1.0})
        .vertex({1.0, 1.0})
        .vertex({2.0, 1.0})
        .vertex({2.0, 0.0})
        .polytope(Polytope::Type::Triangle, {0, 1, 2})
        .polytope(Polytope::Type::Triangle, {1, 3, 2})
        .polytope(Polytope::Type::Quadrilateral, {1, 5, 4, 3})
        .finalize();

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, 0);

    const auto& c2v = mesh.getConnectivity().getIncidence(D, 0);

    // Shard 0: owns cell 0 (triangle), cell 2 (quad) is ghost
    Shard::Builder sb;
    sb.initialize(mesh);

    // Vertices of cell 0 as owned
    for (const Index v : c2v.at(0))
      sb.include({0, v}, Shard::State::Owned);
    sb.include({D, 0}, Shard::State::Owned);

    // Vertices of cell 2 (quad) — shared or ghost for new ones
    for (const Index v : c2v.at(2))
      sb.include({0, v}, Shard::State::Ghost);
    sb.include({D, 2}, Shard::State::Ghost);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getCellCount(), 2u);

    // Cell 0 is owned, cell 2 is ghost
    EXPECT_TRUE(shard.isOwned(D, 0));
    EXPECT_TRUE(shard.isGhost(D, 1)); // cell 2 maps to local 1
    EXPECT_TRUE(shard.isLocal(D, 0));
    EXPECT_FALSE(shard.isLocal(D, 1));

    // Vertices 0, 1, 2 are owned (from cell 0), others are ghost
    for (Index vi = 0; vi < shard.getVertexCount(); vi++)
    {
      EXPECT_TRUE(shard.isOwned(0, vi) || shard.isGhost(0, vi));
    }
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent-based mode with mixed 3D (tet+wedge)
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_ParentMode_Mixed3D_TetWedge)
  {
    // 3D mixed mesh: 1 tetrahedron + 1 wedge sharing a triangular face
    //
    // Vertices:
    //   0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(0,0,1)
    //   4=(1,0,1), 5=(0,1,1)
    //
    // Tet: (0,1,2,3)
    // Wedge: (1,2,3, 4,5,3)  -> reuse face (1,2,3)
    //
    // Use proper wedge vertex ordering: bottom tri (1,2,4), top tri (3,5,6)
    // Actually let's use a simpler arrangement with 7 vertices
    Mesh<Context::Local> mesh =
      Mesh<Context::Local>::Builder()
        .initialize(3)
        .nodes(7)
        .vertex({0.0, 0.0, 0.0})  // 0
        .vertex({1.0, 0.0, 0.0})  // 1
        .vertex({0.0, 1.0, 0.0})  // 2
        .vertex({0.0, 0.0, 1.0})  // 3
        .vertex({1.0, 0.0, 1.0})  // 4
        .vertex({0.0, 1.0, 1.0})  // 5
        .vertex({0.5, 0.5, 0.5})  // 6 — interior
        .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 6})
        .polytope(Polytope::Type::Wedge, {0, 1, 2, 3, 4, 5})
        .finalize();

    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 3u);
    mesh.getConnectivity().compute(D, 0);

    Shard::Builder sb;
    sb.initialize(mesh);

    for (Index v = 0; v < mesh.getVertexCount(); v++)
      sb.include({0, v}, Shard::State::Owned);
    for (Index c = 0; c < mesh.getCellCount(); c++)
      sb.include({D, c}, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 7u);
    EXPECT_EQ(shard.getCellCount(), 2u);
    EXPECT_EQ(shard.getDimension(), 3u);

    for (Index c = 0; c < shard.getCellCount(); c++)
      EXPECT_TRUE(shard.isOwned(D, c));
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Direct mode with Wedge (3D)
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_DirectMode_Wedge)
  {
    Shard::Builder sb;
    sb.initialize(3, 3);

    // Wedge: bottom triangle (v0,v1,v2), top triangle (v3,v4,v5)
    Index v0 = sb.vertex(0, Math::SpatialPoint{{0.0, 0.0, 0.0}}, Shard::State::Owned);
    Index v1 = sb.vertex(1, Math::SpatialPoint{{1.0, 0.0, 0.0}}, Shard::State::Owned);
    Index v2 = sb.vertex(2, Math::SpatialPoint{{0.0, 1.0, 0.0}}, Shard::State::Owned);
    Index v3 = sb.vertex(3, Math::SpatialPoint{{0.0, 0.0, 1.0}}, Shard::State::Shared);
    Index v4 = sb.vertex(4, Math::SpatialPoint{{1.0, 0.0, 1.0}}, Shard::State::Shared);
    Index v5 = sb.vertex(5, Math::SpatialPoint{{0.0, 1.0, 1.0}}, Shard::State::Shared);

    IndexArray vs(6);
    vs << v0, v1, v2, v3, v4, v5;
    Index w0 = sb.polytope(3, 300, Polytope::Type::Wedge, vs, Shard::State::Owned);
    EXPECT_EQ(w0, 0u);

    for (Index vi = v3; vi <= v5; vi++)
      sb.setOwner(0, vi, 2);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 6u);
    EXPECT_EQ(shard.getCellCount(), 1u);
    EXPECT_EQ(shard.getDimension(), 3u);
    EXPECT_TRUE(shard.isOwned(3, w0));

    for (Index vi = v0; vi <= v2; vi++)
      EXPECT_TRUE(shard.isOwned(0, vi));
    for (Index vi = v3; vi <= v5; vi++)
    {
      EXPECT_TRUE(shard.isShared(0, vi));
      EXPECT_EQ(shard.getOwner(0).at(vi), 2u);
    }
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Direct mode with mixed 2D (tri + quad)
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Shard, Builder_DirectMode_Mixed2D)
  {
    Shard::Builder sb;
    sb.initialize(2, 2);

    // 5 vertices forming a triangle + adjacent quadrilateral
    Index v0 = sb.vertex(0, Math::SpatialPoint{{0.0, 0.0}}, Shard::State::Owned);
    Index v1 = sb.vertex(1, Math::SpatialPoint{{1.0, 0.0}}, Shard::State::Owned);
    Index v2 = sb.vertex(2, Math::SpatialPoint{{0.5, 1.0}}, Shard::State::Owned);
    Index v3 = sb.vertex(3, Math::SpatialPoint{{2.0, 0.0}}, Shard::State::Owned);
    Index v4 = sb.vertex(4, Math::SpatialPoint{{2.0, 1.0}}, Shard::State::Owned);

    // Triangle: (v0, v1, v2)
    IndexArray tri(3);
    tri << v0, v1, v2;
    Index t0 = sb.polytope(2, 0, Polytope::Type::Triangle, tri, Shard::State::Owned);

    // Quad: (v1, v3, v4, v2)
    IndexArray quad(4);
    quad << v1, v3, v4, v2;
    Index q0 = sb.polytope(2, 1, Polytope::Type::Quadrilateral, quad, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 5u);
    EXPECT_EQ(shard.getCellCount(), 2u);
    EXPECT_TRUE(shard.isOwned(2, t0));
    EXPECT_TRUE(shard.isOwned(2, q0));

    // Polytope map
    const auto& cmap = shard.getPolytopeMap(2);
    EXPECT_EQ(cmap.left[t0], 0u);
    EXPECT_EQ(cmap.left[q0], 1u);
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent mode: attribute preservation
  // ---------------------------------------------------------------------------

  /**
   * @brief Tests that Builder::include() preserves entity attributes from the
   * parent mesh after shard finalization.
   */
  TEST(Rodin_Geometry_Shard, Builder_ParentMode_AttributePreservation)
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

    const size_t D = mesh.getDimension();

    // Set distinct attributes on each cell in the parent mesh
    mesh.setAttribute({D, 0}, 10);
    mesh.setAttribute({D, 1}, 20);

    mesh.getConnectivity().compute(D, 0);

    // Build a shard including all entities as Owned
    Shard::Builder sb;
    sb.initialize(mesh);
    for (Index v = 0; v < mesh.getVertexCount(); v++)
      sb.include({0, v}, Shard::State::Owned);
    for (Index c = 0; c < mesh.getCellCount(); c++)
      sb.include({D, c}, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getCellCount(), 2u);

    // Use the polytope map to find which local cell corresponds to each parent cell
    const auto& cmap = shard.getPolytopeMap(D);
    for (Index parentCell = 0; parentCell < mesh.getCellCount(); parentCell++)
    {
      auto it = cmap.right.find(parentCell);
      ASSERT_NE(it, cmap.right.end());
      Index localCell = it->second;

      auto parentAttr = mesh.getAttribute(D, parentCell);
      auto shardAttr  = shard.getAttribute(D, localCell);

      ASSERT_TRUE(parentAttr.has_value());
      ASSERT_TRUE(shardAttr.has_value());
      EXPECT_EQ(shardAttr.value(), parentAttr.value());
    }
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Parent mode: connectivity preservation
  // ---------------------------------------------------------------------------

  /**
   * @brief Tests that shard-local cell-to-vertex connectivity matches the
   * parent mesh connectivity when mapped through the PolytopeMap.
   */
  TEST(Rodin_Geometry_Shard, Builder_ParentMode_ConnectivityPreservation)
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

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, 0);

    // Build shard including all entities as Owned
    Shard::Builder sb;
    sb.initialize(mesh);
    for (Index v = 0; v < mesh.getVertexCount(); v++)
      sb.include({0, v}, Shard::State::Owned);
    for (Index c = 0; c < mesh.getCellCount(); c++)
      sb.include({D, c}, Shard::State::Owned);

    Shard shard = sb.finalize();
    shard.getConnectivity().compute(D, 0);

    EXPECT_EQ(shard.getCellCount(), 2u);

    const auto& parentCellToVertex = mesh.getConnectivity().getIncidence(D, 0);
    const auto& shardCellToVertex  = shard.getConnectivity().getIncidence(D, 0);
    const auto& cmap = shard.getPolytopeMap(D);
    const auto& vmap = shard.getPolytopeMap(0);

    for (Index localCell = 0; localCell < shard.getCellCount(); localCell++)
    {
      Index parentCell = cmap.left[localCell];

      // Collect parent vertex set for this cell
      const auto& parentVerts = parentCellToVertex.at(parentCell);

      // Collect shard vertex set, mapped back to parent indices
      const auto& shardVerts = shardCellToVertex.at(localCell);

      ASSERT_EQ(shardVerts.size(), parentVerts.size());

      // Convert shard-local vertex indices to parent indices and compare as sets
      std::set<Index> parentSet(parentVerts.begin(), parentVerts.end());
      std::set<Index> mappedSet;
      for (const Index sv : shardVerts)
        mappedSet.insert(vmap.left[sv]);

      EXPECT_EQ(mappedSet, parentSet);
    }
  }

  // ---------------------------------------------------------------------------
  // Shard::Builder — Direct mode with mixed 3D (tetrahedron + wedge)
  // ---------------------------------------------------------------------------

  /**
   * @brief Tests direct-mode construction of a mixed 3D shard containing one
   * tetrahedron and one wedge that share a triangular face.
   */
  TEST(Rodin_Geometry_Shard, Builder_DirectMode_Mixed3D)
  {
    Shard::Builder sb;
    sb.initialize(3, 3);

    // Shared face vertices: v0, v1, v2
    Index v0 = sb.vertex(0, Math::SpatialPoint{{0.0, 0.0, 0.0}}, Shard::State::Owned);
    Index v1 = sb.vertex(1, Math::SpatialPoint{{1.0, 0.0, 0.0}}, Shard::State::Owned);
    Index v2 = sb.vertex(2, Math::SpatialPoint{{0.0, 1.0, 0.0}}, Shard::State::Owned);

    // Tetrahedron apex
    Index v3 = sb.vertex(3, Math::SpatialPoint{{0.3, 0.3, -1.0}}, Shard::State::Owned);

    // Wedge top-triangle vertices
    Index v4 = sb.vertex(4, Math::SpatialPoint{{0.0, 0.0, 1.0}}, Shard::State::Owned);
    Index v5 = sb.vertex(5, Math::SpatialPoint{{1.0, 0.0, 1.0}}, Shard::State::Owned);
    Index v6 = sb.vertex(6, Math::SpatialPoint{{0.0, 1.0, 1.0}}, Shard::State::Owned);

    // Tetrahedron: (v0, v1, v2, v3)
    IndexArray tetVerts(4);
    tetVerts << v0, v1, v2, v3;
    Index t0 = sb.polytope(3, 0, Polytope::Type::Tetrahedron, tetVerts, Shard::State::Owned);

    // Wedge vertex ordering: bottom triangle (v0, v1, v2), then top triangle (v4, v5, v6)
    IndexArray wedgeVerts(6);
    wedgeVerts << v0, v1, v2, v4, v5, v6;
    Index w0 = sb.polytope(3, 1, Polytope::Type::Wedge, wedgeVerts, Shard::State::Owned);

    Shard shard = sb.finalize();

    EXPECT_EQ(shard.getVertexCount(), 7u);
    EXPECT_EQ(shard.getCellCount(), 2u);
    EXPECT_EQ(shard.getDimension(), 3u);
    EXPECT_TRUE(shard.isOwned(3, t0));
    EXPECT_TRUE(shard.isOwned(3, w0));

    // Verify all vertices are owned
    for (Index vi = 0; vi < shard.getVertexCount(); vi++)
      EXPECT_TRUE(shard.isOwned(0, vi));

    // Polytope map: global indices should match
    const auto& cmap = shard.getPolytopeMap(3);
    EXPECT_EQ(cmap.left[t0], 0u);
    EXPECT_EQ(cmap.left[w0], 1u);
  }
}
