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
#include <Rodin/Geometry/Sharder.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ---------------------------------------------------------------------------
  // Helper: Prepare a mesh with the required incidences for sharding.
  // ---------------------------------------------------------------------------

  static Mesh<Context::Local> makeShardableMesh(Polytope::Type type,
                                                std::initializer_list<size_t> shape)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(type, shape);
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    return mesh;
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Basic sharding produces the correct number of shards.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, BasicSharding_TriangleMesh_2Parts)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    EXPECT_EQ(shards.size(), 2u);

    // Each shard must have at least one cell
    for (const auto& s : shards)
    {
      EXPECT_GT(s.getCellCount(), 0u);
      EXPECT_GT(s.getVertexCount(), 0u);
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Every parent cell is owned by exactly one shard.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, EveryCellOwnedByExactlyOneShard)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(3);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();

    // Collect all parent cell indices marked Owned, exactly once
    std::set<Index> ownedCells;
    for (size_t si = 0; si < shards.size(); si++)
    {
      const auto& shard = shards[si];
      const auto& cmap = shard.getPolytopeMap(D);
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (shard.isOwned(D, ci))
        {
          Index parentIdx = cmap.left[ci];
          auto [it, inserted] = ownedCells.insert(parentIdx);
          EXPECT_TRUE(inserted) << "Cell " << parentIdx
                                << " is owned by more than one shard.";
        }
      }
    }

    // All parent cells must be accounted for
    EXPECT_EQ(ownedCells.size(), mesh.getCellCount());
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Every parent vertex appears in at least one shard.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, EveryVertexPresent)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Triangle, {3, 3});

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();

    std::set<Index> seenVertices;
    for (const auto& shard : shards)
    {
      const auto& vmap = shard.getPolytopeMap(0);
      for (Index vi = 0; vi < shard.getVertexCount(); vi++)
        seenVertices.insert(vmap.left[vi]);
    }

    EXPECT_EQ(seenVertices.size(), mesh.getVertexCount());
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Every vertex is owned by exactly one shard.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, EveryVertexOwnedByExactlyOneShard)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(3);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();

    std::set<Index> ownedVertices;
    for (const auto& shard : shards)
    {
      const auto& vmap = shard.getPolytopeMap(0);
      for (Index vi = 0; vi < shard.getVertexCount(); vi++)
      {
        if (shard.isOwned(0, vi))
        {
          Index parentIdx = vmap.left[vi];
          auto [it, inserted] = ownedVertices.insert(parentIdx);
          EXPECT_TRUE(inserted) << "Vertex " << parentIdx
                                << " is owned by more than one shard.";
        }
      }
    }

    EXPECT_EQ(ownedVertices.size(), mesh.getVertexCount());
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Ghost cells belong to a neighboring partition.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, GhostCellsAreNeighbors)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();

    // For each shard, every ghost cell must be adjacent (in the parent mesh)
    // to at least one owned cell in that shard.
    const auto& parentAdj = mesh.getConnectivity().getIncidence(D, D);

    for (size_t si = 0; si < shards.size(); si++)
    {
      const auto& shard = shards[si];
      const auto& cmap = shard.getPolytopeMap(D);

      // Collect parent-indices of owned cells in this shard
      std::set<Index> ownedParent;
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (shard.isOwned(D, ci))
          ownedParent.insert(cmap.left[ci]);
      }

      // Check ghost cells
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (!shard.isGhost(D, ci))
          continue;

        Index parentGhost = cmap.left[ci];
        const auto& neighbors = parentAdj.at(parentGhost);
        bool adjacentToOwned = false;
        for (const Index nbr : neighbors)
        {
          if (ownedParent.count(nbr))
          {
            adjacentToOwned = true;
            break;
          }
        }
        EXPECT_TRUE(adjacentToOwned)
          << "Ghost cell " << parentGhost << " in shard " << si
          << " is not adjacent to any owned cell.";
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Owner map points to a valid partition for non-owned entities.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, OwnerMapConsistency)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(3);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    const size_t numShards = shards.size();

    for (size_t si = 0; si < numShards; si++)
    {
      const auto& shard = shards[si];
      for (size_t d = 0; d <= D; d++)
      {
        const auto& ownerMap = shard.getOwner(d);
        for (const auto& [localIdx, ownerRank] : ownerMap)
        {
          EXPECT_LT(ownerRank, numShards)
            << "Owner rank " << ownerRank << " out of range for shard " << si
            << " dim " << d << " local idx " << localIdx;
          EXPECT_NE(ownerRank, si)
            << "Non-owned entity in shard " << si << " should not point to itself.";
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Halo map consistency: remote ranks are valid and the entity
  // is present in those remote shards.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, HaloMapConsistency)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(3);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    const size_t numShards = shards.size();

    for (size_t si = 0; si < numShards; si++)
    {
      const auto& shard = shards[si];
      for (size_t d = 0; d <= D; d++)
      {
        const auto& haloMap = shard.getHalo(d);
        const auto& pmap = shard.getPolytopeMap(d);
        for (const auto& [localIdx, remoteRanks] : haloMap)
        {
          // The entity must be owned in this shard
          EXPECT_TRUE(shard.isOwned(d, localIdx))
            << "Halo entry at shard " << si << " dim " << d
            << " local " << localIdx << " is not owned.";

          Index parentIdx = pmap.left[localIdx];

          for (const Index remoteRank : remoteRanks)
          {
            EXPECT_LT(remoteRank, numShards)
              << "Halo rank " << remoteRank << " out of range.";
            EXPECT_NE(remoteRank, si)
              << "Halo should not reference the owning shard itself.";

            // The entity should be present in the remote shard
            const auto& remotePmap = shards[remoteRank].getPolytopeMap(d);
            auto it = remotePmap.right.find(parentIdx);
            EXPECT_NE(it, remotePmap.right.end())
              << "Entity " << parentIdx << " (dim " << d
              << ") referenced in halo of shard " << si
              << " but not found in remote shard " << remoteRank;
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — PolytopeMap is bidirectionally consistent in every shard.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, PolytopeMapConsistency)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Triangle, {3, 3});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();

    for (size_t si = 0; si < shards.size(); si++)
    {
      const auto& shard = shards[si];
      for (size_t d = 0; d <= D; d++)
      {
        const auto& pmap = shard.getPolytopeMap(d);
        EXPECT_EQ(pmap.left.size(), pmap.right.size())
          << "PolytopeMap size mismatch in shard " << si << " dim " << d;

        for (Index li = 0; li < pmap.left.size(); li++)
        {
          Index parentIdx = pmap.left[li];
          auto it = pmap.right.find(parentIdx);
          ASSERT_NE(it, pmap.right.end())
            << "Parent idx " << parentIdx << " not found in right map, shard "
            << si << " dim " << d;
          EXPECT_EQ(it->second, li);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Quadrilateral mesh sharding works.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, QuadrilateralMesh)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Quadrilateral, {4, 4});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    EXPECT_EQ(shards.size(), 2u);

    // All cells must be owned by exactly one shard
    std::set<Index> ownedCells;
    for (const auto& shard : shards)
    {
      const auto& cmap = shard.getPolytopeMap(D);
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (shard.isOwned(D, ci))
        {
          auto [it, ins] = ownedCells.insert(cmap.left[ci]);
          EXPECT_TRUE(ins);
        }
      }
    }
    EXPECT_EQ(ownedCells.size(), mesh.getCellCount());
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Multiple partition counts (2, 3, 4).
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, MultiplePartitionCounts)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Triangle, {6, 6});
    const size_t D = mesh.getDimension();
    const size_t totalCells = mesh.getCellCount();

    for (size_t numParts : {2u, 3u, 4u})
    {
      BalancedCompactPartitioner partitioner(mesh);
      partitioner.partition(numParts);

      Context::Local ctx;
      SharderBase<Context::Local> sharder(ctx);
      sharder.shard(partitioner);

      const auto& shards = sharder.getShards();
      EXPECT_EQ(shards.size(), numParts);

      size_t totalOwned = 0;
      for (const auto& shard : shards)
      {
        for (Index ci = 0; ci < shard.getCellCount(); ci++)
        {
          if (shard.isOwned(D, ci))
            totalOwned++;
        }
      }
      EXPECT_EQ(totalOwned, totalCells)
        << "Total owned cells mismatch with " << numParts << " partitions.";
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — isLocal == (isOwned || isShared) for all entities.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, IsLocalConsistency)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(3);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();

    for (size_t si = 0; si < shards.size(); si++)
    {
      const auto& shard = shards[si];
      for (size_t d = 0; d <= D; d++)
      {
        const auto& state = shard.getState(d);
        for (Index i = 0; i < state.size(); i++)
        {
          bool local = shard.isLocal(d, i);
          bool ownedOrShared = shard.isOwned(d, i) || shard.isShared(d, i);
          EXPECT_EQ(local, ownedOrShared)
            << "isLocal mismatch at shard " << si
            << " dim " << d << " idx " << i;

          // Exactly one state must be true
          int count = static_cast<int>(shard.isOwned(d, i))
                    + static_cast<int>(shard.isShared(d, i))
                    + static_cast<int>(shard.isGhost(d, i));
          EXPECT_EQ(count, 1)
            << "Entity at shard " << si << " dim " << d << " idx " << i
            << " has " << count << " states (expected exactly 1).";
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — getContext() returns the context passed at construction.
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, GetContext)
  {
    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    // Just ensure it compiles and returns a reference
    const auto& ref = sharder.getContext();
    (void)ref;
  }

  // ---------------------------------------------------------------------------
  // SharderBase — 1D Segment mesh sharding
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, SegmentMesh_BasicSharding)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Segment, {8});
    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 1u);

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    EXPECT_EQ(shards.size(), 2u);

    // Every cell owned exactly once
    std::set<Index> ownedCells;
    for (const auto& shard : shards)
    {
      const auto& cmap = shard.getPolytopeMap(D);
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (shard.isOwned(D, ci))
        {
          auto [it, ins] = ownedCells.insert(cmap.left[ci]);
          EXPECT_TRUE(ins);
        }
      }
    }
    EXPECT_EQ(ownedCells.size(), mesh.getCellCount());

    // Every vertex present
    std::set<Index> seenVerts;
    for (const auto& shard : shards)
    {
      const auto& vmap = shard.getPolytopeMap(0);
      for (Index vi = 0; vi < shard.getVertexCount(); vi++)
        seenVerts.insert(vmap.left[vi]);
    }
    EXPECT_EQ(seenVerts.size(), mesh.getVertexCount());
  }

  TEST(Rodin_Geometry_Sharder, SegmentMesh_IsLocalConsistency)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Segment, {10});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(3);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    for (const auto& shard : sharder.getShards())
    {
      for (size_t d = 0; d <= D; d++)
      {
        const auto& state = shard.getState(d);
        for (Index i = 0; i < state.size(); i++)
        {
          bool local = shard.isLocal(d, i);
          bool ownedOrShared = shard.isOwned(d, i) || shard.isShared(d, i);
          EXPECT_EQ(local, ownedOrShared);

          int count = static_cast<int>(shard.isOwned(d, i))
                    + static_cast<int>(shard.isShared(d, i))
                    + static_cast<int>(shard.isGhost(d, i));
          EXPECT_EQ(count, 1);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — 3D Tetrahedron mesh sharding
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, TetrahedronMesh_BasicSharding)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Tetrahedron, {3, 3, 3});
    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 3u);

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    EXPECT_EQ(shards.size(), 2u);

    for (const auto& s : shards)
    {
      EXPECT_GT(s.getCellCount(), 0u);
      EXPECT_GT(s.getVertexCount(), 0u);
    }
  }

  TEST(Rodin_Geometry_Sharder, TetrahedronMesh_CellOwnershipUniqueness)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Tetrahedron, {3, 3, 3});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(3);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    std::set<Index> ownedCells;
    for (const auto& shard : sharder.getShards())
    {
      const auto& cmap = shard.getPolytopeMap(D);
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (shard.isOwned(D, ci))
        {
          auto [it, ins] = ownedCells.insert(cmap.left[ci]);
          EXPECT_TRUE(ins);
        }
      }
    }
    EXPECT_EQ(ownedCells.size(), mesh.getCellCount());
  }

  TEST(Rodin_Geometry_Sharder, TetrahedronMesh_GhostCellsAreNeighbors)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Tetrahedron, {3, 3, 3});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& parentAdj = mesh.getConnectivity().getIncidence(D, D);

    for (size_t si = 0; si < sharder.getShards().size(); si++)
    {
      const auto& shard = sharder.getShards()[si];
      const auto& cmap = shard.getPolytopeMap(D);

      std::set<Index> ownedParent;
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (shard.isOwned(D, ci))
          ownedParent.insert(cmap.left[ci]);
      }

      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (!shard.isGhost(D, ci))
          continue;
        Index parentGhost = cmap.left[ci];
        const auto& neighbors = parentAdj.at(parentGhost);
        bool adjacentToOwned = false;
        for (const Index nbr : neighbors)
        {
          if (ownedParent.count(nbr))
          {
            adjacentToOwned = true;
            break;
          }
        }
        EXPECT_TRUE(adjacentToOwned)
          << "Ghost tet " << parentGhost << " in shard " << si
          << " is not adjacent to any owned cell.";
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — 3D Hexahedron mesh sharding
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, HexahedronMesh_BasicSharding)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Hexahedron, {3, 3, 3});
    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 3u);

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    EXPECT_EQ(shards.size(), 2u);

    std::set<Index> ownedCells;
    for (const auto& shard : shards)
    {
      const auto& cmap = shard.getPolytopeMap(D);
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (shard.isOwned(D, ci))
        {
          auto [it, ins] = ownedCells.insert(cmap.left[ci]);
          EXPECT_TRUE(ins);
        }
      }
    }
    EXPECT_EQ(ownedCells.size(), mesh.getCellCount());
  }

  TEST(Rodin_Geometry_Sharder, HexahedronMesh_OwnerMapConsistency)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Hexahedron, {3, 3, 3});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(3);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    const size_t numShards = shards.size();

    for (size_t si = 0; si < numShards; si++)
    {
      const auto& shard = shards[si];
      for (size_t d = 0; d <= D; d++)
      {
        for (const auto& [localIdx, ownerRank] : shard.getOwner(d))
        {
          EXPECT_LT(ownerRank, numShards);
          EXPECT_NE(ownerRank, si);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — 3D Wedge mesh sharding
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, WedgeMesh_BasicSharding)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Wedge, {3, 3, 3});
    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 3u);

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    EXPECT_EQ(shards.size(), 2u);

    std::set<Index> ownedCells;
    for (const auto& shard : shards)
    {
      const auto& cmap = shard.getPolytopeMap(D);
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (shard.isOwned(D, ci))
        {
          auto [it, ins] = ownedCells.insert(cmap.left[ci]);
          EXPECT_TRUE(ins);
        }
      }
    }
    EXPECT_EQ(ownedCells.size(), mesh.getCellCount());
  }

  TEST(Rodin_Geometry_Sharder, WedgeMesh_HaloMapConsistency)
  {
    auto mesh = makeShardableMesh(Polytope::Type::Wedge, {3, 3, 3});
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();

    for (size_t si = 0; si < shards.size(); si++)
    {
      const auto& shard = shards[si];
      for (size_t d = 0; d <= D; d++)
      {
        const auto& haloMap = shard.getHalo(d);
        const auto& pmap = shard.getPolytopeMap(d);
        for (const auto& [localIdx, remoteRanks] : haloMap)
        {
          EXPECT_TRUE(shard.isOwned(d, localIdx));
          Index parentIdx = pmap.left[localIdx];
          for (const Index remoteRank : remoteRanks)
          {
            EXPECT_LT(remoteRank, shards.size());
            EXPECT_NE(remoteRank, si);
            const auto& remotePmap = shards[remoteRank].getPolytopeMap(d);
            EXPECT_NE(remotePmap.right.find(parentIdx), remotePmap.right.end());
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Mixed 2D mesh (triangles + quadrilaterals) sharding
  // ---------------------------------------------------------------------------

  /**
   * @brief Helper to build a mixed 2D mesh with triangles and quadrilaterals.
   *
   * Layout (6 vertices, 3 cells):
   *   2---3---4
   *   |\  | Q |
   *   | \ |   |
   *   0--1----5
   *
   * T0: (0,1,2), T1: (1,3,2), Q0: (1,5,4,3)
   */
  static Mesh<Context::Local> makeMixed2DMesh()
  {
    auto mesh = Mesh<Context::Local>::Builder()
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
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    return mesh;
  }

  TEST(Rodin_Geometry_Sharder, Mixed2D_BasicSharding)
  {
    auto mesh = makeMixed2DMesh();
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    EXPECT_EQ(shards.size(), 2u);

    // Every cell owned exactly once
    std::set<Index> ownedCells;
    for (const auto& shard : shards)
    {
      const auto& cmap = shard.getPolytopeMap(D);
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (shard.isOwned(D, ci))
        {
          auto [it, ins] = ownedCells.insert(cmap.left[ci]);
          EXPECT_TRUE(ins);
        }
      }
    }
    EXPECT_EQ(ownedCells.size(), mesh.getCellCount());
  }

  TEST(Rodin_Geometry_Sharder, Mixed2D_VertexOwnershipUniqueness)
  {
    auto mesh = makeMixed2DMesh();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    std::set<Index> ownedVerts;
    for (const auto& shard : sharder.getShards())
    {
      const auto& vmap = shard.getPolytopeMap(0);
      for (Index vi = 0; vi < shard.getVertexCount(); vi++)
      {
        if (shard.isOwned(0, vi))
        {
          auto [it, ins] = ownedVerts.insert(vmap.left[vi]);
          EXPECT_TRUE(ins);
        }
      }
    }
    EXPECT_EQ(ownedVerts.size(), mesh.getVertexCount());
  }

  TEST(Rodin_Geometry_Sharder, Mixed2D_PolytopeMapConsistency)
  {
    auto mesh = makeMixed2DMesh();
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    for (size_t si = 0; si < sharder.getShards().size(); si++)
    {
      const auto& shard = sharder.getShards()[si];
      for (size_t d = 0; d <= D; d++)
      {
        const auto& pmap = shard.getPolytopeMap(d);
        EXPECT_EQ(pmap.left.size(), pmap.right.size());
        for (Index li = 0; li < pmap.left.size(); li++)
        {
          Index parentIdx = pmap.left[li];
          auto it = pmap.right.find(parentIdx);
          ASSERT_NE(it, pmap.right.end());
          EXPECT_EQ(it->second, li);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Mixed 3D mesh (tetrahedron + wedge) sharding
  // ---------------------------------------------------------------------------

  /**
   * @brief Helper to build a mixed 3D mesh with tetrahedra and wedges.
   *
   * Tet: (0,1,2,6)
   * Wedge: (0,1,2,3,4,5) — prism extruding the base triangle upward
   */
  static Mesh<Context::Local> makeMixed3DMesh()
  {
    auto mesh = Mesh<Context::Local>::Builder()
      .initialize(3)
      .nodes(7)
      .vertex({0.0, 0.0, 0.0})  // 0
      .vertex({1.0, 0.0, 0.0})  // 1
      .vertex({0.0, 1.0, 0.0})  // 2
      .vertex({0.0, 0.0, 1.0})  // 3
      .vertex({1.0, 0.0, 1.0})  // 4
      .vertex({0.0, 1.0, 1.0})  // 5
      .vertex({0.5, 0.5, 0.5})  // 6 — interior point
      .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 6})
      .polytope(Polytope::Type::Wedge, {0, 1, 2, 3, 4, 5})
      .finalize();

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    return mesh;
  }

  TEST(Rodin_Geometry_Sharder, Mixed3D_BasicSharding)
  {
    auto mesh = makeMixed3DMesh();
    const size_t D = mesh.getDimension();
    EXPECT_EQ(D, 3u);

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    const auto& shards = sharder.getShards();
    EXPECT_EQ(shards.size(), 2u);

    // Every cell owned exactly once
    std::set<Index> ownedCells;
    for (const auto& shard : shards)
    {
      const auto& cmap = shard.getPolytopeMap(D);
      for (Index ci = 0; ci < shard.getCellCount(); ci++)
      {
        if (shard.isOwned(D, ci))
        {
          auto [it, ins] = ownedCells.insert(cmap.left[ci]);
          EXPECT_TRUE(ins);
        }
      }
    }
    EXPECT_EQ(ownedCells.size(), mesh.getCellCount());
  }

  TEST(Rodin_Geometry_Sharder, Mixed3D_VertexCoverage)
  {
    auto mesh = makeMixed3DMesh();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    std::set<Index> seenVerts;
    for (const auto& shard : sharder.getShards())
    {
      const auto& vmap = shard.getPolytopeMap(0);
      for (Index vi = 0; vi < shard.getVertexCount(); vi++)
        seenVerts.insert(vmap.left[vi]);
    }
    EXPECT_EQ(seenVerts.size(), mesh.getVertexCount());
  }

  TEST(Rodin_Geometry_Sharder, Mixed3D_IsLocalConsistency)
  {
    auto mesh = makeMixed3DMesh();
    const size_t D = mesh.getDimension();

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);

    Context::Local ctx;
    SharderBase<Context::Local> sharder(ctx);
    sharder.shard(partitioner);

    for (const auto& shard : sharder.getShards())
    {
      for (size_t d = 0; d <= D; d++)
      {
        const auto& state = shard.getState(d);
        for (Index i = 0; i < state.size(); i++)
        {
          bool local = shard.isLocal(d, i);
          bool ownedOrShared = shard.isOwned(d, i) || shard.isShared(d, i);
          EXPECT_EQ(local, ownedOrShared);

          int count = static_cast<int>(shard.isOwned(d, i))
                    + static_cast<int>(shard.isShared(d, i))
                    + static_cast<int>(shard.isGhost(d, i));
          EXPECT_EQ(count, 1);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — Larger mixed 2D mesh with more cells
  // ---------------------------------------------------------------------------

  /**
   * @brief Builds a larger mixed 2D mesh: 4 triangles + 2 quadrilaterals.
   *
   * Layout (8 vertices):
   *   5---6---7
   *   | Q1| Q0|
   *   2---3---4
   *   |\ T1  /|
   *   |T0\ /T3|
   *   0---1
   *       |\ T2
   *       (uses vertex 4 from above)
   *
   * Actually, simpler to just build a strip:
   *   4---5---6---7
   *   |T1/|  Q1  |
   *   | / | Q0|  |
   *   0---1---2---3
   */
  static Mesh<Context::Local> makeLargerMixed2DMesh()
  {
    auto mesh = Mesh<Context::Local>::Builder()
      .initialize(2)
      .nodes(8)
      .vertex({0.0, 0.0})  // 0
      .vertex({1.0, 0.0})  // 1
      .vertex({2.0, 0.0})  // 2
      .vertex({3.0, 0.0})  // 3
      .vertex({0.0, 1.0})  // 4
      .vertex({1.0, 1.0})  // 5
      .vertex({2.0, 1.0})  // 6
      .vertex({3.0, 1.0})  // 7
      .polytope(Polytope::Type::Triangle, {0, 1, 4})
      .polytope(Polytope::Type::Triangle, {1, 5, 4})
      .polytope(Polytope::Type::Quadrilateral, {1, 2, 6, 5})
      .polytope(Polytope::Type::Quadrilateral, {2, 3, 7, 6})
      .finalize();

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    return mesh;
  }

  TEST(Rodin_Geometry_Sharder, LargerMixed2D_MultiplePartitions)
  {
    auto mesh = makeLargerMixed2DMesh();
    const size_t D = mesh.getDimension();
    const size_t totalCells = mesh.getCellCount();
    EXPECT_EQ(totalCells, 4u);

    for (size_t numParts : {2u, 3u})
    {
      BalancedCompactPartitioner partitioner(mesh);
      partitioner.partition(numParts);

      Context::Local ctx;
      SharderBase<Context::Local> sharder(ctx);
      sharder.shard(partitioner);

      const auto& shards = sharder.getShards();
      EXPECT_EQ(shards.size(), numParts);

      size_t totalOwned = 0;
      for (const auto& shard : shards)
      {
        for (Index ci = 0; ci < shard.getCellCount(); ci++)
        {
          if (shard.isOwned(D, ci))
            totalOwned++;
        }
      }
      EXPECT_EQ(totalOwned, totalCells)
        << "Total owned cells mismatch with " << numParts << " partitions.";
    }
  }

  // ---------------------------------------------------------------------------
  // SharderBase — 3D mesh with multiple partitions across types
  // ---------------------------------------------------------------------------

  TEST(Rodin_Geometry_Sharder, AllUniformGrid3D_MultiplePartitions)
  {
    // Test all 3D uniform grid types with multiple partition counts
    for (Polytope::Type type : { Polytope::Type::Tetrahedron,
                                  Polytope::Type::Hexahedron,
                                  Polytope::Type::Wedge })
    {
      auto mesh = makeShardableMesh(type, {3, 3, 3});
      const size_t D = mesh.getDimension();
      const size_t totalCells = mesh.getCellCount();

      for (size_t numParts : {2u, 3u})
      {
        BalancedCompactPartitioner partitioner(mesh);
        partitioner.partition(numParts);

        Context::Local ctx;
        SharderBase<Context::Local> sharder(ctx);
        sharder.shard(partitioner);

        const auto& shards = sharder.getShards();
        EXPECT_EQ(shards.size(), numParts);

        size_t totalOwned = 0;
        for (const auto& shard : shards)
        {
          for (Index ci = 0; ci < shard.getCellCount(); ci++)
          {
            if (shard.isOwned(D, ci))
              totalOwned++;
          }
        }
        EXPECT_EQ(totalOwned, totalCells)
          << "Type " << static_cast<int>(type) << " with " << numParts
          << " partitions: owned cell count mismatch.";
      }
    }
  }
}
