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
}
