/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <set>
#include <numeric>

#include <gtest/gtest.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/Shard.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>

using namespace Rodin;
using namespace Rodin::Geometry;

// ---------------------------------------------------------------------------
// Global MPI handles initialized in main().
// ---------------------------------------------------------------------------

static boost::mpi::environment* g_env   = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  /**
   * @brief Creates a local mesh with the incidences required for sharding.
   */
  Mesh<Context::Local> makeShardableMesh(Polytope::Type type,
                                         std::initializer_list<size_t> shape)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(type, shape);
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    return mesh;
  }

  /**
   * @brief Creates a mixed 2D mesh (triangles + quadrilaterals) with
   * the required incidences for sharding.
   *
   * Layout (8 vertices, 4 cells: 2 triangles + 2 quads):
   *   4---5---6---7
   *   |T1/|       |
   *   | / | Q0| Q1|
   *   0---1---2---3
   */
  Mesh<Context::Local> makeMixed2DMesh()
  {
    auto mesh = Mesh<Context::Local>::Builder()
      .initialize(2)
      .nodes(8)
      .vertex({0.0, 0.0})
      .vertex({1.0, 0.0})
      .vertex({2.0, 0.0})
      .vertex({3.0, 0.0})
      .vertex({0.0, 1.0})
      .vertex({1.0, 1.0})
      .vertex({2.0, 1.0})
      .vertex({3.0, 1.0})
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

  /**
   * @brief Distributes a uniform-grid mesh from root following the canonical
   * root-only shard/scatter pattern (see examples/MPI/Geometry/Sharder.cpp).
   */
  Mesh<Context::MPI> distributeFromRoot(
      const Context::MPI& ctx,
      Polytope::Type type,
      std::initializer_list<size_t> shape)
  {
    const auto& comm = ctx.getCommunicator();

    Sharder<Context::MPI> sharder(ctx);

    if (comm.rank() == 0)
    {
      auto localMesh = makeShardableMesh(type, shape);
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }

    return sharder.gather(0);
  }
}

namespace Rodin::Tests::Unit
{
  // =========================================================================
  // MPI Sharder — Triangle mesh distribution
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Triangle_EachRankHasCellsAndVertices)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();
    EXPECT_GT(shard.getCellCount(), 0u)
      << "Rank " << world.rank() << " has no cells.";
    EXPECT_GT(shard.getVertexCount(), 0u)
      << "Rank " << world.rank() << " has no vertices.";
  }

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Triangle_GlobalCellCount)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    // Compute expected total on root and broadcast
    size_t totalCells = 0;
    size_t D = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
      totalCells = localMesh.getCellCount();
      D = localMesh.getDimension();
    }
    boost::mpi::broadcast(world, totalCells, 0);
    boost::mpi::broadcast(world, D, 0);

    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();
    size_t ownedLocal = 0;
    for (Index ci = 0; ci < shard.getCellCount(); ci++)
    {
      if (shard.isOwned(D, ci))
        ownedLocal++;
    }

    size_t ownedGlobal = 0;
    boost::mpi::all_reduce(world, ownedLocal, ownedGlobal, std::plus<size_t>());

    EXPECT_EQ(ownedGlobal, totalCells);
  }

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Triangle_GlobalVertexOwnership)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    size_t totalVerts = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
      totalVerts = localMesh.getVertexCount();
    }
    boost::mpi::broadcast(world, totalVerts, 0);

    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();
    size_t ownedLocal = 0;
    for (Index vi = 0; vi < shard.getVertexCount(); vi++)
    {
      if (shard.isOwned(0, vi))
        ownedLocal++;
    }

    size_t ownedGlobal = 0;
    boost::mpi::all_reduce(world, ownedLocal, ownedGlobal, std::plus<size_t>());

    EXPECT_EQ(ownedGlobal, totalVerts);
  }

  // =========================================================================
  // MPI Sharder — isLocal invariant on each rank
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Triangle_IsLocalInvariant)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();
    const size_t D = shard.getDimension();
    for (size_t d = 0; d <= D; d++)
    {
      const auto& state = shard.getState(d);
      for (Index i = 0; i < state.size(); i++)
      {
        bool local = shard.isLocal(d, i);
        bool ownedOrShared = shard.isOwned(d, i) || shard.isShared(d, i);
        EXPECT_EQ(local, ownedOrShared)
          << "Rank " << world.rank() << " dim " << d << " idx " << i;

        int count = static_cast<int>(shard.isOwned(d, i))
                  + static_cast<int>(shard.isShared(d, i))
                  + static_cast<int>(shard.isGhost(d, i));
        EXPECT_EQ(count, 1)
          << "Rank " << world.rank() << " dim " << d << " idx " << i
          << " has " << count << " states (expected 1).";
      }
    }
  }

  // =========================================================================
  // MPI Sharder — PolytopeMap bidirectionality on each rank
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Triangle_PolytopeMapConsistency)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();
    const size_t D = shard.getDimension();
    for (size_t d = 0; d <= D; d++)
    {
      const auto& pmap = shard.getPolytopeMap(d);
      EXPECT_EQ(pmap.left.size(), pmap.right.size())
        << "Rank " << world.rank() << " dim " << d;

      for (Index li = 0; li < pmap.left.size(); li++)
      {
        Index parentIdx = pmap.left[li];
        auto it = pmap.right.find(parentIdx);
        ASSERT_NE(it, pmap.right.end())
          << "Rank " << world.rank() << " dim " << d
          << " parent " << parentIdx << " not in right map.";
        EXPECT_EQ(it->second, li);
      }
    }
  }

  // =========================================================================
  // MPI Sharder — Owner map consistency on each rank
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Triangle_OwnerMapConsistency)
  {
    const auto& world = *g_world;
    if (world.size() < 2 || world.size() > 3)
      GTEST_SKIP() << "Test requires 2 or 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();
    const size_t D = shard.getDimension();
    const size_t myRank = static_cast<size_t>(world.rank());
    const size_t numRanks = static_cast<size_t>(world.size());

    for (size_t d = 0; d <= D; d++)
    {
      const auto& ownerMap = shard.getOwner(d);
      for (const auto& [localIdx, ownerRank] : ownerMap)
      {
        EXPECT_LT(ownerRank, numRanks)
          << "Rank " << myRank << " dim " << d
          << " owner rank out of range: " << ownerRank;
        EXPECT_NE(ownerRank, myRank)
          << "Rank " << myRank << " dim " << d
          << " entity should not be owned by self in owner map.";
      }
    }
  }

  // =========================================================================
  // MPI Sharder — Dimension and space dimension consistency
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Triangle_DimensionConsistency)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    EXPECT_EQ(mpiMesh.getDimension(), 2u);
    EXPECT_EQ(mpiMesh.getSpaceDimension(), 2u);
  }

  // =========================================================================
  // MPI Sharder — Quadrilateral mesh
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Quadrilateral_GlobalCellCount)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    size_t totalCells = 0;
    size_t D = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Quadrilateral, {4, 4});
      totalCells = localMesh.getCellCount();
      D = localMesh.getDimension();
    }
    boost::mpi::broadcast(world, totalCells, 0);
    boost::mpi::broadcast(world, D, 0);

    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Quadrilateral, {4, 4});

    const auto& shard = mpiMesh.getShard();
    size_t ownedLocal = 0;
    for (Index ci = 0; ci < shard.getCellCount(); ci++)
    {
      if (shard.isOwned(D, ci))
        ownedLocal++;
    }

    size_t ownedGlobal = 0;
    boost::mpi::all_reduce(world, ownedLocal, ownedGlobal, std::plus<size_t>());

    EXPECT_EQ(ownedGlobal, totalCells);
  }

  // =========================================================================
  // MPI Sharder — 1D Segment mesh
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Segment_GlobalCellCount)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    size_t totalCells = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Segment, {10});
      totalCells = localMesh.getCellCount();
    }
    boost::mpi::broadcast(world, totalCells, 0);

    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Segment, {10});

    const auto& shard = mpiMesh.getShard();
    const size_t D = shard.getDimension();
    size_t ownedLocal = 0;
    for (Index ci = 0; ci < shard.getCellCount(); ci++)
    {
      if (shard.isOwned(D, ci))
        ownedLocal++;
    }

    size_t ownedGlobal = 0;
    boost::mpi::all_reduce(world, ownedLocal, ownedGlobal, std::plus<size_t>());

    EXPECT_EQ(ownedGlobal, totalCells);
  }

  // =========================================================================
  // MPI Sharder — 3D Tetrahedron mesh
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Tetrahedron_GlobalCellCount)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    size_t totalCells = 0;
    size_t D = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Tetrahedron, {3, 3, 3});
      totalCells = localMesh.getCellCount();
      D = localMesh.getDimension();
    }
    boost::mpi::broadcast(world, totalCells, 0);
    boost::mpi::broadcast(world, D, 0);

    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, {3, 3, 3});

    const auto& shard = mpiMesh.getShard();
    size_t ownedLocal = 0;
    for (Index ci = 0; ci < shard.getCellCount(); ci++)
    {
      if (shard.isOwned(D, ci))
        ownedLocal++;
    }

    size_t ownedGlobal = 0;
    boost::mpi::all_reduce(world, ownedLocal, ownedGlobal, std::plus<size_t>());

    EXPECT_EQ(ownedGlobal, totalCells);
  }

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Tetrahedron_DimensionConsistency)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, {3, 3, 3});

    EXPECT_EQ(mpiMesh.getDimension(), 3u);
    EXPECT_EQ(mpiMesh.getSpaceDimension(), 3u);
  }

  // =========================================================================
  // MPI Sharder — 3D Hexahedron mesh
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Hexahedron_GlobalCellCount)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    size_t totalCells = 0;
    size_t D = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Hexahedron, {3, 3, 3});
      totalCells = localMesh.getCellCount();
      D = localMesh.getDimension();
    }
    boost::mpi::broadcast(world, totalCells, 0);
    boost::mpi::broadcast(world, D, 0);

    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Hexahedron, {3, 3, 3});

    const auto& shard = mpiMesh.getShard();
    size_t ownedLocal = 0;
    for (Index ci = 0; ci < shard.getCellCount(); ci++)
    {
      if (shard.isOwned(D, ci))
        ownedLocal++;
    }

    size_t ownedGlobal = 0;
    boost::mpi::all_reduce(world, ownedLocal, ownedGlobal, std::plus<size_t>());

    EXPECT_EQ(ownedGlobal, totalCells);
  }

  // =========================================================================
  // MPI Sharder — 3D Wedge mesh
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Wedge_GlobalCellCount)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    size_t totalCells = 0;
    size_t D = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Wedge, {3, 3, 3});
      totalCells = localMesh.getCellCount();
      D = localMesh.getDimension();
    }
    boost::mpi::broadcast(world, totalCells, 0);
    boost::mpi::broadcast(world, D, 0);

    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Wedge, {3, 3, 3});

    const auto& shard = mpiMesh.getShard();
    size_t ownedLocal = 0;
    for (Index ci = 0; ci < shard.getCellCount(); ci++)
    {
      if (shard.isOwned(D, ci))
        ownedLocal++;
    }

    size_t ownedGlobal = 0;
    boost::mpi::all_reduce(world, ownedLocal, ownedGlobal, std::plus<size_t>());

    EXPECT_EQ(ownedGlobal, totalCells);
  }

  // =========================================================================
  // MPI Sharder — Mixed 2D mesh (triangles + quadrilaterals)
  // =========================================================================

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Mixed2D_GlobalCellCount)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    size_t totalCells = 0;
    size_t D = 0;
    Mesh<Context::Local> localMesh;
    if (world.rank() == 0)
    {
      localMesh = makeMixed2DMesh();
      totalCells = localMesh.getCellCount();
      D = localMesh.getDimension();
    }
    boost::mpi::broadcast(world, totalCells, 0);
    boost::mpi::broadcast(world, D, 0);

    Sharder<Context::MPI> sharder(ctx);
    if (world.rank() == 0)
    {
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(world.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mpiMesh = sharder.gather(0);

    const auto& shard = mpiMesh.getShard();
    size_t ownedLocal = 0;
    for (Index ci = 0; ci < shard.getCellCount(); ci++)
    {
      if (shard.isOwned(D, ci))
        ownedLocal++;
    }

    size_t ownedGlobal = 0;
    boost::mpi::all_reduce(world, ownedLocal, ownedGlobal, std::plus<size_t>());

    EXPECT_EQ(ownedGlobal, totalCells);
  }

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Mixed2D_GlobalVertexOwnership)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    size_t totalVerts = 0;
    Mesh<Context::Local> localMesh;
    if (world.rank() == 0)
    {
      localMesh = makeMixed2DMesh();
      totalVerts = localMesh.getVertexCount();
    }
    boost::mpi::broadcast(world, totalVerts, 0);

    Sharder<Context::MPI> sharder(ctx);
    if (world.rank() == 0)
    {
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(world.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mpiMesh = sharder.gather(0);

    const auto& shard = mpiMesh.getShard();
    size_t ownedLocal = 0;
    for (Index vi = 0; vi < shard.getVertexCount(); vi++)
    {
      if (shard.isOwned(0, vi))
        ownedLocal++;
    }

    size_t ownedGlobal = 0;
    boost::mpi::all_reduce(world, ownedLocal, ownedGlobal, std::plus<size_t>());

    EXPECT_EQ(ownedGlobal, totalVerts);
  }

  TEST(Rodin_MPI_Geometry_Sharder, Distribute_Mixed2D_IsLocalInvariant)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    Mesh<Context::Local> localMesh;
    if (world.rank() == 0)
      localMesh = makeMixed2DMesh();

    Sharder<Context::MPI> sharder(ctx);
    if (world.rank() == 0)
    {
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(world.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mpiMesh = sharder.gather(0);

    const auto& shard = mpiMesh.getShard();
    const size_t D = shard.getDimension();
    for (size_t d = 0; d <= D; d++)
    {
      const auto& state = shard.getState(d);
      for (Index i = 0; i < state.size(); i++)
      {
        bool local = shard.isLocal(d, i);
        bool ownedOrShared = shard.isOwned(d, i) || shard.isShared(d, i);
        EXPECT_EQ(local, ownedOrShared)
          << "Rank " << world.rank() << " dim " << d << " idx " << i;

        int count = static_cast<int>(shard.isOwned(d, i))
                  + static_cast<int>(shard.isShared(d, i))
                  + static_cast<int>(shard.isGhost(d, i));
        EXPECT_EQ(count, 1)
          << "Rank " << world.rank() << " dim " << d << " idx " << i;
      }
    }
  }

  // =========================================================================
  // MPI Mesh — distributed aggregate queries
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::getPolytopeCount(D) for cells equals the
   * parent mesh cell count (using the MPIMesh API directly, not manual shard
   * iteration).
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetPolytopeCount_Cells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    size_t totalCells = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
      totalCells = localMesh.getCellCount();
    }
    boost::mpi::broadcast(world, totalCells, 0);

    const size_t D = mpiMesh.getDimension();
    size_t distributedCount = mpiMesh.getPolytopeCount(D);
    EXPECT_EQ(distributedCount, totalCells)
      << "Rank " << world.rank()
      << ": getPolytopeCount(D) mismatch with parent cell count.";
  }

  /**
   * @brief Verifies that MPIMesh::getPolytopeCount(0) for vertices equals
   * the parent mesh vertex count.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetPolytopeCount_Vertices)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    size_t totalVerts = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
      totalVerts = localMesh.getVertexCount();
    }
    boost::mpi::broadcast(world, totalVerts, 0);

    size_t distributedCount = mpiMesh.getPolytopeCount(0);
    EXPECT_EQ(distributedCount, totalVerts)
      << "Rank " << world.rank()
      << ": getPolytopeCount(0) mismatch with parent vertex count.";
  }

  /**
   * @brief Verifies that MPIMesh::getPolytopeCount(Polytope::Type) correctly
   * counts cells by geometry type.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetPolytopeCount_ByType)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    size_t totalCells = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
      totalCells = localMesh.getCellCount();
    }
    boost::mpi::broadcast(world, totalCells, 0);

    size_t byType = mpiMesh.getPolytopeCount(Polytope::Type::Triangle);
    EXPECT_EQ(byType, totalCells)
      << "Rank " << world.rank()
      << ": getPolytopeCount(Triangle) mismatch.";
  }

  /**
   * @brief Verifies that getLocalIndex and getGlobalIndex are inverses
   * of each other for all local polytopes.
   */
  TEST(Rodin_MPI_Geometry_Mesh, IndexMapping_RoundTrip)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();
    const size_t D = mpiMesh.getDimension();

    for (size_t d = 0; d <= D; d++)
    {
      for (Index i = 0; i < shard.getPolytopeCount(d); i++)
      {
        Index globalIdx = mpiMesh.getGlobalIndex(d, i);
        auto localOpt = mpiMesh.getLocalIndex(d, globalIdx);
        ASSERT_TRUE(localOpt.has_value())
          << "Rank " << world.rank() << " dim " << d
          << ": getLocalIndex returned nullopt for global " << globalIdx;
        EXPECT_EQ(*localOpt, i)
          << "Rank " << world.rank() << " dim " << d
          << ": round-trip mismatch for local " << i;
      }
    }
  }

  /**
   * @brief Verifies that getLocalIndex returns nullopt for a global index
   * that is not present in the local shard.
   */
  TEST(Rodin_MPI_Geometry_Mesh, IndexMapping_MissingReturnsNullopt)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    if (world.size() < 2)
      GTEST_SKIP() << "Need at least 2 ranks to test missing indices.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const size_t D = mpiMesh.getDimension();

    // Use a very large global index that shouldn't exist
    auto localOpt = mpiMesh.getLocalIndex(D, 999999);
    EXPECT_FALSE(localOpt.has_value())
      << "Rank " << world.rank()
      << ": getLocalIndex should return nullopt for non-existent global index.";
  }

  /**
   * @brief Verifies that MPIMesh::getContext() returns the construction
   * context.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetContext)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& meshCtx = mpiMesh.getContext();
    EXPECT_EQ(meshCtx.getCommunicator().rank(), world.rank());
    EXPECT_EQ(meshCtx.getCommunicator().size(), world.size());
  }

  /**
   * @brief Verifies that local cell and vertex iterators on MPIMesh
   * work correctly and iterate over the expected number of entities.
   */
  TEST(Rodin_MPI_Geometry_Mesh, LocalIterators)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();

    // Count cells via iterator
    size_t cellCount = 0;
    for (auto it = mpiMesh.getCell(); it; ++it)
      cellCount++;
    EXPECT_EQ(cellCount, shard.getCellCount())
      << "Rank " << world.rank()
      << ": cell iterator count mismatch with shard cell count.";

    // Count vertices via iterator
    size_t vertexCount = 0;
    for (auto it = mpiMesh.getVertex(); it; ++it)
      vertexCount++;
    EXPECT_EQ(vertexCount, shard.getVertexCount())
      << "Rank " << world.rank()
      << ": vertex iterator count mismatch with shard vertex count.";
  }

  /**
   * @brief Verifies that MPIMesh::getGeometry() returns correct geometry
   * types for local cells.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetGeometry_Cells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const size_t D = mpiMesh.getDimension();
    const auto& shard = mpiMesh.getShard();

    for (Index i = 0; i < shard.getCellCount(); i++)
    {
      auto geom = mpiMesh.getGeometry(D, i);
      EXPECT_EQ(geom, Polytope::Type::Triangle)
        << "Rank " << world.rank() << " cell " << i
        << ": expected Triangle geometry.";
    }
  }

  /**
   * @brief Verifies that MPIMesh::getVertexCoordinates() returns valid
   * coordinates for all local vertices.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetVertexCoordinates)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();
    const size_t sdim = mpiMesh.getSpaceDimension();

    for (Index v = 0; v < shard.getVertexCount(); v++)
    {
      auto coords = mpiMesh.getVertexCoordinates(v);
      EXPECT_EQ(static_cast<size_t>(coords.size()), sdim)
        << "Rank " << world.rank() << " vertex " << v
        << ": coordinate dimension mismatch.";
      // For a UniformGrid({4, 4}) grid, coordinates are in [0, 4]
      for (int j = 0; j < coords.size(); j++)
      {
        EXPECT_GE(coords(j), 0.0 - 1e-14)
          << "Rank " << world.rank() << " vertex " << v << " coord " << j;
        EXPECT_LE(coords(j), 4.0 + 1e-14)
          << "Rank " << world.rank() << " vertex " << v << " coord " << j;
      }
    }
  }

  /**
   * @brief Verifies that a 3D Tetrahedron mesh distributes correctly via
   * MPIMesh and getPolytopeCount returns the correct global count.
   */
  TEST(Rodin_MPI_Geometry_Mesh, Tetrahedron3D_GlobalCellCount)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, {2, 2, 2});

    size_t totalCells = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Tetrahedron, {2, 2, 2});
      totalCells = localMesh.getCellCount();
    }
    boost::mpi::broadcast(world, totalCells, 0);

    const size_t D = mpiMesh.getDimension();
    size_t distributedCount = mpiMesh.getPolytopeCount(D);
    EXPECT_EQ(distributedCount, totalCells)
      << "Rank " << world.rank()
      << ": 3D distributed cell count mismatch.";
  }

  // =========================================================================
  // MPI Mesh — scale()
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::scale() modifies vertex coordinates
   * by the given factor on each rank.
   */
  TEST(Rodin_MPI_Geometry_Mesh, Scale)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    // Record original coordinates for vertex 0
    auto origCoords = mpiMesh.getVertexCoordinates(0);

    const Real factor = 2.5;
    mpiMesh.scale(factor);

    auto newCoords = mpiMesh.getVertexCoordinates(0);
    for (int j = 0; j < origCoords.size(); j++)
    {
      EXPECT_NEAR(newCoords(j), factor * origCoords(j), 1e-12)
        << "Rank " << world.rank() << " coord " << j
        << ": scale factor not applied correctly.";
    }
  }

  // =========================================================================
  // MPI Mesh — flush()
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::flush() does not crash and the mesh
   * remains valid afterward.
   */
  TEST(Rodin_MPI_Geometry_Mesh, Flush)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    // flush should not crash
    mpiMesh.flush();

    // Mesh should still be valid
    EXPECT_EQ(mpiMesh.getDimension(), 2u);
    EXPECT_GT(mpiMesh.getShard().getCellCount(), 0u);
  }

  // =========================================================================
  // MPI Mesh — isSubMesh()
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::isSubMesh() returns false.
   */
  TEST(Rodin_MPI_Geometry_Mesh, IsSubMesh)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    EXPECT_FALSE(mpiMesh.isSubMesh());
  }

  // =========================================================================
  // MPI Mesh — getConnectivity()
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::getConnectivity() returns the shard's
   * connectivity and that cell-to-vertex incidence is valid.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetConnectivity)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& conn = mpiMesh.getConnectivity();
    const size_t D = mpiMesh.getDimension();

    // Cell-to-vertex incidence should exist (was computed before sharding)
    const auto& d2v = conn.getIncidence(D, 0);
    EXPECT_GT(d2v.size(), 0u)
      << "Rank " << world.rank()
      << ": cell-to-vertex incidence is empty.";

    // Each triangle cell should have 3 vertices
    const auto& shard = mpiMesh.getShard();
    for (Index ci = 0; ci < shard.getCellCount(); ci++)
    {
      const auto& verts = conn.getIncidence({D, 0}, ci);
      EXPECT_EQ(verts.size(), 3u)
        << "Rank " << world.rank() << " cell " << ci
        << ": expected 3 vertices for triangle.";
    }
  }

  // =========================================================================
  // MPI Mesh — getAttribute() / setAttribute()
  // =========================================================================

  /**
   * @brief Verifies that getAttribute returns valid attributes and
   * setAttribute modifies them.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetSetAttribute)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const size_t D = mpiMesh.getDimension();
    const auto& shard = mpiMesh.getShard();

    if (shard.getCellCount() > 0)
    {
      // Set attribute on first local cell
      mpiMesh.setAttribute({D, 0}, 42);
      auto attr = mpiMesh.getAttribute(D, 0);
      ASSERT_TRUE(attr.has_value())
        << "Rank " << world.rank()
        << ": getAttribute returned nullopt after setAttribute.";
      EXPECT_EQ(*attr, 42u)
        << "Rank " << world.rank()
        << ": attribute value mismatch after setAttribute.";
    }
  }

  // =========================================================================
  // MPI Mesh — setVertexCoordinates()
  // =========================================================================

  /**
   * @brief Verifies that setVertexCoordinates modifies vertex coordinates
   * and the change is reflected in getVertexCoordinates.
   */
  TEST(Rodin_MPI_Geometry_Mesh, SetVertexCoordinates)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();
    if (shard.getVertexCount() > 0)
    {
      const size_t sdim = mpiMesh.getSpaceDimension();
      Math::SpatialPoint newCoords(static_cast<Eigen::Index>(sdim));
      newCoords(0) = 99.0;
      newCoords(1) = 88.0;
      mpiMesh.setVertexCoordinates(0, newCoords);

      auto readBack = mpiMesh.getVertexCoordinates(0);
      EXPECT_NEAR(readBack(0), 99.0, 1e-14);
      EXPECT_NEAR(readBack(1), 88.0, 1e-14);
    }
  }

  // =========================================================================
  // MPI Mesh — getArea() for 2D triangle mesh
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::getArea() returns the correct total area
   * for a distributed 2D triangle mesh, matching the parent mesh area.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetArea_Triangle2D)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    Real expectedArea = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
      expectedArea = localMesh.getArea();
    }
    boost::mpi::broadcast(world, expectedArea, 0);

    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    Real distributedArea = mpiMesh.getArea();

    EXPECT_NEAR(distributedArea, expectedArea, 1e-10)
      << "Rank " << world.rank()
      << ": distributed area does not match parent area.";
  }

  // =========================================================================
  // MPI Mesh — getVolume() for 3D tetrahedron mesh
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::getVolume() returns the correct total
   * volume for a distributed 3D tetrahedron mesh.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetVolume_Tetrahedron3D)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    Real expectedVolume = 0;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Tetrahedron, {2, 2, 2});
      expectedVolume = localMesh.getVolume();
    }
    boost::mpi::broadcast(world, expectedVolume, 0);

    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, {2, 2, 2});
    Real distributedVolume = mpiMesh.getVolume();

    EXPECT_NEAR(distributedVolume, expectedVolume, 1e-10)
      << "Rank " << world.rank()
      << ": distributed volume does not match parent volume.";
  }

  // =========================================================================
  // MPI Mesh — getPolytopeTransformation()
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::getPolytopeTransformation() returns valid
   * transformations for all local cells.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetPolytopeTransformation)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const size_t D = mpiMesh.getDimension();
    const auto& shard = mpiMesh.getShard();

    for (Index ci = 0; ci < shard.getCellCount(); ci++)
    {
      const auto& trans = mpiMesh.getPolytopeTransformation(D, ci);
      // Transformation order should be non-negative
      EXPECT_GE(trans.getOrder(), 0u)
        << "Rank " << world.rank() << " cell " << ci
        << ": invalid transformation order.";
    }
  }

  // =========================================================================
  // MPI Mesh — getCell(localIdx)
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::getCell(localIdx) returns a valid polytope
   * for each local cell index.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetCellByIndex)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();
    const size_t D = mpiMesh.getDimension();

    for (Index ci = 0; ci < shard.getCellCount(); ci++)
    {
      auto it = mpiMesh.getCell(ci);
      ASSERT_TRUE(static_cast<bool>(it))
        << "Rank " << world.rank() << " cell " << ci
        << ": getCell(idx) returned invalid iterator.";
      EXPECT_EQ(it->getDimension(), D);
      EXPECT_EQ(it->getIndex(), ci);
    }
  }

  // =========================================================================
  // MPI Mesh — getVertex(localIdx)
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::getVertex(localIdx) returns a valid polytope
   * for each local vertex index.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetVertexByIndex)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const auto& shard = mpiMesh.getShard();

    for (Index vi = 0; vi < shard.getVertexCount(); vi++)
    {
      auto it = mpiMesh.getVertex(vi);
      ASSERT_TRUE(static_cast<bool>(it))
        << "Rank " << world.rank() << " vertex " << vi
        << ": getVertex(idx) returned invalid iterator.";
      EXPECT_EQ(it->getDimension(), 0u);
      EXPECT_EQ(it->getIndex(), vi);
    }
  }

  // =========================================================================
  // MPI Mesh — getMeasure(d)
  // =========================================================================

  /**
   * @brief Verifies that MPIMesh::getMeasure(D) returns the same value as
   * getArea() for a 2D mesh.
   */
  TEST(Rodin_MPI_Geometry_Mesh, GetMeasure_MatchesGetArea)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    const size_t D = mpiMesh.getDimension();
    Real area = mpiMesh.getArea();
    Real measure = mpiMesh.getMeasure(D);

    EXPECT_NEAR(measure, area, 1e-10)
      << "Rank " << world.rank()
      << ": getMeasure(D) does not match getArea() for 2D mesh.";
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env = &env;
  g_world = &world;

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
