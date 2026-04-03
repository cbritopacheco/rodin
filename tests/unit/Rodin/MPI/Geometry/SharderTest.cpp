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
