/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for the distributed (MPI) P0g finite element space.
 *
 * These tests verify that the MPI P0g space is constructed correctly and that
 * its global DOF mappings are consistent:
 *
 *   - Scalar P0g has exactly 1 global DOF.
 *   - Vector P0g has exactly vdim global DOFs.
 *   - getDOFs() always returns the same global DOF indices regardless of cell.
 *   - getGlobalIndex() is consistent with getDOFs().
 *   - getOwnershipRange() on rank 0 covers all DOFs; other ranks are empty.
 *   - getSize() does not block (no collective call).
 *   - getVectorDimension() matches construction argument.
 *
 * Run with: mpirun -n 1/2/3/4 as registered in CMakeLists.txt.
 */

#include <vector>
#include <numeric>
#include <algorithm>

#include <gtest/gtest.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Variational/P0g.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// ---------------------------------------------------------------------------
// Global MPI handles (initialized in main())
// ---------------------------------------------------------------------------
static boost::mpi::environment*  g_env   = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  /**
   * @brief Creates a local mesh with all incidences required for sharding.
   */
  static Mesh<Context::Local> makeShardableMesh(
      Polytope::Type type,
      std::initializer_list<size_t> shape)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(type, shape);
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    return mesh;
  }

  /**
   * @brief Distributes a uniform-grid mesh from rank 0.
   */
  static Mesh<Context::MPI> distributeFromRoot(
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
} // anonymous namespace

namespace Rodin::Tests::Unit
{
  // =========================================================================
  // MPIP0gSpace — scalar P0g
  // =========================================================================

  /**
   * @brief Scalar P0g has exactly 1 global DOF on all ranks.
   */
  TEST(MPIP0gSpace, GetSize_IsOne_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0g<Real, Mesh<Context::MPI>> fes(mpiMesh);

    // getSize() is non-collective for P0g — no all_reduce.
    EXPECT_EQ(fes.getSize(), 1u);
  }

  /**
   * @brief Scalar P0g has vector dimension 1.
   */
  TEST(MPIP0gSpace, VectorDimension_IsOne_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0g<Real, Mesh<Context::MPI>> fes(mpiMesh);
    EXPECT_EQ(fes.getVectorDimension(), 1u);
  }

  /**
   * @brief For scalar P0g every local cell maps to global DOF 0.
   */
  TEST(MPIP0gSpace, AllLocalCells_MapToGlobalDOF0_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0g<Real, Mesh<Context::MPI>> fes(mpiMesh);
    const size_t D = mpiMesh.getDimension();
    const auto& shard = mpiMesh.getShard();

    for (size_t i = 0; i < shard.getCellCount(); ++i)
    {
      const auto& dofs = fes.getDOFs(D, static_cast<Index>(i));
      ASSERT_EQ(dofs.size(), 1u);
      EXPECT_EQ(dofs[0], Index(0));
    }
  }

  /**
   * @brief getGlobalIndex() is consistent with getDOFs() for scalar P0g.
   */
  TEST(MPIP0gSpace, GetGlobalIndex_ConsistentWithGetDOFs_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0g<Real, Mesh<Context::MPI>> fes(mpiMesh);
    const size_t D = mpiMesh.getDimension();
    const auto& shard = mpiMesh.getShard();

    for (size_t i = 0; i < shard.getCellCount(); ++i)
    {
      const auto idx = std::make_pair(D, static_cast<Index>(i));
      const auto& dofs = fes.getDOFs(D, static_cast<Index>(i));
      for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
      {
        EXPECT_EQ(fes.getGlobalIndex(idx, local), dofs[static_cast<size_t>(local)]);
      }
    }
  }

  /**
   * @brief Ownership range on rank 0 covers DOF 0; all other ranks are empty.
   */
  TEST(MPIP0gSpace, OwnershipRange_Rank0OwnsAll_OthersEmpty_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0g<Real, Mesh<Context::MPI>> fes(mpiMesh);

    Index begin = 0;
    Index end   = 0;
    fes.getOwnershipRange(begin, end);

    if (world.rank() == 0)
    {
      EXPECT_EQ(begin, Index(0));
      EXPECT_EQ(end,   Index(1));
    }
    else
    {
      EXPECT_EQ(begin, Index(1));
      EXPECT_EQ(end,   Index(1));
    }
  }

  /**
   * @brief For single-rank run, MPI P0g matches the sequential P0g space.
   */
  TEST(MPIP0gSpace, SingleRank_MatchesSequentialSpace_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() != 1)
      GTEST_SKIP() << "This test is designed for exactly 1 MPI rank.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });
    auto localMesh = makeShardableMesh(Polytope::Type::Triangle, { 4, 4 });

    P0g<Real, Mesh<Context::MPI>> mpiFes(mpiMesh);
    P0g<Real, Mesh<Context::Local>> seqFes(localMesh);

    EXPECT_EQ(mpiFes.getSize(), seqFes.getSize());
    EXPECT_EQ(mpiFes.getVectorDimension(), seqFes.getVectorDimension());
  }

  // =========================================================================
  // MPIP0gSpace — vector P0g
  // =========================================================================

  /**
   * @brief Vector P0g(mesh, 2) has exactly 2 global DOFs.
   */
  TEST(MPIP0gSpace, VectorP0g_GetSize_Is2_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0g<Math::SpatialVector<Real>, Mesh<Context::MPI>> fes(mpiMesh, 2);
    EXPECT_EQ(fes.getSize(), 2u);
    EXPECT_EQ(fes.getVectorDimension(), 2u);
  }

  /**
   * @brief Vector P0g: every local cell maps to global DOFs {0, 1, ..., vdim-1}.
   */
  TEST(MPIP0gSpace, VectorP0g_AllLocalCells_MapToGlobalDOFs_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    const size_t vdim = 3;
    P0g<Math::SpatialVector<Real>, Mesh<Context::MPI>> fes(mpiMesh, vdim);
    const size_t D = mpiMesh.getDimension();
    const auto& shard = mpiMesh.getShard();

    for (size_t i = 0; i < shard.getCellCount(); ++i)
    {
      const auto& dofs = fes.getDOFs(D, static_cast<Index>(i));
      ASSERT_EQ(dofs.size(), vdim);
      for (size_t k = 0; k < vdim; ++k)
        EXPECT_EQ(dofs[k], static_cast<Index>(k));
    }
  }

  /**
   * @brief Vector P0g ownership range: rank 0 owns [0, vdim), others empty.
   */
  TEST(MPIP0gSpace, VectorP0g_OwnershipRange_Rank0OwnsAll_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    const size_t vdim = 2;
    P0g<Math::SpatialVector<Real>, Mesh<Context::MPI>> fes(mpiMesh, vdim);

    Index begin = 0;
    Index end   = 0;
    fes.getOwnershipRange(begin, end);

    if (world.rank() == 0)
    {
      EXPECT_EQ(begin, Index(0));
      EXPECT_EQ(end,   static_cast<Index>(vdim));
    }
    else
    {
      EXPECT_EQ(begin, static_cast<Index>(vdim));
      EXPECT_EQ(end,   static_cast<Index>(vdim));
    }
  }

  /**
   * @brief Scalar P0g works on a 3D Tetrahedron mesh.
   */
  TEST(MPIP0gSpace, GetSize_IsOne_Tetrahedron)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, { 3, 3, 3 });

    P0g<Real, Mesh<Context::MPI>> fes(mpiMesh);
    EXPECT_EQ(fes.getSize(), 1u);
  }
} // namespace Rodin::Tests::Unit

// ---------------------------------------------------------------------------
// main() — initializes MPI environment used by all tests.
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env   = &env;
  g_world = &world;

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
