/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for the distributed (MPI) P0 finite element space.
 *
 * These tests verify that the MPI P0 space is constructed correctly and that
 * its local/global DOF mappings are consistent:
 *
 *   - getSize() returns the global cell count.
 *   - Every local cell (owned and ghost) has a valid global DOF index.
 *   - Global DOF indices are unique across all ranks.
 *   - Owned DOFs form a contiguous range per rank.
 *   - getDOFs(D, i) returns a single global index for each local cell.
 *   - getGlobalIndex(p, 0) is consistent with getDOFs(D, i).
 *   - For a single-rank run, the distributed space matches the sequential space.
 *
 * Run with: mpirun -n 1/2/3/4 as registered in CMakeLists.txt.
 */

#include <set>
#include <vector>
#include <limits>
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
#include <Rodin/MPI/Variational/P0.h>
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

  static const char* polytopeName(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Tetrahedron: return "Tetrahedron";
      case Polytope::Type::Hexahedron:  return "Hexahedron";
      case Polytope::Type::Pyramid:     return "Pyramid";
      case Polytope::Type::Wedge:       return "Wedge";
      default:                          return "Other";
    }
  }
}

namespace Rodin::Tests::Unit
{
  // =========================================================================
  // P0 MPI: getSize() returns global cell count
  // =========================================================================

  /**
   * @brief getSize() of the distributed P0 space equals the global cell count.
   */
  TEST(MPIP0Space, GetSize_EqualsGlobalCellCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0<Real, Mesh<Context::MPI>> fes(mpiMesh);

    const size_t globalCells = mpiMesh.getCellCount();
    EXPECT_EQ(fes.getSize(), globalCells);
  }

  /**
   * @brief getSize() of the distributed P0 space equals the global cell count
   * for every supported 3D cell mesh.
   */
  TEST(MPIP0Space, GetSize_EqualsGlobalCellCount_All3D)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    for (auto type :
         { Polytope::Type::Tetrahedron,
           Polytope::Type::Hexahedron,
           Polytope::Type::Pyramid,
           Polytope::Type::Wedge })
    {
      SCOPED_TRACE(polytopeName(type));
      auto mpiMesh = distributeFromRoot(ctx, type, { 4, 3, 3 });

      P0<Real, Mesh<Context::MPI>> fes(mpiMesh);

      const size_t globalCells = mpiMesh.getCellCount();
      EXPECT_EQ(fes.getSize(), globalCells);
    }
  }

  // =========================================================================
  // P0 MPI: all local cells have a valid global DOF
  // =========================================================================

  /**
   * @brief Every local cell (owned and ghost) has a valid global DOF index
   * after construction.
   */
  TEST(MPIP0Space, AllLocalCells_HaveValidGlobalDOF_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0<Real, Mesh<Context::MPI>> fes(mpiMesh);

    const auto& shard = mpiMesh.getShard();
    const size_t D    = mpiMesh.getDimension();
    const size_t localCells = shard.getCellCount();

    for (size_t i = 0; i < localCells; ++i)
    {
      const auto& dofs = fes.getDOFs(D, static_cast<Index>(i));
      ASSERT_EQ(dofs.size(), 1u);
      EXPECT_LT(dofs[0], fes.getSize())
          << "Global DOF " << dofs[0]
          << " out of range for cell " << i;
    }
  }

  // =========================================================================
  // P0 MPI: global DOFs are unique across all ranks
  // =========================================================================

  /**
   * @brief Owned global DOF indices are unique across all ranks on a Triangle mesh.
   */
  TEST(MPIP0Space, OwnedDOFs_GloballyUnique_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0<Real, Mesh<Context::MPI>> fes(mpiMesh);

    const auto& shard = mpiMesh.getShard();
    const size_t D    = mpiMesh.getDimension();
    const size_t localCells = shard.getCellCount();

    // Collect owned global DOFs on this rank.
    std::vector<Index> ownedDofs;
    ownedDofs.reserve(localCells);
    for (size_t i = 0; i < localCells; ++i)
    {
      if (shard.isOwned(D, i))
      {
        const auto& dofs = fes.getDOFs(D, static_cast<Index>(i));
        ASSERT_EQ(dofs.size(), 1u);
        ownedDofs.push_back(dofs[0]);
      }
    }

    // Gather all owned DOFs to rank 0 and check uniqueness.
    // getCellCount() calls all_reduce so it must be called by ALL ranks before
    // any rank-gated conditional below.
    const size_t globalCells = mpiMesh.getCellCount();

    std::vector<std::vector<Index>> allDofs;
    boost::mpi::gather(world, ownedDofs, allDofs, 0);

    if (world.rank() == 0)
    {
      std::vector<Index> combined;
      for (const auto& v : allDofs)
        combined.insert(combined.end(), v.begin(), v.end());

      std::sort(combined.begin(), combined.end());
      const size_t uniqueCount =
          static_cast<size_t>(std::unique(combined.begin(), combined.end()) - combined.begin());

      EXPECT_EQ(uniqueCount, combined.size())
          << "Duplicate global DOF indices detected across MPI ranks.";

      // Verify total count matches global cell count.
      EXPECT_EQ(combined.size(), globalCells);
    }
  }

  // =========================================================================
  // P0 MPI: ownership range is contiguous and consistent
  // =========================================================================

  /**
   * @brief getOwnershipRange returns a contiguous range consistent with the
   * number of owned cells.
   */
  TEST(MPIP0Space, OwnershipRange_Contiguous_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0<Real, Mesh<Context::MPI>> fes(mpiMesh);

    const auto& shard = mpiMesh.getShard();
    const size_t D    = mpiMesh.getDimension();
    const size_t localCells = shard.getCellCount();

    // Count owned cells locally.
    size_t ownedCount = 0;
    for (size_t i = 0; i < localCells; ++i)
    {
      if (shard.isOwned(D, i))
        ++ownedCount;
    }

    Index begin = 0, end = 0;
    fes.getOwnershipRange(begin, end);

    EXPECT_EQ(static_cast<size_t>(end - begin), ownedCount)
        << "Ownership range width must equal owned cell count.";
    EXPECT_LE(begin, end);
  }

  // =========================================================================
  // P0 MPI: getDOFs and getGlobalIndex are consistent
  // =========================================================================

  /**
   * @brief getDOFs(D, i)[0] == getGlobalIndex({D, i}, 0) for all local cells.
   */
  TEST(MPIP0Space, GetDOFs_ConsistentWithGetGlobalIndex_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0<Real, Mesh<Context::MPI>> fes(mpiMesh);

    const auto& shard = mpiMesh.getShard();
    const size_t D    = mpiMesh.getDimension();
    const size_t localCells = shard.getCellCount();

    for (size_t i = 0; i < localCells; ++i)
    {
      const auto& dofs = fes.getDOFs(D, static_cast<Index>(i));
      ASSERT_EQ(dofs.size(), 1u);

      const Index via_ggi =
          fes.getGlobalIndex({ D, static_cast<Index>(i) }, 0);

      EXPECT_EQ(dofs[0], via_ggi)
          << "getDOFs and getGlobalIndex disagree for local cell " << i;
    }
  }

  // =========================================================================
  // P0 MPI: single-rank run matches sequential P0 space
  // =========================================================================

  /**
   * @brief For a single MPI rank the distributed P0 space matches the
   * sequential P0 space: same global size and same DOF indices for all cells.
   */
  TEST(MPIP0Space, SingleRank_MatchesSequentialSpace_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() != 1)
      GTEST_SKIP() << "This test is designed for exactly 1 MPI rank.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0<Real, Mesh<Context::MPI>> mpiFes(mpiMesh);

    // Sequential reference.
    auto localMesh = makeShardableMesh(Polytope::Type::Triangle, { 4, 4 });
    P0 seqFes(localMesh);

    EXPECT_EQ(mpiFes.getSize(), seqFes.getSize());

    const size_t D = mpiMesh.getDimension();
    const size_t globalCells = mpiMesh.getCellCount();

    // For a 1-rank run every local cell is owned; local index == global cell ID.
    // The sequential space also uses cell index as the DOF index.
    for (size_t i = 0; i < globalCells; ++i)
    {
      const auto& mpiDofs = mpiFes.getDOFs(D, static_cast<Index>(i));
      const auto& seqDofs = seqFes.getDOFs(D, static_cast<Index>(i));
      ASSERT_EQ(mpiDofs.size(), 1u);
      ASSERT_EQ(seqDofs.size(), 1u);
      EXPECT_EQ(mpiDofs[0], seqDofs[0])
          << "DOF mismatch for cell " << i;
    }
  }

  // =========================================================================
  // P0 MPI: getVectorDimension
  // =========================================================================

  /**
   * @brief Scalar P0 has vector dimension 1.
   */
  TEST(MPIP0Space, VectorDimension_IsOne_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P0<Real, Mesh<Context::MPI>> fes(mpiMesh);
    EXPECT_EQ(fes.getVectorDimension(), 1u);
  }

  // =========================================================================
  // P0 MPI: Quadrilateral mesh
  // =========================================================================

  /**
   * @brief Distributed P0 on a Quadrilateral mesh: global size and unique DOFs.
   */
  TEST(MPIP0Space, QuadrilateralMesh_GlobalSizeAndUniqueDOFs)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Quadrilateral, { 4, 4 });

    P0<Real, Mesh<Context::MPI>> fes(mpiMesh);

    // Global size.
    EXPECT_EQ(fes.getSize(), mpiMesh.getCellCount());

    // Collect owned DOFs.
    const auto& shard = mpiMesh.getShard();
    const size_t D    = mpiMesh.getDimension();
    const size_t localCells = shard.getCellCount();

    std::vector<Index> ownedDofs;
    for (size_t i = 0; i < localCells; ++i)
    {
      if (shard.isOwned(D, i))
      {
        const auto& dofs = fes.getDOFs(D, static_cast<Index>(i));
        ASSERT_EQ(dofs.size(), 1u);
        ownedDofs.push_back(dofs[0]);
      }
    }

    // getCellCount() calls all_reduce so it must be called by ALL ranks before
    // any rank-gated conditional below.
    const size_t globalCells = mpiMesh.getCellCount();

    std::vector<std::vector<Index>> allDofs;
    boost::mpi::gather(world, ownedDofs, allDofs, 0);

    if (world.rank() == 0)
    {
      std::vector<Index> combined;
      for (const auto& v : allDofs)
        combined.insert(combined.end(), v.begin(), v.end());

      std::sort(combined.begin(), combined.end());
      const size_t uniqueCount =
          static_cast<size_t>(std::unique(combined.begin(), combined.end()) - combined.begin());

      EXPECT_EQ(uniqueCount, combined.size());
      EXPECT_EQ(combined.size(), globalCells);
    }
  }
}

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
