/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for P0 finite element spaces on SubMesh<Context::MPI>.
 *
 * These tests verify that a distributed P0 FES can be constructed on
 * distributed SubMeshes (both boundary skin and full-cell sub-regions)
 * and that the resulting DOF mappings are consistent across all MPI ranks:
 *
 *   - getSize() matches the expected global DOF count.
 *   - getVectorDimension() reports 1.
 *   - getOwnershipRange() width equals the number of locally owned DOFs.
 *   - getDOFs(D, i) and getGlobalIndex({D, i}, 0) are consistent.
 *   - All local polytopes have a valid (in-range) global DOF.
 *   - Owned DOFs are globally unique across all ranks.
 *
 * NOTE: P1 on distributed SubMesh is not yet fully supported — the
 * distributed P1 constructor can crash at np>=2 when ghost vertices in the
 * SubMesh shard are not present in the owner rank's SubMesh shard.
 * P1 tests on SubMesh<Context::MPI> are excluded until that limitation is
 * resolved.
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
#include <Rodin/MPI/Geometry/SubMesh.h>
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

  static SubMesh<Context::MPI> makeBoundarySubMesh(
      const Mesh<Context::MPI>& mesh)
  {
    const size_t faceDim = mesh.getDimension() - 1;
    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(mesh);
    for (auto it = mesh.getBoundary(); it; ++it)
      builder.include(faceDim, it->getIndex());
    return builder.finalize();
  }

  static SubMesh<Context::MPI> makeCellSubMesh(
      const Mesh<Context::MPI>& mesh)
  {
    const size_t cellDim = mesh.getDimension();
    const auto& shard = mesh.getShard();
    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(mesh);
    for (Index i = 0; i < static_cast<Index>(shard.getCellCount()); ++i)
    {
      if (shard.isOwned(cellDim, i))
        builder.include(cellDim, i);
    }
    return builder.finalize();
  }
}

namespace Rodin::Tests::Unit
{
  // ==========================================================================
  // Group 1 — P0 FES on boundary SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief Scalar P0 FES on boundary SubMesh: global DOF count equals the
   *        number of boundary faces.
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_GetSize_EqualsFaceCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    // P0 DOFs = number of faces (cells of the SubMesh).
    const size_t globalFaces = sub.getPolytopeCount(sub.getDimension());
    EXPECT_EQ(fes.getSize(), globalFaces);
  }

  /**
   * @brief Scalar P0 on boundary SubMesh has vector dimension 1.
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_VectorDimension_IsOne_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);
    EXPECT_EQ(fes.getVectorDimension(), 1u);
  }

  /**
   * @brief Ownership range of P0 on boundary SubMesh has width equal to
   *        the number of locally owned boundary faces.
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_OwnershipRange_MatchesOwnedFaces_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard  = sub.getShard();
    const size_t fDim  = sub.getDimension();
    const size_t nFace = shard.getCellCount();

    size_t ownedCount = 0;
    for (size_t i = 0; i < nFace; ++i)
    {
      if (shard.isOwned(fDim, i))
        ++ownedCount;
    }

    Index begin = 0, end = 0;
    fes.getOwnershipRange(begin, end);
    EXPECT_EQ(static_cast<size_t>(end - begin), ownedCount);
  }

  /**
   * @brief getDOFs and getGlobalIndex are consistent for P0 on boundary SubMesh.
   *
   * For P0 each element has exactly one DOF.  getDOFs(D, i)[0] must equal
   * getGlobalIndex({D, i}, 0).
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_GetDOFs_ConsistentWithGetGlobalIndex_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard  = sub.getShard();
    const size_t fDim  = sub.getDimension();
    const size_t nFace = shard.getCellCount();

    for (size_t i = 0; i < nFace; ++i)
    {
      const auto& dofs = fes.getDOFs(fDim, static_cast<Index>(i));
      ASSERT_EQ(dofs.size(), 1u);
      const Index via_dofs   = dofs[0];
      const Index via_global = fes.getGlobalIndex({ fDim, static_cast<Index>(i) }, 0);
      EXPECT_EQ(via_dofs, via_global)
          << "Mismatch for face " << i;
    }
  }

  /**
   * @brief All local boundary faces have a valid global P0 DOF index.
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_AllLocalFaces_HaveValidGlobalDOF_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard  = sub.getShard();
    const size_t fDim  = sub.getDimension();
    const size_t nFace = shard.getCellCount();

    // Cache getSize() once — it does all_reduce internally and must not be
    // called inside a rank-local loop where iteration counts differ per rank.
    const size_t globalSize = fes.getSize();

    for (size_t i = 0; i < nFace; ++i)
    {
      const auto& dofs = fes.getDOFs(fDim, static_cast<Index>(i));
      for (const Index d : dofs)
      {
        EXPECT_LT(d, globalSize)
            << "Global P0 DOF " << d << " out of range for local face " << i;
      }
    }
  }

  // ==========================================================================
  // Group 2 — P0 FES on full-cell SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief P0 FES on full-cell SubMesh: global DOF count equals the parent
   *        mesh global cell count (SubMesh covers the whole domain).
   */
  TEST(MPIP0FESSubMesh, CellSubMesh_GetSize_EqualsParentCellCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeCellSubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    // Cell SubMesh is built from all owned cells, so global size = global cell count.
    // Collective call — issued by all ranks.
    const size_t parentCells = mesh.getCellCount();
    EXPECT_EQ(fes.getSize(), parentCells);
  }

  /**
   * @brief Owned P0 DOFs on the full-cell SubMesh are globally unique.
   *
   * Gather all owned DOFs on rank 0 and verify no duplicates appear.
   */
  TEST(MPIP0FESSubMesh, CellSubMesh_OwnedDOFs_GloballyUnique_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeCellSubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard     = sub.getShard();
    const size_t cellDim  = sub.getDimension();
    const size_t nCells   = shard.getCellCount();

    std::vector<Index> ownedDofs;
    for (size_t i = 0; i < nCells; ++i)
    {
      if (shard.isOwned(cellDim, i))
      {
        const auto& dofs = fes.getDOFs(cellDim, static_cast<Index>(i));
        ASSERT_EQ(dofs.size(), 1u);
        ownedDofs.push_back(dofs[0]);
      }
    }

    // Collective call — must be issued by all ranks.
    const size_t globalCells = mesh.getCellCount();

    std::vector<std::vector<Index>> allDofs;
    boost::mpi::gather(world, ownedDofs, allDofs, 0);

    if (world.rank() == 0)
    {
      std::vector<Index> combined;
      for (const auto& v : allDofs)
        combined.insert(combined.end(), v.begin(), v.end());

      std::sort(combined.begin(), combined.end());
      const size_t uniqueCount = static_cast<size_t>(
          std::unique(combined.begin(), combined.end()) - combined.begin());

      EXPECT_EQ(uniqueCount, combined.size())
          << "Duplicate global P0 DOF indices on cell SubMesh.";
      EXPECT_EQ(combined.size(), globalCells);
    }
  }

  // ==========================================================================
  // Group 3 — Tetrahedron mesh variant
  // ==========================================================================

  /**
   * @brief P0 on boundary SubMesh of a Tetrahedron mesh: size equals face count.
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_GetSize_EqualsFaceCount_Tetrahedron)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, {3, 3, 3});
    auto sub  = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    const size_t globalFaces = sub.getPolytopeCount(sub.getDimension());
    EXPECT_EQ(fes.getSize(), globalFaces);
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
