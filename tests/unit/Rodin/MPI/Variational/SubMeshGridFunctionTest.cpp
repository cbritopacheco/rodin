/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for P0 and P1 finite element spaces on SubMesh<Context::MPI>.
 *
 * These tests verify that distributed FES can be constructed on distributed
 * SubMeshes (both boundary skin and cell-attribute sub-regions) and that
 * the resulting DOF mappings are consistent:
 *
 *   - getSize() matches the expected global DOF count for P0 / P1 on
 *     a boundary SubMesh and a cell-attribute SubMesh.
 *   - getVectorDimension() reports the correct dimension.
 *   - getOwnershipRange() width equals the number of locally owned DOFs.
 *   - getDOFs(D, i) and getGlobalIndex({D, i}, 0) are consistent for P0.
 *   - All local cells/vertices have a valid (in-range) global DOF.
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
#include <Rodin/MPI/Variational/P1.h>
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
  // -------------------------------------------------------------------------
  // Helper: build a fully-connected shardable local mesh.
  // -------------------------------------------------------------------------
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

  // -------------------------------------------------------------------------
  // Helper: distribute a uniform-grid mesh from rank 0.
  // -------------------------------------------------------------------------
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
    }
    sharder.scatter(0);
    return sharder.gather(0);
  }

  // -------------------------------------------------------------------------
  // Helper: build a boundary SubMesh from a distributed parent mesh.
  // -------------------------------------------------------------------------
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

  // -------------------------------------------------------------------------
  // Helper: build a cell SubMesh keeping only owned cells.
  // -------------------------------------------------------------------------
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
   * @brief P0 FES on boundary SubMesh: global DOF size equals global face count.
   *
   * For a P0 space on the boundary SubMesh each face is one DOF.
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

    // Global DOF count == number of boundary faces
    const size_t globalFaces = sub.getPolytopeCount(sub.getDimension());
    EXPECT_EQ(fes.getSize(), globalFaces);
  }

  /**
   * @brief Scalar P0 on SubMesh has vector dimension 1.
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
   * @brief Ownership range width equals locally owned face count on SubMesh P0.
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
   * @brief getDOFs and getGlobalIndex are consistent for every local face on SubMesh P0.
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

      const Index viaGGI =
          fes.getGlobalIndex({fDim, static_cast<Index>(i)}, 0);
      EXPECT_EQ(dofs[0], viaGGI)
          << "getDOFs and getGlobalIndex disagree for local face " << i;
    }
  }

  /**
   * @brief All local faces have a valid global DOF index (within [0, size)).
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

    for (size_t i = 0; i < nFace; ++i)
    {
      const auto& dofs = fes.getDOFs(fDim, static_cast<Index>(i));
      ASSERT_EQ(dofs.size(), 1u);
      EXPECT_LT(dofs[0], fes.getSize())
          << "Global DOF " << dofs[0] << " out of range for local face " << i;
    }
  }

  // ==========================================================================
  // Group 2 — P0 FES on full-cell SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief P0 FES on cell SubMesh: global DOF count equals the parent mesh
   *        global cell count.
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

    // getCellCount() is collective — call it outside any rank-gated conditional.
    const size_t parentCells = mesh.getCellCount();
    EXPECT_EQ(fes.getSize(), parentCells);
  }

  /**
   * @brief Owned DOFs on cell SubMesh P0 are globally unique and cover the
   *        entire index range.
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
    const size_t nLocal   = shard.getCellCount();

    std::vector<Index> ownedDofs;
    for (size_t i = 0; i < nLocal; ++i)
    {
      if (shard.isOwned(cellDim, i))
      {
        const auto& dofs = fes.getDOFs(cellDim, static_cast<Index>(i));
        ASSERT_EQ(dofs.size(), 1u);
        ownedDofs.push_back(dofs[0]);
      }
    }

    // Collective: must be called by all ranks before the rank-gated conditional.
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
          << "Duplicate global DOF indices on cell SubMesh P0.";
      EXPECT_EQ(combined.size(), globalCells);
    }
  }

  // ==========================================================================
  // Group 3 — P1 FES on boundary SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief Scalar P1 FES on boundary SubMesh: global DOF count equals the
   *        number of boundary vertices.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_GetSize_EqualsBoundaryVertexCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    // P1 DOFs = number of vertices in the boundary SubMesh
    const size_t globalVerts = sub.getPolytopeCount(0);
    EXPECT_EQ(fes.getSize(), globalVerts);
  }

  /**
   * @brief Scalar P1 on SubMesh has vector dimension 1.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_VectorDimension_IsOne_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);
    EXPECT_EQ(fes.getVectorDimension(), 1u);
  }

  /**
   * @brief Ownership range of P1 on boundary SubMesh has width equal to
   *        the number of locally owned boundary vertices.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_OwnershipRange_MatchesOwnedVertices_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard  = sub.getShard();
    const size_t nVert = shard.getVertexCount();

    size_t ownedCount = 0;
    for (size_t i = 0; i < nVert; ++i)
    {
      if (shard.isOwned(0, i))
        ++ownedCount;
    }

    Index begin = 0, end = 0;
    fes.getOwnershipRange(begin, end);
    EXPECT_EQ(static_cast<size_t>(end - begin), ownedCount);
  }

  /**
   * @brief All local boundary vertices have a valid global P1 DOF index.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_AllLocalVertices_HaveValidGlobalDOF_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard  = sub.getShard();
    const size_t fDim  = sub.getDimension();
    const size_t nFace = shard.getCellCount();

    for (size_t i = 0; i < nFace; ++i)
    {
      const auto& dofs = fes.getDOFs(fDim, static_cast<Index>(i));
      for (const Index d : dofs)
      {
        EXPECT_LT(d, fes.getSize())
            << "Global P1 DOF " << d << " out of range for local face " << i;
      }
    }
  }

  // ==========================================================================
  // Group 4 — P1 FES on full-cell SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief P1 FES on cell SubMesh: global DOF count equals the parent mesh
   *        global vertex count (same vertices, same P1 DOFs).
   */
  TEST(MPIP1FESSubMesh, CellSubMesh_GetSize_EqualsParentVertexCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeCellSubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const size_t parentVerts = mesh.getPolytopeCount(0);
    EXPECT_EQ(fes.getSize(), parentVerts);
  }

  /**
   * @brief Vector-valued P1 (vdim=2) on cell SubMesh: global DOF count
   *        equals 2 × parent vertex count.
   */
  TEST(MPIP1FESSubMesh, CellSubMesh_VectorP1_GetSize_EqualsTwoTimesVertexCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    constexpr size_t vdim = 2;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeCellSubMesh(mesh);

    P1<Math::SpatialVector<Real>, Mesh<Context::MPI>> fes(sub, vdim);

    const size_t parentVerts = mesh.getPolytopeCount(0);
    EXPECT_EQ(fes.getSize(), vdim * parentVerts);
    EXPECT_EQ(fes.getVectorDimension(), vdim);
  }

  /**
   * @brief Owned P1 DOFs on cell SubMesh are globally unique.
   */
  TEST(MPIP1FESSubMesh, CellSubMesh_OwnedDOFs_GloballyUnique_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeCellSubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard     = sub.getShard();
    const size_t nVert    = shard.getVertexCount();

    std::vector<Index> ownedDofs;
    for (size_t v = 0; v < nVert; ++v)
    {
      if (shard.isOwned(0, v))
      {
        const auto& dofs = fes.getDOFs(0, static_cast<Index>(v));
        ASSERT_EQ(dofs.size(), 1u);
        ownedDofs.push_back(dofs[0]);
      }
    }

    // Collective call — must be issued by all ranks.
    const size_t globalVerts = mesh.getPolytopeCount(0);

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
          << "Duplicate global P1 DOF indices on cell SubMesh.";
      EXPECT_EQ(combined.size(), globalVerts);
    }
  }

  // ==========================================================================
  // Group 5 — Tetrahedron mesh variants
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

  /**
   * @brief P1 on boundary SubMesh of a Tetrahedron mesh: size equals
   *        boundary vertex count.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_GetSize_EqualsBoundaryVertexCount_Tetrahedron)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, {3, 3, 3});
    auto sub  = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const size_t globalVerts = sub.getPolytopeCount(0);
    EXPECT_EQ(fes.getSize(), globalVerts);
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
