/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for P0 and P1 finite element spaces on SubMesh<Context::MPI>.
 *
 * These tests verify that distributed P0 and P1 FES can be constructed on
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
 * Run with: mpirun -n 1/2/3/4 as registered in CMakeLists.txt.
 */

#include <set>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>
#include <cassert>

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
#include <Rodin/MPI/Variational/P0g.h>
#include <Rodin/MPI/Variational/P0.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/MPI/Variational/H1.h>
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
    if (D > 1)
    {
      mesh.getConnectivity().compute(D, 1);
      mesh.getConnectivity().compute(D - 1, 1);
    }
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

  static SubMesh<Context::MPI> makeSparseBoundarySubMesh(
      const Mesh<Context::MPI>& mesh,
      const boost::mpi::communicator& comm)
  {
    const size_t faceDim = mesh.getDimension() - 1;
    const auto& shard = mesh.getShard();

    bool hasOwnedBoundary = false;
    Index selectedFace = std::numeric_limits<Index>::max();
    for (auto it = mesh.getBoundary(); it; ++it)
    {
      const Index faceIdx = it->getIndex();
      if (shard.isOwned(faceDim, faceIdx))
      {
        hasOwnedBoundary = true;
        selectedFace = faceIdx;
        break;
      }
    }

    std::vector<int> hasBoundaryByRank;
    boost::mpi::all_gather(comm, hasOwnedBoundary ? 1 : 0, hasBoundaryByRank);

    int selectedRank = -1;
    for (int r = 0; r < static_cast<int>(hasBoundaryByRank.size()); ++r)
    {
      if (hasBoundaryByRank[static_cast<size_t>(r)])
      {
        selectedRank = r;
        break;
      }
    }
    assert(selectedRank >= 0);

    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(mesh);
    if (comm.rank() == selectedRank)
      builder.include(faceDim, selectedFace);
    return builder.finalize();
  }
}

namespace Rodin::Tests::Unit
{
  TEST(MPISubMeshSparseSelection, EmptyRanks_ConstructP0P1H1_TriangleBoundary)
  {
    const auto& world = *g_world;
    if (world.size() < 2)
      GTEST_SKIP() << "Requires at least two MPI ranks to exercise empty submesh ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeSparseBoundarySubMesh(mesh, world);

    const size_t faceDim = mesh.getDimension() - 1;
    ASSERT_EQ(sub.getDimension(), faceDim);
    ASSERT_EQ(sub.getPolytopeCount(faceDim), 1u);

    P0<Real, Mesh<Context::MPI>> p0(sub);
    EXPECT_EQ(p0.getSize(), 1u);

    P1<Real, Mesh<Context::MPI>> p1(sub);
    const size_t globalVertices = sub.getPolytopeCount(0);
    EXPECT_EQ(p1.getSize(), globalVertices);

    H1<2, Real, Mesh<Context::MPI>> h1(std::integral_constant<size_t, 2>{}, sub);
    EXPECT_EQ(h1.getSize(), globalVertices + 1u);

    Index begin = 0;
    Index end = 0;
    h1.getOwnershipRange(begin, end);
    const size_t localOwned = static_cast<size_t>(end - begin);
    size_t globalOwned = 0;
    boost::mpi::all_reduce(world, localOwned, globalOwned, std::plus<size_t>());
    EXPECT_EQ(globalOwned, h1.getSize());
  }

  TEST(MPISubMeshP0g, BoundarySubMesh_ScalarGlobalConstantDOFs_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P0g<Real, Mesh<Context::MPI>> fes(sub);
    EXPECT_EQ(fes.getSize(), 1u);
    EXPECT_EQ(fes.getVectorDimension(), 1u);

    Index begin = 0;
    Index end = 0;
    fes.getOwnershipRange(begin, end);
    if (world.rank() == 0)
    {
      EXPECT_EQ(begin, Index(0));
      EXPECT_EQ(end, Index(1));
    }
    else
    {
      EXPECT_EQ(begin, Index(1));
      EXPECT_EQ(end, Index(1));
    }

    const size_t d = sub.getDimension();
    const auto& shard = sub.getShard();
    for (Index i = 0; i < static_cast<Index>(shard.getPolytopeCount(d)); ++i)
    {
      const auto& dofs = fes.getDOFs(d, i);
      ASSERT_EQ(dofs.size(), 1u);
      EXPECT_EQ(dofs[0], Index(0));
      EXPECT_EQ(fes.getGlobalIndex({d, i}, 0), Index(0));
    }
  }

  TEST(MPISubMeshP0g, CellSubMesh_VectorGlobalConstantDOFs_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeCellSubMesh(mesh);

    const size_t vdim = 3;
    P0g<Math::SpatialVector<Real>, Mesh<Context::MPI>> fes(sub, vdim);
    EXPECT_EQ(fes.getSize(), vdim);
    EXPECT_EQ(fes.getVectorDimension(), vdim);

    Index begin = 0;
    Index end = 0;
    fes.getOwnershipRange(begin, end);
    if (world.rank() == 0)
    {
      EXPECT_EQ(begin, Index(0));
      EXPECT_EQ(end, static_cast<Index>(vdim));
    }
    else
    {
      EXPECT_EQ(begin, static_cast<Index>(vdim));
      EXPECT_EQ(end, static_cast<Index>(vdim));
    }

    const size_t d = sub.getDimension();
    const auto& shard = sub.getShard();
    for (Index i = 0; i < static_cast<Index>(shard.getPolytopeCount(d)); ++i)
    {
      const auto& dofs = fes.getDOFs(d, i);
      ASSERT_EQ(dofs.size(), vdim);
      for (size_t k = 0; k < vdim; ++k)
      {
        EXPECT_EQ(dofs[k], static_cast<Index>(k));
        EXPECT_EQ(fes.getGlobalIndex({d, i}, static_cast<Index>(k)),
                  static_cast<Index>(k));
      }
    }
  }

  TEST(MPISubMeshP0g, SparseBoundarySubMesh_ConstructsOnEmptyRanks)
  {
    const auto& world = *g_world;
    if (world.size() < 2)
      GTEST_SKIP() << "Requires at least two MPI ranks to exercise empty submesh ranks.";
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeSparseBoundarySubMesh(mesh, world);

    const size_t d = sub.getDimension();
    ASSERT_EQ(sub.getPolytopeCount(d), 1u);

    P0g<Real, Mesh<Context::MPI>> fes(sub);
    EXPECT_EQ(fes.getSize(), 1u);

    const auto& shard = sub.getShard();
    for (Index i = 0; i < static_cast<Index>(shard.getPolytopeCount(d)); ++i)
    {
      const auto& dofs = fes.getDOFs(d, i);
      ASSERT_EQ(dofs.size(), 1u);
      EXPECT_EQ(dofs[0], Index(0));
    }
  }

  TEST(MPISubMeshH1, BoundarySubMesh_QuadraticH1_Tetrahedron)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, {2, 2, 2});
    auto sub  = makeBoundarySubMesh(mesh);

    ASSERT_EQ(sub.getDimension(), 2u);

    H1<2, Real, Mesh<Context::MPI>> h1(std::integral_constant<size_t, 2>{}, sub);

    const size_t globalVertices = sub.getPolytopeCount(0);
    const size_t globalEdges = sub.getPolytopeCount(1);
    EXPECT_EQ(h1.getSize(), globalVertices + globalEdges);

    const auto& shard = sub.getShard();
    for (Index i = 0; i < static_cast<Index>(shard.getPolytopeCount(2)); ++i)
    {
      const auto& dofs = h1.getDOFs(2, i);
      EXPECT_EQ(dofs.size(), 6u);
      for (const Index dof : dofs)
        EXPECT_LT(dof, static_cast<Index>(h1.getSize()));
    }
  }

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

  // ==========================================================================
  // Group 4 — P1 FES on boundary SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief Scalar P1 FES on boundary SubMesh: global DOF count equals the
   *        number of boundary vertices.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_GetSize_EqualsVertexCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const size_t globalVerts = sub.getPolytopeCount(0);
    EXPECT_EQ(fes.getSize(), globalVerts);
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

    const auto& shard = sub.getShard();
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

    const auto& shard = sub.getShard();
    const size_t nVert = shard.getVertexCount();

    const size_t globalSize = fes.getSize();

    for (size_t i = 0; i < nVert; ++i)
    {
      const auto& dofs = fes.getDOFs(0, static_cast<Index>(i));
      for (const Index d : dofs)
      {
        EXPECT_LT(d, globalSize)
            << "Global P1 DOF " << d << " out of range for local vertex " << i;
      }
    }
  }

  // ==========================================================================
  // Group 5 — P1 FES on full-cell SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief P1 FES on full-cell SubMesh: global DOF count equals the parent
   *        mesh global vertex count.
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

    const size_t parentVerts = mesh.getVertexCount();
    EXPECT_EQ(fes.getSize(), parentVerts);
  }

  /**
   * @brief Owned P1 DOFs on the full-cell SubMesh are globally unique.
   *
   * Gather all owned DOFs on rank 0 and verify no duplicates appear.
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

    const auto& shard    = sub.getShard();
    const size_t nVerts  = shard.getVertexCount();

    std::vector<Index> ownedDofs;
    for (size_t i = 0; i < nVerts; ++i)
    {
      if (shard.isOwned(0, i))
      {
        const auto& dofs = fes.getDOFs(0, static_cast<Index>(i));
        for (const Index d : dofs)
          ownedDofs.push_back(d);
      }
    }

    const size_t globalVerts = mesh.getVertexCount();

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
  // Group 6 — P1 FES on boundary SubMesh: owned DOF uniqueness
  // ==========================================================================

  /**
   * @brief Owned P1 DOFs on boundary SubMesh are globally unique and cover
   *        the expected range exactly.
   *
   * This is a regression test for the MPI tag-collision bug: before the fix,
   * P0 cell SubMesh sends with tag=0 would accumulate and be consumed by a
   * subsequent P1 boundary SubMesh irecv on the same tag, leaving some
   * Shared vertices with invalid DOFs.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_OwnedDOFs_GloballyUnique_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard   = sub.getShard();
    const size_t nVerts = shard.getVertexCount();

    std::vector<Index> ownedDofs;
    for (size_t i = 0; i < nVerts; ++i)
    {
      if (shard.isOwned(0, i))
      {
        const auto& dofs = fes.getDOFs(0, static_cast<Index>(i));
        for (const Index d : dofs)
          ownedDofs.push_back(d);
      }
    }

    // Collective gather — must be called by all ranks.
    const size_t globalBoundaryVerts = sub.getVertexCount();

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
          << "Duplicate global P1 DOF indices on boundary SubMesh.";
      EXPECT_EQ(combined.size(), globalBoundaryVerts)
          << "Owned P1 DOF count does not match global boundary vertex count.";
    }
  }

  // ==========================================================================
  // Group 7 — Regression: P0 cell SubMesh then P1 boundary SubMesh (tag-leakage)
  // ==========================================================================

  /**
   * @brief Regression test for MPI message tag collision.
   *
   * Previously, P0 on a cell SubMesh used tag=0 for its DOF exchange and
   * could leave unmatched sends in the MPI buffer.  A subsequent P1 boundary
   * SubMesh construction would irecv on tag=0 and consume those stale P0
   * messages, corrupting Shared-vertex DOF assignments.
   *
   * After the fix (symmetric drain, distinct tags), the P1 boundary SubMesh
   * must produce valid DOFs even when constructed after P0 cell and P0
   * boundary SubMesh spaces.
   */
  TEST(MPIP1FESSubMesh,
       Regression_P0CellThenP1Boundary_AllLocalVertices_HaveValidGlobalDOF)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    // Build both SubMeshes.
    auto cellSub     = makeCellSubMesh(mesh);
    auto boundarySub = makeBoundarySubMesh(mesh);

    // Construct P0 on cell SubMesh first (this was the source of stale sends).
    P0<Real, Mesh<Context::MPI>> p0cell(cellSub);
    (void) p0cell;

    // Construct P0 on boundary SubMesh (additional potential source of stale sends).
    P0<Real, Mesh<Context::MPI>> p0bnd(boundarySub);
    (void) p0bnd;

    // Now construct P1 on boundary SubMesh.  After the symmetric-drain fix,
    // this must not consume any stale P0 messages and must produce valid DOFs.
    P1<Real, Mesh<Context::MPI>> p1bnd(boundarySub);

    const auto& shard   = boundarySub.getShard();
    const size_t nVerts = shard.getVertexCount();
    const size_t globalSize = p1bnd.getSize();

    for (size_t i = 0; i < nVerts; ++i)
    {
      const auto& dofs = p1bnd.getDOFs(0, static_cast<Index>(i));
      for (const Index d : dofs)
      {
        EXPECT_LT(d, globalSize)
            << "P1 boundary SubMesh DOF " << d
            << " out of range for local vertex " << i
            << " (constructed after P0 cell and boundary SubMesh spaces).";
      }
    }
  }

  /**
   * @brief Same regression test but with a Tetrahedron mesh (3D variant).
   */
  TEST(MPIP1FESSubMesh,
       Regression_P0CellThenP1Boundary_AllLocalVertices_HaveValidGlobalDOF_Tetrahedron)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, {3, 3, 3});

    auto cellSub     = makeCellSubMesh(mesh);
    auto boundarySub = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> p0cell(cellSub);
    (void) p0cell;

    P0<Real, Mesh<Context::MPI>> p0bnd(boundarySub);
    (void) p0bnd;

    P1<Real, Mesh<Context::MPI>> p1bnd(boundarySub);

    const auto& shard   = boundarySub.getShard();
    const size_t nVerts = shard.getVertexCount();
    const size_t globalSize = p1bnd.getSize();

    for (size_t i = 0; i < nVerts; ++i)
    {
      const auto& dofs = p1bnd.getDOFs(0, static_cast<Index>(i));
      for (const Index d : dofs)
      {
        EXPECT_LT(d, globalSize)
            << "P1 boundary SubMesh DOF " << d
            << " out of range for local vertex " << i;
      }
    }
  }

  // ==========================================================================
  // Group 8 — Regression: SubMesh orphaned vertex (2-round ownership fix)
  // ==========================================================================

  /**
   * @brief Regression test for the SubMesh orphaned-vertex ownership bug.
   *
   * A vertex v with Shared(owner=A) in the volume partition is added to the
   * boundary SubMesh on ranks B (and possibly C,D,...) when those ranks have
   * boundary faces adjacent to v.  If rank A has no boundary face adjacent
   * to v, the SubMesh finalize() 2-round exchange must promote the
   * minimum-rank querier to Owned so that P1 can assign a valid global DOF.
   *
   * Passing this test requires the 2-round ownership fix in SubMesh::finalize.
   */
  TEST(MPIP1FESSubMesh,
       Regression_OrphanVertex_AllLocalVertices_HaveValidGlobalDOF_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    // A finer 8×8 mesh increases the probability of hitting corner vertices
    // where the volume-partition owner is not adjacent to any boundary face.
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {8, 8});
    auto sub  = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard   = sub.getShard();
    const size_t nVerts = shard.getVertexCount();
    const size_t globalSize = fes.getSize();

    for (size_t i = 0; i < nVerts; ++i)
    {
      const auto& dofs = fes.getDOFs(0, static_cast<Index>(i));
      for (const Index d : dofs)
      {
        EXPECT_LT(d, globalSize)
            << "P1 DOF " << d << " out of range for local vertex " << i
            << " (orphan-vertex regression, 8x8 triangle mesh).";
      }
    }

    // Also verify owned DOFs are globally unique.
    std::vector<Index> ownedDofs;
    for (size_t i = 0; i < nVerts; ++i)
    {
      if (shard.isOwned(0, i))
      {
        const auto& dofs = fes.getDOFs(0, static_cast<Index>(i));
        for (const Index d : dofs)
          ownedDofs.push_back(d);
      }
    }

    const size_t globalBoundaryVerts = sub.getVertexCount();

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
          << "Duplicate owned P1 DOFs on 8x8 boundary SubMesh (orphan-vertex regression).";
      EXPECT_EQ(combined.size(), globalBoundaryVerts)
          << "Owned P1 DOF count does not match global boundary vertex count.";
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
