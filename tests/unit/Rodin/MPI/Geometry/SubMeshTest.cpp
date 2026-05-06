/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for SubMesh<Context::MPI>.
 *
 * These tests verify that the distributed submesh infrastructure works correctly:
 * - Building a boundary submesh from a distributed mesh
 * - Building a cell submesh by attribute
 * - Polytope maps are consistent (parent-local ↔ submesh-local)
 * - isSubMesh() returns true
 * - getParent() returns the correct mesh
 * - Global polytope counts match across ranks
 * - restriction() maps points from parent to submesh
 *
 * Run with mpirun -n 1/2/3/4 as registered in CMakeLists.txt.
 */
#include <set>
#include <limits>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/Math/SpatialVector.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/SubMesh.h>

using namespace Rodin;
using namespace Rodin::Geometry;

static boost::mpi::environment*  g_env   = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  Mesh<Context::Local> makeShardableMesh(
      Polytope::Type type, std::initializer_list<size_t> shape)
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
    }
    sharder.scatter(0);
    return sharder.gather(0);
  }
}

// ---------------------------------------------------------------------------

TEST(MPI_Geometry_SubMesh, IsSubMesh)
{
  Context::MPI ctx(*g_env, *g_world);
  auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

  SubMesh<Context::MPI>::Builder builder;
  builder.initialize(mesh);
  for (auto it = mesh.getBoundary(); it; ++it)
    builder.include(mesh.getDimension() - 1, it->getIndex());

  SubMesh<Context::MPI> boundary = builder.finalize();

  EXPECT_TRUE(boundary.isSubMesh());
}

TEST(MPI_Geometry_SubMesh, GetParent)
{
  Context::MPI ctx(*g_env, *g_world);
  auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

  SubMesh<Context::MPI>::Builder builder;
  builder.initialize(mesh);
  for (auto it = mesh.getBoundary(); it; ++it)
    builder.include(mesh.getDimension() - 1, it->getIndex());

  SubMesh<Context::MPI> boundary = builder.finalize();

  EXPECT_EQ(&boundary.getParent(), &mesh);
}

TEST(MPI_Geometry_SubMesh, AncestorChain)
{
  Context::MPI ctx(*g_env, *g_world);
  auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

  SubMesh<Context::MPI>::Builder builder;
  builder.initialize(mesh);
  for (auto it = mesh.getBoundary(); it; ++it)
    builder.include(mesh.getDimension() - 1, it->getIndex());

  SubMesh<Context::MPI> boundary = builder.finalize();

  const auto& ancestors = boundary.getAncestors();
  ASSERT_GE(ancestors.size(), static_cast<size_t>(1));
  EXPECT_EQ(&ancestors.front().get(), &mesh);
}

TEST(MPI_Geometry_SubMesh, BoundaryPolytopeMapConsistency)
{
  Context::MPI ctx(*g_env, *g_world);
  auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
  const size_t faceDim = mesh.getDimension() - 1;

  SubMesh<Context::MPI>::Builder builder;
  builder.initialize(mesh);
  for (auto it = mesh.getBoundary(); it; ++it)
    builder.include(faceDim, it->getIndex());

  SubMesh<Context::MPI> boundary = builder.finalize();

  const auto& pmap = boundary.getPolytopeMap(faceDim);

  // left[i] → parent local index; right[p] → submesh local index
  EXPECT_EQ(pmap.left.size(), pmap.right.size());

  for (size_t i = 0; i < pmap.left.size(); ++i)
  {
    const Index parentLocal = pmap.left[i];
    auto it = pmap.right.find(parentLocal);
    ASSERT_NE(it, pmap.right.end());
    EXPECT_EQ(it->second, static_cast<Index>(i));
  }
}

TEST(MPI_Geometry_SubMesh, VertexPolytopeMapConsistency)
{
  Context::MPI ctx(*g_env, *g_world);
  auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
  const size_t faceDim = mesh.getDimension() - 1;

  SubMesh<Context::MPI>::Builder builder;
  builder.initialize(mesh);
  for (auto it = mesh.getBoundary(); it; ++it)
    builder.include(faceDim, it->getIndex());

  SubMesh<Context::MPI> boundary = builder.finalize();

  const auto& vmap = boundary.getPolytopeMap(0);
  EXPECT_EQ(vmap.left.size(), vmap.right.size());

  for (size_t i = 0; i < vmap.left.size(); ++i)
  {
    const Index parentLocal = vmap.left[i];
    auto it = vmap.right.find(parentLocal);
    ASSERT_NE(it, vmap.right.end());
    EXPECT_EQ(it->second, static_cast<Index>(i));
  }
}

TEST(MPI_Geometry_SubMesh, GlobalBoundaryCellCount)
{
  Context::MPI ctx(*g_env, *g_world);
  auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
  const size_t faceDim = mesh.getDimension() - 1;

  // Compute expected boundary face count from sequential mesh.
  size_t expectedCount = 0;
  if (g_world->rank() == 0)
  {
    auto localMesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
    localMesh.getConnectivity().compute(faceDim, localMesh.getDimension());
    for (auto it = localMesh.getBoundary(); it; ++it)
      ++expectedCount;
  }
  boost::mpi::broadcast(*g_world, expectedCount, 0);

  SubMesh<Context::MPI>::Builder builder;
  builder.initialize(mesh);
  for (auto it = mesh.getBoundary(); it; ++it)
    builder.include(faceDim, it->getIndex());

  SubMesh<Context::MPI> boundary = builder.finalize();

  // getPolytopeCount uses MPI reduction over owned entities.
  const size_t globalCount = boundary.getPolytopeCount(faceDim);
  EXPECT_EQ(globalCount, expectedCount);
}

TEST(MPI_Geometry_SubMesh, SparseBoundarySelection_DimensionConsistentOnEmptyRanks)
{
  const auto& world = *g_world;
  if (world.size() < 2)
    GTEST_SKIP() << "Requires at least two MPI ranks to exercise empty submesh ranks.";

  Context::MPI ctx(*g_env, world);
  auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
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
  boost::mpi::all_gather(world, hasOwnedBoundary ? 1 : 0, hasBoundaryByRank);

  int selectedRank = -1;
  for (int r = 0; r < static_cast<int>(hasBoundaryByRank.size()); ++r)
  {
    if (hasBoundaryByRank[static_cast<size_t>(r)])
    {
      selectedRank = r;
      break;
    }
  }
  ASSERT_GE(selectedRank, 0);

  SubMesh<Context::MPI>::Builder builder;
  builder.initialize(mesh);
  if (world.rank() == selectedRank)
    builder.include(faceDim, selectedFace);

  SubMesh<Context::MPI> sub = builder.finalize();

  EXPECT_EQ(sub.getDimension(), faceDim);
  EXPECT_EQ(sub.getPolytopeCount(faceDim), 1u);

  const size_t localFaces = sub.getShard().getPolytopeCount(faceDim);
  const size_t localEmpty = localFaces == 0 ? 1u : 0u;
  size_t emptyRanks = 0;
  boost::mpi::all_reduce(world, localEmpty, emptyRanks, std::plus<size_t>());
  EXPECT_GT(emptyRanks, 0u);
}

TEST(MPI_Geometry_SubMesh, CellSubmeshByDimension)
{
  Context::MPI ctx(*g_env, *g_world);
  auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
  const size_t cellDim = mesh.getDimension();
  const auto& shard = mesh.getShard();

  SubMesh<Context::MPI>::Builder builder;
  builder.initialize(mesh);

  // Include all locally owned cells.
  for (Index i = 0; i < static_cast<Index>(shard.getCellCount()); ++i)
  {
    if (shard.isOwned(cellDim, i))
      builder.include(cellDim, i);
  }

  SubMesh<Context::MPI> sub = builder.finalize();

  EXPECT_TRUE(sub.isSubMesh());

  // The global cell count of the submesh should equal the global cell count
  // of the parent mesh.
  const size_t parentGlobalCells = mesh.getPolytopeCount(cellDim);
  const size_t subGlobalCells    = sub.getPolytopeCount(cellDim);
  EXPECT_EQ(subGlobalCells, parentGlobalCells);
}

TEST(MPI_Geometry_SubMesh, RestrictionBoundaryPoint)
{
  Context::MPI ctx(*g_env, *g_world);
  auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
  const size_t faceDim = mesh.getDimension() - 1;
  const auto& shard = mesh.getShard();

  SubMesh<Context::MPI>::Builder builder;
  builder.initialize(mesh);
  for (auto it = mesh.getBoundary(); it; ++it)
    builder.include(faceDim, it->getIndex());

  SubMesh<Context::MPI> boundary = builder.finalize();

  // For each owned boundary face of the parent, build a quadrature point
  // and restrict it to the submesh.
  for (auto it = mesh.getBoundary(); it; ++it)
  {
    const Index faceIdx = it->getIndex();
    if (!shard.isOwned(faceDim, faceIdx))
      continue;

    // Reference coordinates at the origin of the reference element.
    Math::SpatialPoint refCoords = Math::SpatialPoint::Zero(faceDim);

    // Build a Point attached to the parent mesh polytope (physical coords
    // are computed lazily from the geometric transformation).
    Point parentPoint(*mesh.getPolytope(faceDim, faceIdx), refCoords);

    Optional<Point> subPoint = boundary.restriction(parentPoint);
    ASSERT_TRUE(subPoint.has_value())
      << "restriction() returned nullopt for owned boundary face " << faceIdx;

    // The restricted point must be in the submesh.
    EXPECT_TRUE(boundary.isLocalPoint(*subPoint));
  }
}

/**
 * @brief Regression test for restriction() with quadrature-produced Points.
 *
 * When Mesh<Context::MPI>::getQuadrature() delegates to the shard, the
 * resulting Points carry the shard (a Mesh<Context::Local>) as their mesh
 * rather than the MPI wrapper.  restriction() must still map such points
 * to the submesh correctly.
 *
 * This test calls polytope.getQuadrature() on a parent MPI mesh polytope and
 * then calls boundary.restriction() on the resulting quadrature points.  It
 * would have returned nullopt before the shard-aware match was added.
 */
TEST(MPI_Geometry_SubMesh, Regression_RestrictionQuadraturePoints_ShardOwnedMesh)
{
  Context::MPI ctx(*g_env, *g_world);
  auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
  const size_t faceDim = mesh.getDimension() - 1;
  const auto& shard = mesh.getShard();

  SubMesh<Context::MPI>::Builder builder;
  builder.initialize(mesh);
  for (auto it = mesh.getBoundary(); it; ++it)
    builder.include(faceDim, it->getIndex());

  SubMesh<Context::MPI> boundary = builder.finalize();

  // Use the default order-1 quadrature formula for faces.
  // Only execute this block if the shard has at least one face.
  if (shard.getPolytopeCount(faceDim) == 0)
    return;

  const auto& qf = QF::PolytopeQuadratureFormula::get(
      1, shard.getGeometry(faceDim, 0));

  // For each owned boundary face, get quadrature points via the mesh's
  // getQuadrature() (which now delegates to the shard) and restrict each
  // point to the boundary submesh.
  for (auto it = mesh.getBoundary(); it; ++it)
  {
    const Index faceIdx = it->getIndex();
    if (!shard.isOwned(faceDim, faceIdx))
      continue;

    // These points have p.getPolytope().getMesh() == shard, not MPIMesh.
    const auto& quadrature = mesh.getQuadrature(faceDim, faceIdx, qf);

    for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
    {
      const auto& p = quadrature.getPoint(qp);

      Optional<Point> subPoint = boundary.restriction(p);
      EXPECT_TRUE(subPoint.has_value())
          << "restriction() returned nullopt for quadrature point " << qp
          << " on boundary face " << faceIdx
          << " (shard-owned-mesh regression).";

      if (subPoint) {
        EXPECT_TRUE(boundary.isLocalPoint(*subPoint));
      }
    }
  }
}

// ---------------------------------------------------------------------------

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env   = &env;
  g_world = &world;

  return RUN_ALL_TESTS();
}
