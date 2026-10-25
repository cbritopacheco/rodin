/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <limits>
#include <string>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Checked dimension queries preserve the realized submesh maps.
  TEST(Geometry_SubMesh, CheckedMapDimensionUsesMemberContext)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);
    SubMesh<Context::Local>::Builder builder;
    builder.initialize(mesh);
    for (auto it = mesh.getBoundary(); it; ++it)
      builder.include(1, it->getIndex());
    const auto boundary = builder.finalize();
    const auto& map = boundary.getPolytopeMap(1);
    const auto left = map.left;
    const auto right = map.right;
    EXPECT_TRUE(boundary.getPolytopeMap(mesh.getDimension()).left.empty());
    for (size_t dimension : {mesh.getDimension() + 1, std::numeric_limits<size_t>::max()})
    {
      try
      {
        boundary.getPolytopeMap(dimension);
        ADD_FAILURE() << "An unavailable map dimension must be rejected.";
      }
      catch (const Alert::Exception& error)
      {
        EXPECT_NE(std::string(error.what()).find("getPolytopeMap"), std::string::npos);
        EXPECT_NE(
          std::string(error.what()).find("In member function"), std::string::npos);
      }
      catch (...)
      {
        ADD_FAILURE() << "Checked map access must raise an Alert exception.";
      }
    }
    EXPECT_EQ(&boundary.getPolytopeMap(1), &map);
    EXPECT_EQ(map.left, left);
    EXPECT_EQ(map.right, right);
    EXPECT_EQ(&boundary.getParent(), &mesh);
  }

  // ---- Builder basic construction ----

  /// @brief Verifies build from boundary faces for geometry sub mesh by checking true predicates.
  TEST(Geometry_SubMesh, BuildFromBoundaryFaces)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    SubMesh<Context::Local>::Builder builder;
    builder.initialize(mesh);
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      builder.include(mesh.getDimension() - 1, it->getIndex());

    SubMesh<Context::Local> boundary = builder.finalize();

    // Boundary of a 4x4 triangle grid should have edges on all 4 sides
    EXPECT_GT(boundary.getCellCount(), 0);
    EXPECT_TRUE(boundary.isSubMesh());
  }

  /// @brief Verifies get parent for geometry sub mesh by checking exact expected values.
  TEST(Geometry_SubMesh, GetParent)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(1, 2);

    SubMesh<Context::Local>::Builder builder;
    builder.initialize(mesh);
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      builder.include(mesh.getDimension() - 1, it->getIndex());

    SubMesh<Context::Local> boundary = builder.finalize();

    // Parent should be the original mesh
    EXPECT_EQ(&boundary.getParent(), &mesh);
  }

  /// @brief Verifies get ancestors for geometry sub mesh by checking exact expected values.
  TEST(Geometry_SubMesh, GetAncestors)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(1, 2);

    SubMesh<Context::Local>::Builder builder;
    builder.initialize(mesh);
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      builder.include(mesh.getDimension() - 1, it->getIndex());

    SubMesh<Context::Local> boundary = builder.finalize();

    const auto& ancestors = boundary.getAncestors();
    EXPECT_GE(ancestors.size(), 1);
    // First ancestor should be the parent
    EXPECT_EQ(&ancestors.front().get(), &mesh);
  }

  /// @brief Verifies polytope map for geometry sub mesh by checking exact expected values.
  TEST(Geometry_SubMesh, PolytopeMap)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    SubMesh<Context::Local>::Builder builder;
    builder.initialize(mesh);
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      builder.include(mesh.getDimension() - 1, it->getIndex());

    SubMesh<Context::Local> boundary = builder.finalize();

    // Check that polytopeMap is consistent for cells of the submesh
    const size_t subDim = boundary.getDimension();
    const auto& pmap = boundary.getPolytopeMap(subDim);

    // left[i] gives parent index for submesh cell i
    EXPECT_EQ(pmap.left.size(), boundary.getCellCount());

    // Every left entry should have a right entry
    for (size_t i = 0; i < pmap.left.size(); ++i)
    {
      Index parentIdx = pmap.left[i];
      auto found = pmap.right.find(parentIdx);
      EXPECT_NE(found, pmap.right.end());
      EXPECT_EQ(found->second, static_cast<Index>(i));
    }
  }

  /// @brief Verifies vertex polytope map for geometry sub mesh by checking exact expected values.
  TEST(Geometry_SubMesh, VertexPolytopeMap)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    SubMesh<Context::Local>::Builder builder;
    builder.initialize(mesh);
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      builder.include(mesh.getDimension() - 1, it->getIndex());

    SubMesh<Context::Local> boundary = builder.finalize();

    // Check vertex map consistency
    const auto& vmap = boundary.getPolytopeMap(0);
    EXPECT_EQ(vmap.left.size(), boundary.getVertexCount());

    for (size_t i = 0; i < vmap.left.size(); ++i)
    {
      Index parentIdx = vmap.left[i];
      auto found = vmap.right.find(parentIdx);
      EXPECT_NE(found, vmap.right.end());
      EXPECT_EQ(found->second, static_cast<Index>(i));
    }
  }

  // ---- Mesh operations producing SubMeshes ----

  /// @brief Verifies skin produces submesh for geometry sub mesh by checking exact expected values, true predicates.
  TEST(Geometry_SubMesh, SkinProducesSubmesh)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    SubMesh<Context::Local> skin = mesh.skin();

    EXPECT_TRUE(skin.isSubMesh());
    EXPECT_GT(skin.getCellCount(), 0);
    EXPECT_EQ(&skin.getParent(), &mesh);
    // Boundary of a square domain should form a closed loop
    // For 4x4 grid: boundary has edges on all 4 sides
    EXPECT_GT(skin.getVertexCount(), 0);
  }

  /// @brief Verifies trim single attribute for geometry sub mesh by checking true predicates.
  TEST(Geometry_SubMesh, TrimSingleAttribute)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});

    // Set some cells to attribute 2
    size_t halfCells = mesh.getCellCount() / 2;
    for (size_t i = 0; i < halfCells; ++i)
      mesh.setAttribute({mesh.getDimension(), static_cast<Index>(i)}, 2);

    SubMesh<Context::Local> trimmed = mesh.trim(2);

    EXPECT_TRUE(trimmed.isSubMesh());
    // Trimmed mesh should have fewer cells than original
    EXPECT_LT(trimmed.getCellCount(), mesh.getCellCount());
    EXPECT_GT(trimmed.getCellCount(), 0);
  }

  /// @brief Verifies keep single attribute for geometry sub mesh by checking exact expected values, true predicates.
  TEST(Geometry_SubMesh, KeepSingleAttribute)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});

    // Set some cells to attribute 2
    size_t halfCells = mesh.getCellCount() / 2;
    for (size_t i = 0; i < halfCells; ++i)
      mesh.setAttribute({mesh.getDimension(), static_cast<Index>(i)}, 2);

    SubMesh<Context::Local> kept = mesh.keep(2);

    EXPECT_TRUE(kept.isSubMesh());
    EXPECT_EQ(kept.getCellCount(), halfCells);
  }

  /// @brief Verifies trim and keep are complementary for geometry sub mesh by checking exact expected values.
  TEST(Geometry_SubMesh, TrimAndKeepAreComplementary)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});

    // Mark first half as attribute 2, rest stays default
    size_t halfCells = mesh.getCellCount() / 2;
    for (size_t i = 0; i < halfCells; ++i)
      mesh.setAttribute({mesh.getDimension(), static_cast<Index>(i)}, 2);

    SubMesh<Context::Local> trimmed = mesh.trim(2);
    SubMesh<Context::Local> kept = mesh.keep(2);

    // Together they should account for all cells
    EXPECT_EQ(trimmed.getCellCount() + kept.getCellCount(), mesh.getCellCount());
  }

  // ---- Move semantics ----

  /// @brief Verifies move construction for geometry sub mesh by checking exact expected values, true predicates, move semantics.
  TEST(Geometry_SubMesh, MoveConstruction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(1, 2);

    SubMesh<Context::Local> skin = mesh.skin();
    size_t origCells = skin.getCellCount();
    size_t origVerts = skin.getVertexCount();

    SubMesh<Context::Local> moved(std::move(skin));
    EXPECT_EQ(moved.getCellCount(), origCells);
    EXPECT_EQ(moved.getVertexCount(), origVerts);
    EXPECT_TRUE(moved.isSubMesh());
  }

  /// @brief Verifies is sub mesh flag for geometry sub mesh by checking true predicates, false predicates.
  TEST(Geometry_SubMesh, IsSubMeshFlag)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(1, 2);

    SubMesh<Context::Local> skin = mesh.skin();
    EXPECT_TRUE(skin.isSubMesh());

    // The base Mesh should not be a submesh
    EXPECT_FALSE(mesh.isSubMesh());
  }

  /// @brief Verifies as sub mesh for geometry sub mesh by checking exact expected values, true predicates.
  TEST(Geometry_SubMesh, AsSubMesh)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(1, 2);

    SubMesh<Context::Local> skin = mesh.skin();

    // Can upcast to MeshBase& and still access submesh interface
    MeshBase& meshRef = skin;
    EXPECT_TRUE(meshRef.isSubMesh());

    const SubMeshBase& smBase = meshRef.asSubMesh();
    EXPECT_EQ(&smBase.getParent(), &mesh);
  }
}

/// @brief Reinitialization clears counters and maps even without completing the old build.
TEST(SubMeshBuilder, ReuseStartsFresh)
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  auto parent = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  SubMesh<Context::Local>::Builder builder;
  for (size_t repeat = 0; repeat < 3; ++repeat)
  {
    builder.initialize(parent).include(2, 0);
    if (repeat == 1)
      builder.initialize(parent).include(2, 1);
    const auto submesh = builder.finalize();
    EXPECT_EQ(submesh.getCellCount(), 1u);
    EXPECT_EQ(submesh.getVertexCount(), 3u);
    const Index selected = repeat == 1 ? 1 : 0;
    const auto& vertices = parent.getPolytope(2, selected)->getVertices();
    for (auto vertex = submesh.getVertex(); vertex; ++vertex)
    {
      bool found = false;
      for (const Index original : vertices)
      {
        found |=
          (vertex->getCoordinates() - parent.getVertexCoordinates(original)).norm() <
          1e-14;
      }
      EXPECT_TRUE(found);
    }
  }
}

/// @brief Move assignment replaces the ancestry used by point restriction.
TEST(Geometry_SubMesh, MoveAssignmentReplacesAncestors)
{
  auto first = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  auto second = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  auto destination =
    SubMesh<Context::Local>::Builder().initialize(first).include(1, 0).finalize();
  auto source =
    SubMesh<Context::Local>::Builder().initialize(second).include(1, 0).finalize();
  destination = std::move(source);
  ASSERT_FALSE(destination.getAncestors().empty());
  EXPECT_EQ(&destination.getParent(), &second);
  EXPECT_EQ(&destination.getAncestors().front().get(), &second);
  const Point point(*second.getCell(), Math::SpatialPoint{0.5});
  EXPECT_TRUE(destination.restriction(point));
  const Point oldPoint(*first.getCell(), Math::SpatialPoint{0.5});
  EXPECT_FALSE(destination.restriction(oldPoint));
}
