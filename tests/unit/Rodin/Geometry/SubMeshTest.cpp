/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ---- Builder basic construction ----

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

  TEST(Geometry_SubMesh, IsSubMeshFlag)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(1, 2);

    SubMesh<Context::Local> skin = mesh.skin();
    EXPECT_TRUE(skin.isSubMesh());

    // The base Mesh should not be a submesh
    EXPECT_FALSE(mesh.isSubMesh());
  }

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
