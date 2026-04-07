/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ==================================================================
  // Mesh::Box tests
  // Box takes a FACE type and produces a mesh of dimension = face_dim+1
  // with only the boundary faces (like a "box" surface)
  // ==================================================================

  TEST(Geometry_MeshBox, Segment_2Nodes)
  {
    // Box with segment faces in 2D space
    auto mesh = LocalMesh::Box(Polytope::Type::Segment, {2, 2});
    EXPECT_EQ(mesh.getSpaceDimension(), 2);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  TEST(Geometry_MeshBox, Segment_4x4Nodes)
  {
    auto mesh = LocalMesh::Box(Polytope::Type::Segment, {4, 4});
    EXPECT_EQ(mesh.getSpaceDimension(), 2);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  TEST(Geometry_MeshBox, Triangle_3Nodes)
  {
    // Box with triangular faces in 3D space
    auto mesh = LocalMesh::Box(Polytope::Type::Triangle, {2, 2, 2});
    EXPECT_EQ(mesh.getSpaceDimension(), 3);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  TEST(Geometry_MeshBox, Quadrilateral_3Nodes)
  {
    // Box with quad faces in 3D space
    auto mesh = LocalMesh::Box(Polytope::Type::Quadrilateral, {2, 2, 2});
    EXPECT_EQ(mesh.getSpaceDimension(), 3);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  TEST(Geometry_MeshBox, Point_2Nodes)
  {
    // Box with point faces in 1D space (just two endpoints)
    auto mesh = LocalMesh::Box(Polytope::Type::Point, {2});
    EXPECT_EQ(mesh.getSpaceDimension(), 1);
    EXPECT_GT(mesh.getVertexCount(), 0);
  }

  TEST(Geometry_MeshBox, Box2D_CoordinatesPositive)
  {
    auto mesh = LocalMesh::Box(Polytope::Type::Segment, {4, 4});
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      auto coords = it->getCoordinates();
      EXPECT_GE(coords(0), -1e-12);
      EXPECT_GE(coords(1), -1e-12);
    }
  }

  TEST(Geometry_MeshBox, Box3D_CoordinatesPositive)
  {
    auto mesh = LocalMesh::Box(Polytope::Type::Triangle, {2, 2, 2});
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      auto coords = it->getCoordinates();
      EXPECT_GE(coords(0), -1e-12);
      EXPECT_GE(coords(1), -1e-12);
      EXPECT_GE(coords(2), -1e-12);
    }
  }

  // ==================================================================
  // Mesh::ccl tests (connected component labeling)
  // ==================================================================

  TEST(Geometry_MeshOps, CCL_SingleComponent)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);

    auto result = mesh.ccl(D, [](const Polytope&, const Polytope&) { return true; });
    const auto& components = result.getComponents();
    // A uniform grid should have exactly one connected component
    EXPECT_EQ(components.size(), 1);
    EXPECT_EQ(components[0].size(), mesh.getCellCount());
  }

  TEST(Geometry_MeshOps, CCL_WithPredicates)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);

    // All cells connected, all included
    auto result = mesh.ccl(
      D,
      [](const Polytope&, const Polytope&) { return true; },
      [](const Polytope&) { return true; }
    );
    EXPECT_EQ(result.getComponents().size(), 1);
  }

  TEST(Geometry_MeshOps, CCL_ExcludeAll)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);

    // Exclude all cells via the inclusion predicate
    auto result = mesh.ccl(
      D,
      [](const Polytope&, const Polytope&) { return true; },
      [](const Polytope&) { return false; }
    );
    EXPECT_EQ(result.getComponents().size(), 0);
  }

  TEST(Geometry_MeshOps, CCL_DisconnectByPredicate)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);

    // Never connect adjacent cells => each cell is its own component
    auto result = mesh.ccl(
      D,
      [](const Polytope&, const Polytope&) { return false; }
    );
    EXPECT_EQ(result.getComponents().size(), mesh.getCellCount());
  }

  // ==================================================================
  // Mesh::trace tests
  // ==================================================================

  TEST(Geometry_MeshOps, Trace_ModifiesMesh)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    // trace requires a Map<Pair<Attribute, Attribute>, Attribute> mapping
    // Get boundary attributes first
    FlatSet<Attribute> bndAttrs;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      auto a = it->getAttribute();
      if (a) bndAttrs.insert(*a);
    }

    // Create a simple trace map
    Map<Pair<Attribute, Attribute>, Attribute> tmap;
    for (auto a : bndAttrs)
      tmap[{a, a}] = a + 100;

    // trace() modifies the mesh in-place
    mesh.trace(tmap);
    // Just verify no crash
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  // ==================================================================
  // Mesh::isSurface tests
  // ==================================================================

  TEST(Geometry_MeshOps, IsSurface_2DMesh)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    // A 2D mesh in 2D space is not a surface mesh
    EXPECT_FALSE(mesh.isSurface());
  }

  TEST(Geometry_MeshOps, IsSurface_Boundary)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);
    auto submesh = mesh.skin();
    // Skin of 2D mesh is 1D embedded in 2D, which IS a surface
    EXPECT_TRUE(submesh.isSurface());
  }

  // ==================================================================
  // Mesh::getPolytope tests
  // ==================================================================

  TEST(Geometry_MeshOps, GetPolytope_Cell)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    EXPECT_FALSE(it.end());
    EXPECT_EQ(it->getDimension(), mesh.getDimension());
  }

  TEST(Geometry_MeshOps, GetPolytope_Vertex)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    auto it = mesh.getPolytope(0, 0);
    EXPECT_FALSE(it.end());
    EXPECT_EQ(it->getDimension(), 0);
  }
}

