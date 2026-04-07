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
  // ---- Vertex coordinate access for 2D meshes ----

  TEST(Geometry_VertexExt, XY_2DMesh)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    auto it = mesh.getVertex();
    ASSERT_FALSE(it.end());
    Vertex v = *it;
    // 2D mesh: x() and y() should work
    Real x = v.x();
    Real y = v.y();
    EXPECT_GE(x, 0.0);
    EXPECT_GE(y, 0.0);
  }

  TEST(Geometry_VertexExt, Operator_ParenIndex_2D)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    auto it = mesh.getVertex();
    ASSERT_FALSE(it.end());
    Vertex v = *it;
    Real c0 = v(0);
    Real c1 = v(1);
    EXPECT_EQ(c0, v.x());
    EXPECT_EQ(c1, v.y());
  }

  TEST(Geometry_VertexExt, GetCoordinates_2D)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    auto it = mesh.getVertex();
    ASSERT_FALSE(it.end());
    Vertex v = *it;
    auto coords = v.getCoordinates();
    EXPECT_EQ(coords.size(), 2);
    EXPECT_EQ(coords(0), v.x());
    EXPECT_EQ(coords(1), v.y());
  }

  // ---- Vertex coordinate access for 3D meshes ----

  TEST(Geometry_VertexExt, XYZ_3DMesh)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    auto it = mesh.getVertex();
    ASSERT_FALSE(it.end());
    Vertex v = *it;
    Real x = v.x();
    Real y = v.y();
    Real z = v.z();
    EXPECT_GE(x, 0.0);
    EXPECT_GE(y, 0.0);
    EXPECT_GE(z, 0.0);
  }

  TEST(Geometry_VertexExt, GetCoordinates_3D)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    auto it = mesh.getVertex();
    ASSERT_FALSE(it.end());
    Vertex v = *it;
    auto coords = v.getCoordinates();
    EXPECT_EQ(coords.size(), 3);
    EXPECT_EQ(coords(0), v.x());
    EXPECT_EQ(coords(1), v.y());
    EXPECT_EQ(coords(2), v.z());
  }

  // ---- Multiple vertex traversal ----

  TEST(Geometry_VertexExt, TraverseAll_Triangle)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    size_t count = 0;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Vertex v = *it;
      auto coords = v.getCoordinates();
      EXPECT_EQ(coords.size(), 2);
      ++count;
    }
    EXPECT_EQ(count, mesh.getVertexCount());
  }

  TEST(Geometry_VertexExt, TraverseAll_Tetrahedron)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    size_t count = 0;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Vertex v = *it;
      auto coords = v.getCoordinates();
      EXPECT_EQ(coords.size(), 3);
      ++count;
    }
    EXPECT_EQ(count, mesh.getVertexCount());
  }

  // ---- Face access ----

  TEST(Geometry_FaceExt, BoundaryFace_Triangle)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);
    auto it = mesh.getBoundary();
    ASSERT_FALSE(it.end());
    Face face = *it;
    EXPECT_TRUE(face.isBoundary());
    EXPECT_FALSE(face.isInterface());
    EXPECT_EQ(face.getDimension(), 1);
  }

  TEST(Geometry_FaceExt, InterfaceFace_Triangle)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);
    auto it = mesh.getInterface();
    ASSERT_FALSE(it.end());
    Face face = *it;
    EXPECT_TRUE(face.isInterface());
    EXPECT_FALSE(face.isBoundary());
  }

  // ---- Cell access ----

  TEST(Geometry_CellExt, CellIteration_Triangle)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    size_t count = 0;
    for (auto it = mesh.getCell(); !it.end(); ++it)
    {
      Cell cell = *it;
      EXPECT_EQ(cell.getDimension(), mesh.getDimension());
      EXPECT_TRUE(cell.isCell());
      ++count;
    }
    EXPECT_EQ(count, mesh.getCellCount());
  }

  TEST(Geometry_CellExt, CellIteration_Quadrilateral)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {4, 4});
    size_t count = 0;
    for (auto it = mesh.getCell(); !it.end(); ++it)
    {
      Cell cell = *it;
      EXPECT_EQ(cell.getGeometry(), Polytope::Type::Quadrilateral);
      ++count;
    }
    EXPECT_EQ(count, mesh.getCellCount());
  }

  TEST(Geometry_CellExt, CellAdjacency_Triangle)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(2, 2);
    auto it = mesh.getCell();
    ASSERT_FALSE(it.end());
    Cell cell = *it;
    auto adj = cell.getAdjacent();
    // First cell should have some adjacent cells
    // (just verify no crash, exact count depends on cell position)
    size_t adjCount = 0;
    for (; !adj.end(); ++adj)
      ++adjCount;
    EXPECT_GE(adjCount, 0);
  }

  // ---- Polytope geometry type ----

  TEST(Geometry_PolytopeExt, GetGeometry_AllTypes)
  {
    // Triangle mesh
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
      auto it = mesh.getCell();
      ASSERT_FALSE(it.end());
      EXPECT_EQ(it->getGeometry(), Polytope::Type::Triangle);
    }
    // Quadrilateral mesh
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {3, 3});
      auto it = mesh.getCell();
      ASSERT_FALSE(it.end());
      EXPECT_EQ(it->getGeometry(), Polytope::Type::Quadrilateral);
    }
    // Tetrahedron mesh
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
      auto it = mesh.getCell();
      ASSERT_FALSE(it.end());
      EXPECT_EQ(it->getGeometry(), Polytope::Type::Tetrahedron);
    }
  }

  // ---- Polytope isCell/isFace/isVertex ----

  TEST(Geometry_PolytopeExt, IsCell_IsFace_IsVertex)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    // Cells
    {
      auto it = mesh.getCell();
      ASSERT_FALSE(it.end());
      EXPECT_TRUE(it->isCell());
      EXPECT_FALSE(it->isFace());
      EXPECT_FALSE(it->isVertex());
    }
    // Faces
    {
      auto it = mesh.getBoundary();
      ASSERT_FALSE(it.end());
      EXPECT_FALSE(it->isCell());
      EXPECT_TRUE(it->isFace());
      EXPECT_FALSE(it->isVertex());
    }
    // Vertices
    {
      auto it = mesh.getVertex();
      ASSERT_FALSE(it.end());
      EXPECT_FALSE(it->isCell());
      EXPECT_FALSE(it->isFace());
      EXPECT_TRUE(it->isVertex());
    }
  }

  // ---- Polytope getVertices ----

  TEST(Geometry_PolytopeExt, GetVertices_Triangle)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    auto it = mesh.getCell();
    ASSERT_FALSE(it.end());
    auto& key = it->getVertices();
    EXPECT_EQ(key.size(), 3); // Triangle has 3 vertices
  }

  TEST(Geometry_PolytopeExt, GetVertices_Quadrilateral)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {3, 3});
    auto it = mesh.getCell();
    ASSERT_FALSE(it.end());
    auto& key = it->getVertices();
    EXPECT_EQ(key.size(), 4); // Quad has 4 vertices
  }

  TEST(Geometry_PolytopeExt, GetVertices_Tetrahedron)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    auto it = mesh.getCell();
    ASSERT_FALSE(it.end());
    auto& key = it->getVertices();
    EXPECT_EQ(key.size(), 4); // Tet has 4 vertices
  }
}
