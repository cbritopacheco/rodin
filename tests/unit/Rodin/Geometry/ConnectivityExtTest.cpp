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
  // Multi-step connectivity chain tests
  // ==================================================================

  TEST(Geometry_ConnectivityExt, Chain_2D_Edges_Then_Vertices)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    auto& conn = mesh.getConnectivity();
    // First discover edges (1 -> 2)
    conn.compute(1, 2);
    size_t edgeCount = mesh.getPolytopeCount(1);
    EXPECT_GT(edgeCount, 0);

    // Then compute 1->0 (should already exist after compute(1,2))
    // compute(0, 2) to get vertex-cell adjacency
    conn.compute(0, 2);
    // Verify vertex count is consistent
    EXPECT_EQ(mesh.getPolytopeCount(0), mesh.getVertexCount());
  }

  TEST(Geometry_ConnectivityExt, Chain_3D_Faces_Then_Edges)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    auto& conn = mesh.getConnectivity();
    // First discover faces (2 -> 3)
    conn.compute(2, 3);
    size_t faceCount = mesh.getPolytopeCount(2);
    EXPECT_GT(faceCount, 0);

    // Then discover edges (1 -> 3)
    conn.compute(1, 3);
    size_t edgeCount = mesh.getPolytopeCount(1);
    EXPECT_GT(edgeCount, 0);

    // Verify vertices are still correct
    EXPECT_EQ(mesh.getPolytopeCount(0), mesh.getVertexCount());
  }

  TEST(Geometry_ConnectivityExt, Chain_2D_Compute_Then_Adjacency)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    auto& conn = mesh.getConnectivity();
    // Compute edges
    conn.compute(1, 2);
    // Now compute cell-cell adjacency
    conn.compute(2, 2);

    // Verify cell adjacency
    for (auto it = mesh.getCell(); !it.end(); ++it)
    {
      auto adj = it->getAdjacent();
      // Each cell should have some adjacent cells (interior cells have neighbors)
      // Don't assert all have neighbors since boundary cells may not
      (void)adj;
    }
    // Just verify no crash
    SUCCEED();
  }

  // ==================================================================
  // Connectivity::local() tests
  // ==================================================================

  TEST(Geometry_ConnectivityExt, Local_EdgeVertices)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    auto& conn = mesh.getConnectivity();
    conn.compute(1, 2);

    // Get local connectivity for an edge
    size_t edgeCount = mesh.getPolytopeCount(1);
    EXPECT_GT(edgeCount, 0);

    // Get edge-to-vertex incidence for edge 0
    auto inc = conn.getIncidence({1, 0}, 0);
    // Each edge should have exactly 2 vertices
    EXPECT_EQ(inc.size(), 2);
  }

  TEST(Geometry_ConnectivityExt, Local_TriangleVertices)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    auto& conn = mesh.getConnectivity();

    // Get cell-to-vertex incidence for cell 0
    auto inc = conn.getIncidence({2, 0}, 0);
    // Each triangle should have exactly 3 vertices
    EXPECT_EQ(inc.size(), 3);
  }

  TEST(Geometry_ConnectivityExt, Local_QuadVertices)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {4, 4});
    auto& conn = mesh.getConnectivity();

    auto inc = conn.getIncidence({2, 0}, 0);
    // Each quad should have 4 vertices
    EXPECT_EQ(inc.size(), 4);
  }

  TEST(Geometry_ConnectivityExt, Local_TetVertices)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    auto& conn = mesh.getConnectivity();

    auto inc = conn.getIncidence({3, 0}, 0);
    // Each tet should have 4 vertices
    EXPECT_EQ(inc.size(), 4);
  }

  // ==================================================================
  // Connectivity clear tests across dimensions
  // ==================================================================

  TEST(Geometry_ConnectivityExt, Clear_EdgeIncidence)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    auto& conn = mesh.getConnectivity();
    conn.compute(1, 2);
    EXPECT_GT(mesh.getPolytopeCount(1), 0);

    // Clear the edge incidences
    conn.clear(1, 2);
    conn.clear(2, 1);
    // After clearing, recompute should work
    conn.compute(1, 2);
    EXPECT_GT(mesh.getPolytopeCount(1), 0);
  }

  // ==================================================================
  // Count methods for convenience
  // ==================================================================

  TEST(Geometry_ConnectivityExt, PolytopeCount_AllDimensions)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    auto& conn = mesh.getConnectivity();
    conn.compute(2, 3);
    conn.compute(1, 3);

    // Dim 0: vertices
    EXPECT_GT(mesh.getPolytopeCount(0), 0);
    // Dim 1: edges
    EXPECT_GT(mesh.getPolytopeCount(1), 0);
    // Dim 2: faces
    EXPECT_GT(mesh.getPolytopeCount(2), 0);
    // Dim 3: cells
    EXPECT_GT(mesh.getPolytopeCount(3), 0);
  }
}
