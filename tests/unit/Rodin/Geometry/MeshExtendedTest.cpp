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
  // ---- UniformGrid for all geometry types ----

  TEST(Geometry_UniformGrid, Quadrilateral_2x2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {2, 2});
    EXPECT_EQ(mesh.getVertexCount(), 4);
    EXPECT_EQ(mesh.getCellCount(), 1);
    EXPECT_EQ(mesh.getDimension(), 2);
    EXPECT_EQ(mesh.getSpaceDimension(), 2);
  }

  TEST(Geometry_UniformGrid, Quadrilateral_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {4, 4});
    EXPECT_EQ(mesh.getVertexCount(), 16);
    EXPECT_EQ(mesh.getCellCount(), 9);
    EXPECT_EQ(mesh.getDimension(), 2);
  }

  TEST(Geometry_UniformGrid, Quadrilateral_3x5)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {3, 5});
    EXPECT_EQ(mesh.getVertexCount(), 15);
    // Grid of (3-1) x (5-1) = 2 x 4 = 8 quads
    EXPECT_EQ(mesh.getCellCount(), 8);
  }

  TEST(Geometry_UniformGrid, Tetrahedron_2x2x2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    EXPECT_EQ(mesh.getVertexCount(), 8);
    EXPECT_GT(mesh.getCellCount(), 0);
    EXPECT_EQ(mesh.getDimension(), 3);
    EXPECT_EQ(mesh.getSpaceDimension(), 3);
  }

  TEST(Geometry_UniformGrid, Tetrahedron_3x3x3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    EXPECT_EQ(mesh.getVertexCount(), 27);
    EXPECT_GT(mesh.getCellCount(), 0);
    EXPECT_EQ(mesh.getDimension(), 3);
  }

  TEST(Geometry_UniformGrid, Hexahedron_2x2x2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Hexahedron, {2, 2, 2});
    EXPECT_EQ(mesh.getVertexCount(), 8);
    EXPECT_EQ(mesh.getCellCount(), 1);
    EXPECT_EQ(mesh.getDimension(), 3);
    EXPECT_EQ(mesh.getSpaceDimension(), 3);
  }

  TEST(Geometry_UniformGrid, Hexahedron_3x3x3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Hexahedron, {3, 3, 3});
    EXPECT_EQ(mesh.getVertexCount(), 27);
    EXPECT_EQ(mesh.getCellCount(), 8);
    EXPECT_EQ(mesh.getDimension(), 3);
  }

  TEST(Geometry_UniformGrid, Hexahedron_4x3x2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Hexahedron, {4, 3, 2});
    EXPECT_EQ(mesh.getVertexCount(), 24);
    // (4-1) x (3-1) x (2-1) = 3 x 2 x 1 = 6 hexahedra
    EXPECT_EQ(mesh.getCellCount(), 6);
  }

  // ---- Geometric measures ----

  TEST(Geometry_MeshMeasures, TriangleUnitSquare_Area)
  {
    // 16x16 uniform grid spans [0,15]^2, area = 15^2 = 225
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {16, 16});
    Real area = mesh.getArea();
    EXPECT_NEAR(area, 225.0, 1e-10);
  }

  TEST(Geometry_MeshMeasures, TriangleUnitSquare_Perimeter)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {16, 16});
    mesh.getConnectivity().compute(1, 2);
    Real perim = mesh.getPerimeter();
    // Perimeter of [0,15]^2 = 4 * 15 = 60
    EXPECT_NEAR(perim, 60.0, 1e-10);
  }

  TEST(Geometry_MeshMeasures, QuadrilateralUnitSquare_Area)
  {
    // 8x8 uniform grid spans [0,7]^2, area = 49
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {8, 8});
    Real area = mesh.getArea();
    EXPECT_NEAR(area, 49.0, 1e-10);
  }

  TEST(Geometry_MeshMeasures, GetMeasure_Dimension2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {16, 16});
    Real measure = mesh.getMeasure(2);
    EXPECT_NEAR(measure, 225.0, 1e-10);
  }

  TEST(Geometry_MeshMeasures, GetMeasure_Dimension1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {16, 16});
    mesh.getConnectivity().compute(1, 2);
    Real edgeMeasure = mesh.getMeasure(1);
    // Total length of all edges in the mesh; should be positive
    EXPECT_GT(edgeMeasure, 0.0);
  }

  TEST(Geometry_MeshMeasures, ScaledMesh_Area)
  {
    // 8x8 grid has area = 7^2 = 49; scaled by 2 -> area * 4 = 196
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {8, 8});
    Real origArea = mesh.getArea();
    mesh.scale(2.0);
    Real scaledArea = mesh.getArea();
    EXPECT_NEAR(scaledArea, origArea * 4.0, 1e-10);
  }

  TEST(Geometry_MeshMeasures, TetrahedronUnitCube_Volume)
  {
    // 3x3x3 grid spans [0,2]^3, volume = 8
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    Real vol = mesh.getVolume();
    EXPECT_NEAR(vol, 8.0, 1e-10);
  }

  TEST(Geometry_MeshMeasures, HexahedronUnitCube_Volume)
  {
    // 3x3x3 grid spans [0,2]^3, volume = 8
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Hexahedron, {3, 3, 3});
    Real vol = mesh.getVolume();
    EXPECT_NEAR(vol, 8.0, 1e-10);
  }

  TEST(Geometry_MeshMeasures, AreaByAttribute)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {8, 8});
    Real totalArea = mesh.getArea();
    // Set half the cells to attribute 2
    size_t halfCells = mesh.getCellCount() / 2;
    for (size_t i = 0; i < halfCells; ++i)
      mesh.setAttribute({mesh.getDimension(), static_cast<Index>(i)}, 2);

    Real area2 = mesh.getArea(static_cast<Attribute>(2));
    EXPECT_GT(area2, 0.0);
    EXPECT_LT(area2, totalArea);
  }

  // ---- Scale ----

  TEST(Geometry_MeshOps, Scale)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    // Get original last vertex coordinates
    auto origCoords = mesh.getVertexCoordinates(mesh.getVertexCount() - 1);
    Real origX = origCoords(0);
    Real origY = origCoords(1);

    mesh.scale(3.0);
    EXPECT_EQ(mesh.getSpaceDimension(), 2);
    auto coords = mesh.getVertexCoordinates(mesh.getVertexCount() - 1);
    EXPECT_NEAR(coords(0), origX * 3.0, 1e-10);
    EXPECT_NEAR(coords(1), origY * 3.0, 1e-10);
  }

  // ---- Copy / Move ----

  TEST(Geometry_MeshOps, CopyConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    Mesh copy(mesh);
    EXPECT_EQ(copy.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(copy.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(copy.getDimension(), mesh.getDimension());
    EXPECT_EQ(copy.getSpaceDimension(), mesh.getSpaceDimension());
  }

  TEST(Geometry_MeshOps, MoveConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    size_t origVerts = mesh.getVertexCount();
    size_t origCells = mesh.getCellCount();
    Mesh moved(std::move(mesh));
    EXPECT_EQ(moved.getVertexCount(), origVerts);
    EXPECT_EQ(moved.getCellCount(), origCells);
  }

  TEST(Geometry_MeshOps, MoveAssignment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    size_t origVerts = mesh.getVertexCount();
    size_t origCells = mesh.getCellCount();
    Mesh other;
    other = std::move(mesh);
    EXPECT_EQ(other.getVertexCount(), origVerts);
    EXPECT_EQ(other.getCellCount(), origCells);
  }

  // ---- isEmpty / isSurface ----

  TEST(Geometry_MeshOps, IsEmpty_Default)
  {
    Mesh mesh;
    EXPECT_TRUE(mesh.isEmpty());
  }

  TEST(Geometry_MeshOps, IsEmpty_AfterBuild)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    EXPECT_FALSE(mesh.isEmpty());
  }

  // ---- Name ----

  TEST(Geometry_MeshOps, NameSetGet)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    EXPECT_FALSE(mesh.getName().has_value());
    mesh.setName("TestMesh");
    EXPECT_TRUE(mesh.getName().has_value());
    EXPECT_EQ(*mesh.getName(), "TestMesh");
  }

  // ---- setAttribute / getAttribute ----

  TEST(Geometry_MeshOps, SetAttribute)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.setAttribute({mesh.getDimension(), 0}, 42);
    auto attr = mesh.getAttribute(mesh.getDimension(), 0);
    EXPECT_TRUE(attr.has_value());
    EXPECT_EQ(*attr, 42);
  }

  // ---- setVertexCoordinates ----

  TEST(Geometry_MeshOps, SetVertexCoordinates_SingleComponent)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.setVertexCoordinates(0, 99.0, 0);
    auto coords = mesh.getVertexCoordinates(0);
    EXPECT_NEAR(coords(0), 99.0, 1e-14);
  }

  TEST(Geometry_MeshOps, SetVertexCoordinates_Vector)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    Math::SpatialPoint p(2);
    p(0) = 10.0;
    p(1) = 20.0;
    mesh.setVertexCoordinates(0, p);
    auto coords = mesh.getVertexCoordinates(0);
    EXPECT_NEAR(coords(0), 10.0, 1e-14);
    EXPECT_NEAR(coords(1), 20.0, 1e-14);
  }

  // ---- Flush ----

  TEST(Geometry_MeshOps, FlushClearsCache)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    // Access transformation to populate cache
    mesh.getPolytopeTransformation(mesh.getDimension(), 0);
    // Flush should not crash
    mesh.flush();
    // Should still be able to get the transformation after flush
    const auto& trans = mesh.getPolytopeTransformation(mesh.getDimension(), 0);
    (void)trans;  // just verify it doesn't crash
  }
}
