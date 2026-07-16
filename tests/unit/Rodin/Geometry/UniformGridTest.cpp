#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Verifies that invalid grid dimensions raise a member-function exception.
  TEST(Rodin_Geometry_Mesh_UniformGrid, InvalidDimensions)
  {
    EXPECT_THROW(LocalMesh::UniformGrid(Polytope::Type::Triangle, {2}), Alert::Exception);
    EXPECT_THROW(
      LocalMesh::UniformGrid(Polytope::Type::Triangle, {1, 4}), Alert::Exception);
    EXPECT_THROW(
      LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 0, 2}), Alert::Exception);
  }

  /// @brief Verifies sanity test for geometry mesh uniform grid by checking exact expected values.
  TEST(Rodin_Geometry_Mesh_UniformGrid, SanityTest)
  {
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
      EXPECT_EQ(mesh.getVertexCount(), 4);
      EXPECT_EQ(mesh.getCellCount(), 2);
    }

    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 3 });
      EXPECT_EQ(mesh.getVertexCount(), 6);
      EXPECT_EQ(mesh.getCellCount(), 4);
    }

    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 2 });
      EXPECT_EQ(mesh.getVertexCount(), 6);
      EXPECT_EQ(mesh.getCellCount(), 4);
    }

    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 5, 10 });
      EXPECT_EQ(mesh.getVertexCount(), 50);
      EXPECT_EQ(mesh.getCellCount(), 72);
    }

    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 9 });
      EXPECT_EQ(mesh.getVertexCount(), 36);
      EXPECT_EQ(mesh.getCellCount(), 48);
    }

    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 16, 16 });
      EXPECT_EQ(mesh.getVertexCount(), 256);
      EXPECT_EQ(mesh.getCellCount(), 450);
    }
  }

  /// @brief Verifies pyramid one brick for geometry mesh uniform grid by checking exact expected values.
  TEST(Rodin_Geometry_Mesh_UniformGrid, Pyramid_OneBrick)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Pyramid, { 2, 2, 2 });

    EXPECT_EQ(mesh.getVertexCount(), 9);
    EXPECT_EQ(mesh.getCellCount(), 6);

    auto& connectivity = mesh.getConnectivity();
    connectivity.compute(3, 0);
    EXPECT_EQ(connectivity.getCount(3), 6);
    EXPECT_EQ(connectivity.getCount(Polytope::Type::Pyramid), 6);

    connectivity.compute(2, 0);
    EXPECT_EQ(connectivity.getCount(2), 18);
    EXPECT_EQ(connectivity.getCount(Polytope::Type::Quadrilateral), 6);
    EXPECT_EQ(connectivity.getCount(Polytope::Type::Triangle), 12);

    connectivity.compute(1, 0);
    EXPECT_EQ(connectivity.getCount(1), 20);
  }
}
