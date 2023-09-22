#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

TEST(Rodin_Geometry_Mesh_UniformGrid, SanityTest)
{
  {
    Mesh mesh = SerialMesh::UniformGrid(Polytope::Type::Triangle, 2, 2);
    EXPECT_EQ(mesh.getVertexCount(), 4);
    EXPECT_EQ(mesh.getElementCount(), 2);
  }

  {
    Mesh mesh = SerialMesh::UniformGrid(Polytope::Type::Triangle, 2, 3);
    EXPECT_EQ(mesh.getVertexCount(), 6);
    EXPECT_EQ(mesh.getElementCount(), 4);
  }

  {
    Mesh mesh = SerialMesh::UniformGrid(Polytope::Type::Triangle, 3, 2);
    EXPECT_EQ(mesh.getVertexCount(), 6);
    EXPECT_EQ(mesh.getElementCount(), 4);
  }

  {
    Mesh mesh = SerialMesh::UniformGrid(Polytope::Type::Triangle, 5, 10);
    EXPECT_EQ(mesh.getVertexCount(), 50);
    EXPECT_EQ(mesh.getElementCount(), 72);
  }

  {
    Mesh mesh = SerialMesh::UniformGrid(Polytope::Type::Triangle, 4, 9);
    EXPECT_EQ(mesh.getVertexCount(), 36);
    EXPECT_EQ(mesh.getElementCount(), 48);
  }

  {
    Mesh mesh = SerialMesh::UniformGrid(Polytope::Type::Triangle, 16, 16);
    EXPECT_EQ(mesh.getVertexCount(), 256);
    EXPECT_EQ(mesh.getElementCount(), 450);
  }
}
