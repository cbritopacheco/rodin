#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

TEST(Rodin_Geometry_Mesh_UniformGrid, SanityTest)
{
  {
    Mesh mesh = SequentialMesh::UniformGrid(Polytope::Type::Triangle, 2, 2);
    EXPECT_EQ(mesh.getVertexCount(), 4);
    EXPECT_EQ(mesh.getCellCount(), 2);
  }

  {
    Mesh mesh = SequentialMesh::UniformGrid(Polytope::Type::Triangle, 2, 3);
    EXPECT_EQ(mesh.getVertexCount(), 6);
    EXPECT_EQ(mesh.getCellCount(), 4);
  }

  {
    Mesh mesh = SequentialMesh::UniformGrid(Polytope::Type::Triangle, 3, 2);
    EXPECT_EQ(mesh.getVertexCount(), 6);
    EXPECT_EQ(mesh.getCellCount(), 4);
  }

  {
    Mesh mesh = SequentialMesh::UniformGrid(Polytope::Type::Triangle, 5, 10);
    EXPECT_EQ(mesh.getVertexCount(), 50);
    EXPECT_EQ(mesh.getCellCount(), 72);
  }

  {
    Mesh mesh = SequentialMesh::UniformGrid(Polytope::Type::Triangle, 4, 9);
    EXPECT_EQ(mesh.getVertexCount(), 36);
    EXPECT_EQ(mesh.getCellCount(), 48);
  }

  {
    Mesh mesh = SequentialMesh::UniformGrid(Polytope::Type::Triangle, 16, 16);
    EXPECT_EQ(mesh.getVertexCount(), 256);
    EXPECT_EQ(mesh.getCellCount(), 450);
  }
}
