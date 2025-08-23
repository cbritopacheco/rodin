#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_Real_P1_TestFunction, SanityTest_Build)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);

    EXPECT_EQ(v.Space, TestSpace);
  }

  TEST(Rodin_Variational_Real_P1_TestFunction, FuzzyTest_UniformGrid_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TestFunction v(fes);

    EXPECT_EQ(v.Space, TestSpace);
  }

  TEST(Rodin_Variational_Real_P1_TestFunction, CopyTest)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    
    auto copied = v.copy();
    EXPECT_NE(copied, nullptr);
    EXPECT_EQ(copied->Space, v.Space);
    
    delete copied;
  }

  TEST(Rodin_Variational_Real_P1_TestFunction, CopyConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    TestFunction v_copy(v);

    EXPECT_EQ(v_copy.Space, v.Space);
  }

  TEST(Rodin_Variational_Real_P1_TestFunction, MoveConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    
    TestFunction v_moved(std::move(v));
    EXPECT_EQ(v_moved.Space, TestSpace);
  }

  TEST(Rodin_Variational_Vector_P1_TestFunction, SanityTest_Build)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TestFunction v(fes);

    EXPECT_EQ(v.Space, TestSpace);
  }

  TEST(Rodin_Variational_Vector_P1_TestFunction, ComponentAccess_2D)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TestFunction v(fes);

    auto x_comp = v.x();
    auto y_comp = v.y();

    // Components should have scalar range shape
  }

  TEST(Rodin_Variational_Vector_P1_TestFunction, ComponentAccess_3D)
  {
    constexpr size_t vdim = 3;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TestFunction v(fes);

    auto x_comp = v.x();
    auto y_comp = v.y();
    auto z_comp = v.z();

    // Components should have scalar range shape
  }

  TEST(Rodin_Variational_Real_P1_TestFunction, GetLeaf)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);

    const auto& leaf = v.getLeaf();
    EXPECT_EQ(&leaf, &v);
  }
}