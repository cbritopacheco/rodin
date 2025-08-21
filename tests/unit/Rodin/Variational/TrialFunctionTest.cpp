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
  TEST(Rodin_Variational_Real_P1_TrialFunction, SanityTest_Build)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);

    EXPECT_EQ(u.Space, TrialSpace);
  }

  TEST(Rodin_Variational_Real_P1_TrialFunction, FuzzyTest_UniformGrid_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TrialFunction u(fes);

    EXPECT_EQ(u.Space, TrialSpace);
  }

  TEST(Rodin_Variational_Real_P1_TrialFunction, CopyConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction u_copy(u);

    EXPECT_EQ(u_copy.Space, u.Space);
  }

  TEST(Rodin_Variational_Real_P1_TrialFunction, MoveConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    
    TrialFunction u_moved(std::move(u));
    EXPECT_EQ(u_moved.Space, TrialSpace);
  }

  TEST(Rodin_Variational_Real_P1_TrialFunction, GetSolution)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);

    const auto& solution = u.getSolution();
    auto& mutable_solution = u.getSolution();
    
  }

  TEST(Rodin_Variational_Vector_P1_TrialFunction, SanityTest_Build)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);

    EXPECT_EQ(u.Space, TrialSpace);
  }

  TEST(Rodin_Variational_Vector_P1_TrialFunction, ComponentAccess_2D)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);

    auto x_comp = u.x();
    auto y_comp = u.y();

    // Components should have scalar range shape
  }

  TEST(Rodin_Variational_Vector_P1_TrialFunction, ComponentAccess_3D)
  {
    constexpr size_t vdim = 3;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);

    auto x_comp = u.x();
    auto y_comp = u.y();
    auto z_comp = u.z();

    // Components should have scalar range shape
  }

  TEST(Rodin_Variational_Real_P1_TrialFunction, GetLeaf)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);

    const auto& leaf = u.getLeaf();
    EXPECT_EQ(&leaf, &u);
  }

  TEST(Rodin_Variational_Vector_P1_TrialFunction, GetSolutionVector)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);

    const auto& solution = u.getSolution();
  }
}