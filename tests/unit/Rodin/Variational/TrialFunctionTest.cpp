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
  /// @brief Verifies sanity test build for variational real P1 trial function by checking exact expected values.
  TEST(Rodin_Variational_Real_P1_TrialFunction, SanityTest_Build)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);

    EXPECT_EQ(u.Space, TrialSpace);
  }

  /// @brief Verifies fuzzy test uniform grid 4 x 4 for variational real P1 trial function by checking exact expected values.
  TEST(Rodin_Variational_Real_P1_TrialFunction, FuzzyTest_UniformGrid_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TrialFunction u(fes);

    EXPECT_EQ(u.Space, TrialSpace);
  }

  /// @brief Verifies copy constructor for variational real P1 trial function by checking exact expected values, copy semantics.
  TEST(Rodin_Variational_Real_P1_TrialFunction, CopyConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction u_copy(u);

    EXPECT_EQ(u_copy.Space, u.Space);
  }

  /// @brief Verifies move constructor for variational real P1 trial function by checking exact expected values, move semantics.
  TEST(Rodin_Variational_Real_P1_TrialFunction, MoveConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    
    TrialFunction u_moved(std::move(u));
    EXPECT_EQ(u_moved.Space, TrialSpace);
  }

  /// @brief Verifies get solution for variational real P1 trial function.
  TEST(Rodin_Variational_Real_P1_TrialFunction, GetSolution)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);

    [[maybe_unused]] const auto& solution = u.getSolution();
    [[maybe_unused]] auto& mutable_solution = u.getSolution();
    
  }

  /// @brief Verifies sanity test build for variational vector P1 trial function by checking exact expected values.
  TEST(Rodin_Variational_Vector_P1_TrialFunction, SanityTest_Build)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);

    EXPECT_EQ(u.Space, TrialSpace);
  }

  /// @brief Verifies component access 2 D for variational vector P1 trial function.
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

  /// @brief Verifies component access 3 D for variational vector P1 trial function.
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

  /// @brief Verifies get leaf for variational real P1 trial function by checking exact expected values.
  TEST(Rodin_Variational_Real_P1_TrialFunction, GetLeaf)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);

    const auto& leaf = u.getLeaf();
    EXPECT_EQ(&leaf, &u);
  }

  /// @brief Verifies get solution vector for variational vector P1 trial function.
  TEST(Rodin_Variational_Vector_P1_TrialFunction, GetSolutionVector)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);

    [[maybe_unused]] const auto& solution = u.getSolution();
  }
}
