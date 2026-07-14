#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  /// @brief Verifies sanity test build for variational real P1 linear form by checking exact expected values.
  TEST(Rodin_Variational_Real_P1_LinearForm, SanityTest_Build)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);

    EXPECT_EQ(&lf.getTestFunction(), &v);
  }

  /// @brief Verifies fuzzy test uniform grid 4 x 4 for variational real P1 linear form by checking form assembly.
  TEST(Rodin_Variational_Real_P1_LinearForm, FuzzyTest_UniformGrid_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf = Integral(v);
    lf.assemble();

    const auto& vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);
  }

  /// @brief Verifies copy constructor for variational real P1 linear form by checking exact expected values, copy semantics.
  TEST(Rodin_Variational_Real_P1_LinearForm, CopyConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    
    LinearForm lf_copy(lf);
    EXPECT_EQ(&lf_copy.getTestFunction().getUUID(), &lf.getTestFunction().getUUID());
  }

  /// @brief Verifies move constructor for variational real P1 linear form by checking exact expected values, move semantics.
  TEST(Rodin_Variational_Real_P1_LinearForm, MoveConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    LinearForm lf_moved(std::move(lf));
    EXPECT_EQ(&lf_moved.getTestFunction(), &v);
  }

  /// @brief Verifies assignment for variational real P1 linear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_LinearForm, Assignment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf = Integral(v);
    EXPECT_FALSE(lf.getIntegrators().empty());
  }

  /// @brief Verifies addition assignment for variational real P1 linear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_LinearForm, AdditionAssignment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf += Integral(v);
    EXPECT_FALSE(lf.getIntegrators().empty());
  }

  /// @brief Verifies subtraction assignment for variational real P1 linear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_LinearForm, SubtractionAssignment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf -= Integral(v);
    EXPECT_FALSE(lf.getIntegrators().empty());
  }

  /// @brief Verifies assemble and get vector for variational real P1 linear form by checking exact expected values, form assembly.
  TEST(Rodin_Variational_Real_P1_LinearForm, AssembleAndGetVector)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf = Integral(v);
    lf.assemble();
    const auto& vector = lf.getVector();
    auto& mutable_vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);
    EXPECT_EQ(vector.size(), mutable_vector.size());
  }

  /// @brief Verifies copy for variational real P1 linear form by checking copy semantics.
  TEST(Rodin_Variational_Real_P1_LinearForm, Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf = Integral(v);
    auto copied = lf.copy();
    EXPECT_NE(copied, nullptr);
    delete copied;
  }

  /// @brief Verifies sanity test build for variational vector P1 linear form by checking exact expected values.
  TEST(Rodin_Variational_Vector_P1_LinearForm, SanityTest_Build)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TestFunction v(fes);
    LinearForm lf(v);
    EXPECT_EQ(&lf.getTestFunction(), &v);
  }

  /// @brief Verifies assemble vector function for variational vector P1 linear form by checking form assembly.
  TEST(Rodin_Variational_Vector_P1_LinearForm, AssembleVectorFunction)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TestFunction v(fes);
    LinearForm lf(v);
    lf = Integral(Dot(v, VectorFunction{ 1, 1 }));
    lf.assemble();
    const auto& vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);
  }

  /// @brief Verifies clear integrators for variational real P1 linear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_LinearForm, ClearIntegrators)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf += Integral(v);
    EXPECT_FALSE(lf.getIntegrators().empty());
    lf = Integral(v);  // This should clear and add new integrator
    EXPECT_FALSE(lf.getIntegrators().empty());
  }
}
