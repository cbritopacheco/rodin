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
  /// @brief Verifies sanity test build for variational real P1 bilinear form by checking exact expected values.
  TEST(Rodin_Variational_Real_P1_BilinearForm, SanityTest_Build)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    EXPECT_EQ(&bf.getTrialFunction(), &u);
    EXPECT_EQ(&bf.getTestFunction(), &v);
  }

  /// @brief Verifies copy constructor for variational real P1 bilinear form by checking exact expected values, copy semantics.
  TEST(Rodin_Variational_Real_P1_BilinearForm, CopyConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    BilinearForm bf_copy(bf);
    EXPECT_EQ(&bf_copy.getTrialFunction().getUUID(), &bf.getTrialFunction().getUUID());
    EXPECT_EQ(&bf_copy.getTestFunction().getUUID(), &bf.getTestFunction().getUUID());
  }

  /// @brief Verifies move constructor for variational real P1 bilinear form by checking exact expected values, move semantics.
  TEST(Rodin_Variational_Real_P1_BilinearForm, MoveConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    BilinearForm bf_moved(std::move(bf));
    EXPECT_EQ(&bf_moved.getTrialFunction(), &u);
    EXPECT_EQ(&bf_moved.getTestFunction(), &v);
  }

  /// @brief Verifies assignment for variational real P1 bilinear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_BilinearForm, Assignment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));
    EXPECT_FALSE(bf.getLocalIntegrators().empty());
  }

  /// @brief Verifies addition assignment for variational real P1 bilinear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_BilinearForm, AdditionAssignment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    bf += Integral(Grad(u), Grad(v));
    EXPECT_FALSE(bf.getLocalIntegrators().empty());
  }

  /// @brief Verifies subtraction assignment for variational real P1 bilinear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_BilinearForm, SubtractionAssignment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    bf -= Integral(Grad(u), Grad(v));
    EXPECT_FALSE(bf.getLocalIntegrators().empty());
  }

  /// @brief Verifies assemble and get operator for variational real P1 bilinear form by checking exact expected values, form assembly.
  TEST(Rodin_Variational_Real_P1_BilinearForm, AssembleAndGetOperator)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));
    bf.assemble();
    const auto& op = bf.getOperator();
    auto& mutable_op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
    EXPECT_EQ(op.rows(), mutable_op.rows());
    EXPECT_EQ(op.cols(), mutable_op.cols());
  }

  /// @brief Verifies copy for variational real P1 bilinear form by checking copy semantics.
  TEST(Rodin_Variational_Real_P1_BilinearForm, Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));
    auto copied = bf.copy();
    EXPECT_NE(copied, nullptr);
    delete copied;
  }

  /// @brief Verifies sanity test build for variational vector P1 bilinear form by checking exact expected values.
  TEST(Rodin_Variational_Vector_P1_BilinearForm, SanityTest_Build)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    EXPECT_EQ(&bf.getTrialFunction(), &u);
    EXPECT_EQ(&bf.getTestFunction(), &v);
  }

  /// @brief Verifies elasticity integrator for variational vector P1 bilinear form by checking form assembly.
  TEST(Rodin_Variational_Vector_P1_BilinearForm, ElasticityIntegrator)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    // Elasticity-like integrator
    bf = Integral(Dot(Jacobian(u), Jacobian(v)));
    bf.assemble();
    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
  }

  /// @brief Verifies mass matrix for variational real P1 bilinear form by checking form assembly.
  TEST(Rodin_Variational_Real_P1_BilinearForm, MassMatrix)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);
    // Mass matrix integrator
    bf = Integral(u, v);
    bf.assemble();
    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
  }

  /// @brief Verifies sparse bilinear-form action preserves trial/test argument order.
  TEST(Rodin_Variational_Real_P1_BilinearForm, SparseActionUsesTrialThenTest)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {{0, 1, 2}})
      .finalize();

    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm form(u, v);

    auto& A = form.getOperator();
    A.resize(3, 3);
    A.insert(0, 0) = 1.0;
    A.insert(0, 1) = 2.0;
    A.insert(1, 2) = -3.0;
    A.insert(2, 0) = 4.0;
    A.makeCompressed();

    GridFunction trial(fes);
    trial.getData() << 1.0, 2.0, -1.0;
    GridFunction test(fes);
    test.getData() << -2.0, 1.0, 3.0;

    const Real expected = (A * trial.getData()).dot(test.getData());
    const Real reversed = (A * test.getData()).dot(trial.getData());
    EXPECT_NE(expected, reversed);
    EXPECT_DOUBLE_EQ(form(trial, test), expected);
  }

  /// @brief Verifies dense bilinear-form action preserves trial/test argument order.
  TEST(Rodin_Variational_Real_P1_BilinearForm, DenseActionUsesTrialThenTest)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {{0, 1, 2}})
      .finalize();

    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    using Solution = typename decltype(u)::SolutionType;
    BilinearForm<Solution, decltype(fes), decltype(fes), Math::Matrix<Real>> form(u, v);

    auto& A = form.getOperator();
    A.resize(3, 3);
    A << 1.0, 2.0, 0.0,
         0.0, 0.0, -3.0,
         4.0, 0.0, 0.0;

    GridFunction trial(fes);
    trial.getData() << 1.0, 2.0, -1.0;
    GridFunction test(fes);
    test.getData() << -2.0, 1.0, 3.0;

    const Real expected = (A * trial.getData()).dot(test.getData());
    const Real reversed = (A * test.getData()).dot(trial.getData());
    EXPECT_NE(expected, reversed);
    EXPECT_DOUBLE_EQ(form(trial, test), expected);
  }
}
