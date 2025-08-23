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
}
