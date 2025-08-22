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
  TEST(Rodin_Variational_Integral, TestFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);

    auto integral_v = Integral(v);

    EXPECT_EQ(integral_v.getRegion(), Region::Cells);
  }

  TEST(Rodin_Variational_Integral, TestFunction_Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);

    auto integral_v = Integral(v);
    auto copied = integral_v.copy();

    EXPECT_NE(copied, nullptr);
    EXPECT_EQ(copied->getRegion(), integral_v.getRegion());

    delete copied;
  }

  TEST(Rodin_Variational_Integral, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    auto integral_gf = Integral(gf);

    EXPECT_EQ(integral_gf.getRegion(), Region::Cells);
  }

  TEST(Rodin_Variational_Integral, GridFunction_Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    auto integral_gf = Integral(gf);
    auto copied = integral_gf.copy();

    EXPECT_NE(copied, nullptr);
    EXPECT_EQ(copied->getRegion(), integral_gf.getRegion());

    delete copied;
  }

  TEST(Rodin_Variational_Integral, Dot_ShapeFunctions_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    auto integral_dot = Integral(Grad(u), Grad(v));

    EXPECT_EQ(integral_dot.getRegion(), Region::Cells);
  }

  TEST(Rodin_Variational_Integral, Dot_ShapeFunctions_Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    auto integral_dot = Integral(Grad(u), Grad(v));
    auto copied = integral_dot.copy();

    EXPECT_NE(copied, nullptr);
    EXPECT_EQ(copied->getRegion(), integral_dot.getRegion());

    delete copied;
  }

  TEST(Rodin_Variational_Integral, InLinearForm_Assembly)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);

    lf = Integral(v);
    lf.assemble();

    const auto& vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);

    // For unit integrand, the integral should be positive (area of domain)
    Real total_integral = vector.sum();
    EXPECT_GT(total_integral, 0.0);
  }

  TEST(Rodin_Variational_Integral, InBilinearForm_MassMatrix)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);

    // Mass matrix: ∫ u v dx
    bf = Integral(u, v);
    bf.assemble();

    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
    EXPECT_EQ(op.rows(), op.cols());

    // Mass matrix should be positive definite (diagonal entries positive)
    for (Index i = 0; i < op.rows(); i++)
    {
      EXPECT_GT(op.coeff(i, i), 0.0);
    }
  }

  TEST(Rodin_Variational_Integral, InBilinearForm_StiffnessMatrix)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);

    // Stiffness matrix: ∫ ∇u · ∇v dx
    bf = Integral(Grad(u), Grad(v));
    bf.assemble();

    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
    EXPECT_EQ(op.rows(), op.cols());
  }

  TEST(Rodin_Variational_Integral, GridFunction_ConstantFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    gf = RealFunction(5.0);

    auto integral_gf = Integral(gf);

    // This integral represents ∫ 5 dx over the domain
    // The integration happens during assembly in linear/bilinear forms
    EXPECT_EQ(integral_gf.getRegion(), Region::Cells);
  }

  TEST(Rodin_Variational_Integral, VectorFunction_Components)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);

    // Test integration of vector function components
    auto u_x = u.x();
    auto v_x = v.x();
    auto integral_comp = Integral(u_x, v_x);

    EXPECT_EQ(integral_comp.getRegion(), Region::Cells);
  }

  TEST(Rodin_Variational_Integral, MoveConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);

    auto integral_v = Integral(v);
    auto region = integral_v.getRegion();

    auto integral_moved = std::move(integral_v);
    EXPECT_EQ(integral_moved.getRegion(), region);
  }

  TEST(Rodin_Variational_Integral, CopyConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);

    auto integral_v = Integral(v);
    auto integral_copy(integral_v);

    EXPECT_EQ(integral_copy.getRegion(), integral_v.getRegion());
  }

  TEST(Rodin_Variational_Integral, WithRealFunction_InLinearForm)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);

    RealFunction source(10.0);
    lf = Integral(source, v);
    lf.assemble();

    const auto& vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);

    // All entries should be positive since source is positive
    for (Index i = 0; i < static_cast<Index>(vector.size()); i++)
    {
      EXPECT_GT(vector(i), 0.0);
    }
  }

  TEST(Rodin_Variational_Integral, MultipleIntegrators)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);

    // Combine mass and stiffness: ∫ (u v + ∇u · ∇v) dx
    bf = Integral(u, v) + Integral(Grad(u), Grad(v));
    bf.assemble();

    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);

    // All diagonal entries should be positive
    for (Index i = 0; i < op.rows(); i++)
    {
      EXPECT_GT(op.coeff(i, i), 0.0);
    }
  }

  TEST(Rodin_Variational_Integral, LinearForm_WithVectorFunction)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TestFunction v(fes);
    LinearForm lf(v);

    VectorFunction force{1.0, 2.0};
    lf = Integral(Dot(v, force));
    lf.assemble();

    const auto& vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);
  }
}
