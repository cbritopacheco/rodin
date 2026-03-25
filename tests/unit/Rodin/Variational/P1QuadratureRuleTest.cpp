#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_P1QuadratureRule, LinearForm_ScalarLinearCoefficient_UsesMultiPointQuadrature)
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
    TestFunction v(fes);
    LinearForm lf(v);

    RealFunction c([](const Geometry::Point& p) { return p.x() + p.y(); });
    lf = Integral(c, v);
    lf.assemble();

    const auto& b = lf.getVector();
    ASSERT_EQ(b.size(), 3);

    // Exact values on the reference triangle:
    // ∫(x+y)φ0 = 1/12, ∫(x+y)φ1 = 1/8, ∫(x+y)φ2 = 1/8.
    EXPECT_NEAR(b(0), 1.0 / 12.0, 1e-12);
    EXPECT_NEAR(b(1), 1.0 / 8.0, 1e-12);
    EXPECT_NEAR(b(2), 1.0 / 8.0, 1e-12);
    EXPECT_NEAR(b.sum(), 1.0 / 3.0, 1e-12);
  }

  TEST(Rodin_Variational_P1QuadratureRule, BilinearForm_LinearCoefficientMass_UsesMultiPointQuadrature)
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
    BilinearForm bf(u, v);

    RealFunction c([](const Geometry::Point& p) { return p.x() + p.y(); });
    bf = Integral(Dot(Mult(c, u), v));
    bf.assemble();

    const auto& A = bf.getOperator();
    ASSERT_EQ(A.rows(), 3);
    ASSERT_EQ(A.cols(), 3);

    // Exact entries on the reference triangle.
    EXPECT_NEAR(A.coeff(0, 0), 1.0 / 30.0, 1e-12);
    EXPECT_NEAR(A.coeff(0, 1), 1.0 / 40.0, 1e-12);
    EXPECT_NEAR(A.coeff(0, 2), 1.0 / 40.0, 1e-12);
    EXPECT_NEAR(A.coeff(1, 1), 1.0 / 15.0, 1e-12);
    EXPECT_NEAR(A.coeff(2, 2), 1.0 / 15.0, 1e-12);
  }

  TEST(Rodin_Variational_P1QuadratureRule, MixedSpaces_GradGrad_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fesTr(mesh);
    P1 fesTe(mesh); // distinct space instance to trigger mixed-space path

    TrialFunction u(fesTr);
    TestFunction v(fesTe);

    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));
    bf.assemble();

    const auto& mat = bf.getOperator();
    EXPECT_EQ(mat.rows(), static_cast<Eigen::Index>(fesTe.getSize()));
    EXPECT_EQ(mat.cols(), static_cast<Eigen::Index>(fesTr.getSize()));
    EXPECT_GT(mat.norm(), 0.0);
  }

  TEST(Rodin_Variational_P1QuadratureRule, MixedSpaces_VectorMass_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const size_t vdim = 2;
    P1<Math::Vector<Real>> fesTr(mesh, vdim);
    P1<Math::Vector<Real>> fesTe(mesh, vdim);

    TrialFunction u(fesTr);
    TestFunction v(fesTe);

    RealFunction coeff(2.0);

    BilinearForm bf(u, v);
    bf = Integral(Dot(Mult(coeff, u), v));
    bf.assemble();

    const auto& mat = bf.getOperator();
    EXPECT_EQ(mat.rows(), static_cast<Eigen::Index>(fesTe.getSize()));
    EXPECT_EQ(mat.cols(), static_cast<Eigen::Index>(fesTr.getSize()));
    EXPECT_GT(mat.norm(), 0.0);
  }

  TEST(Rodin_Variational_P1QuadratureRule, MixedSpaces_VectorMass_NoCoeff_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const size_t vdim = 2;
    P1<Math::Vector<Real>> fesTr(mesh, vdim);
    P1<Math::Vector<Real>> fesTe(mesh, vdim);

    TrialFunction u(fesTr);
    TestFunction v(fesTe);

    BilinearForm bf(u, v);
    bf = Integral(Dot(u, v));
    bf.assemble();

    const auto& mat = bf.getOperator();
    EXPECT_EQ(mat.rows(), static_cast<Eigen::Index>(fesTe.getSize()));
    EXPECT_EQ(mat.cols(), static_cast<Eigen::Index>(fesTr.getSize()));
    EXPECT_GT(mat.norm(), 0.0);
  }

  TEST(Rodin_Variational_P1QuadratureRule, DivergencePressureCoupling_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const size_t vdim = 2;
    P1<Math::Vector<Real>> vel(mesh, vdim);
    P1<Real> pres(mesh);

    TrialFunction u(vel);
    TestFunction q(pres);

    BilinearForm b(u, q);
    b = Integral(Div(u), q);
    b.assemble();

    const auto& mat = b.getOperator();
    EXPECT_EQ(mat.rows(), static_cast<Eigen::Index>(pres.getSize()));
    EXPECT_EQ(mat.cols(), static_cast<Eigen::Index>(vel.getSize()));
    EXPECT_NE(mat.norm(), 0.0);
  }

  TEST(Rodin_Variational_P1QuadratureRule, PressureDivergence_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const size_t vdim = 2;
    P1<Real> pres(mesh);
    P1<Math::Vector<Real>> vel(mesh, vdim);

    TrialFunction p(pres);
    TestFunction v(vel);

    BilinearForm b(p, v);
    b = Integral(p, Div(v));
    b.assemble();

    const auto& mat = b.getOperator();
    EXPECT_EQ(mat.rows(), static_cast<Eigen::Index>(vel.getSize()));
    EXPECT_EQ(mat.cols(), static_cast<Eigen::Index>(pres.getSize()));
    EXPECT_NE(mat.norm(), 0.0);
  }
}
