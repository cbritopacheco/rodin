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
    P1<Math::SpatialVector<Real>> fesTr(mesh, vdim);
    P1<Math::SpatialVector<Real>> fesTe(mesh, vdim);

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
    P1<Math::SpatialVector<Real>> fesTr(mesh, vdim);
    P1<Math::SpatialVector<Real>> fesTe(mesh, vdim);

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

  TEST(Rodin_Variational_P1QuadratureRule, VectorMass_MatchesGridFunctionLinearForm)
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

    const size_t vdim = 2;
    P1<Math::SpatialVector<Real>> fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);

    GridFunction gf(fes);
    auto& x = gf.getData();
    ASSERT_EQ(x.size(), static_cast<Eigen::Index>(fes.getSize()));
    for (Eigen::Index i = 0; i < x.size(); ++i)
      x(i) = static_cast<Real>(i + 1);

    BilinearForm mass(u, v);
    mass = Integral(Dot(u, v));
    mass.assemble();

    LinearForm load(v);
    load = Integral(Dot(gf, v));
    load.assemble();

    const auto residual = mass.getOperator() * x - load.getVector();
    EXPECT_NEAR(residual.norm(), 0.0, 1e-14);
  }

  TEST(Rodin_Variational_P1QuadratureRule, ScalarWeightedVectorMass_MatchesGridFunctionLinearForm)
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

    const size_t vdim = 2;
    P1<Math::SpatialVector<Real>> fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);

    GridFunction gf(fes);
    auto& x = gf.getData();
    ASSERT_EQ(x.size(), static_cast<Eigen::Index>(fes.getSize()));
    for (Eigen::Index i = 0; i < x.size(); ++i)
      x(i) = static_cast<Real>(i + 1);

    RealFunction coeff([](const Geometry::Point& p) { return 1.0 + p.x() + 2.0 * p.y(); });

    BilinearForm mass(u, v);
    mass = Integral(Dot(Mult(coeff, u), v));
    mass.assemble();

    LinearForm load(v);
    load = Integral(Dot(coeff * gf, v));
    load.assemble();

    const auto residual = mass.getOperator() * x - load.getVector();
    EXPECT_NEAR(residual.norm(), 0.0, 1e-14);
  }

  TEST(Rodin_Variational_P1QuadratureRule, OuterScalarWeightedVectorMass_MatchesGridFunctionLinearForm)
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

    const size_t vdim = 2;
    P1<Math::SpatialVector<Real>> fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);

    GridFunction gf(fes);
    auto& x = gf.getData();
    ASSERT_EQ(x.size(), static_cast<Eigen::Index>(fes.getSize()));
    for (Eigen::Index i = 0; i < x.size(); ++i)
      x(i) = static_cast<Real>(i + 1);

    RealFunction coeff([](const Geometry::Point& p) { return 1.0 + p.x() + 2.0 * p.y(); });

    BilinearForm mass(u, v);
    mass = Integral(coeff * Dot(u, v));
    mass.assemble();

    LinearForm load(v);
    load = Integral(coeff * Dot(gf, v));
    load.assemble();

    const auto residual = mass.getOperator() * x - load.getVector();
    EXPECT_NEAR(residual.norm(), 0.0, 1e-14);
  }

  TEST(Rodin_Variational_P1QuadratureRule, MatrixWeightedVectorMass_MatchesGridFunctionLinearForm)
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

    const size_t vdim = 2;
    P1<Math::SpatialVector<Real>> fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);

    GridFunction gf(fes);
    auto& x = gf.getData();
    ASSERT_EQ(x.size(), static_cast<Eigen::Index>(fes.getSize()));
    for (Eigen::Index i = 0; i < x.size(); ++i)
      x(i) = static_cast<Real>(i + 1);

    Math::Matrix<Real> matrix(2, 2);
    matrix(0, 0) = 2.0;
    matrix(0, 1) = -0.25;
    matrix(1, 0) = 0.75;
    matrix(1, 1) = 1.5;
    MatrixFunction coeff(matrix);

    BilinearForm mass(u, v);
    mass = Integral(Dot(coeff * u, v));
    mass.assemble();

    LinearForm load(v);
    load = Integral(Dot(coeff * gf, v));
    load.assemble();

    const auto residual = mass.getOperator() * x - load.getVector();
    EXPECT_NEAR(residual.norm(), 0.0, 1e-14);
  }

  TEST(Rodin_Variational_P1QuadratureRule, VectorPotentialConstantKernel_UsesSeparateTrialAndTestBasisValues)
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

    const size_t vdim = 2;
    P1<Math::SpatialVector<Real>> fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);

    P1<Real> scalarFes(mesh);
    TrialFunction scalarU(scalarFes);
    TestFunction scalarV(scalarFes);

    auto scalarKernel = [](const Point&, const Point&) { return 1.0; };
    DenseProblem scalarProblem(scalarU, scalarV);
    scalarProblem = Integral(Potential(scalarKernel, scalarU), scalarV);
    scalarProblem.assemble();
    const auto& scalarA = scalarProblem.getLinearSystem().getOperator();

    auto kernel =
      [](Math::SpatialMatrix<Real>& out, const Point&, const Point&)
      {
        out.resize(2, 2);
        out(0, 0) = 2.0;
        out(0, 1) = -3.0;
        out(1, 0) = 5.0;
        out(1, 1) = 7.0;
      };

    DenseProblem problem(u, v);
    problem = Integral(Potential(kernel, u), v);
    problem.assemble();

    const auto& A = problem.getLinearSystem().getOperator();
    ASSERT_EQ(A.rows(), static_cast<Eigen::Index>(fes.getSize()));
    ASSERT_EQ(A.cols(), static_cast<Eigen::Index>(fes.getSize()));
    ASSERT_EQ(scalarA.rows(), 3);
    ASSERT_EQ(scalarA.cols(), 3);

    const Real k[2][2] = {{2.0, -3.0}, {5.0, 7.0}};
    for (size_t testVertex = 0; testVertex < 3; ++testVertex)
    {
      for (size_t trialVertex = 0; trialVertex < 3; ++trialVertex)
      {
        for (size_t testComp = 0; testComp < vdim; ++testComp)
        {
          for (size_t trialComp = 0; trialComp < vdim; ++trialComp)
          {
            const Eigen::Index row =
              static_cast<Eigen::Index>(testVertex + testComp * 3);
            const Eigen::Index col =
              static_cast<Eigen::Index>(trialVertex + trialComp * 3);
            EXPECT_NEAR(
              A(row, col),
              scalarA(
                static_cast<Eigen::Index>(testVertex),
                static_cast<Eigen::Index>(trialVertex)) * k[testComp][trialComp],
              1e-12);
          }
        }
      }
    }
  }

  TEST(Rodin_Variational_P1QuadratureRule, DivergencePressureCoupling_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const size_t vdim = 2;
    P1<Math::SpatialVector<Real>> vel(mesh, vdim);
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
    P1<Math::SpatialVector<Real>> vel(mesh, vdim);

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
