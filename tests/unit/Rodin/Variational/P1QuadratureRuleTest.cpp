/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <type_traits>

#include "Rodin/Variational.h"
#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief A constant nonlocal kernel has an independently known rank-one matrix.
  TEST(Rodin_Variational_P1QuadratureRule, NonlocalDenseAndSparseFormsMatchExactMoments)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    P1<Real> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    BilinearForm sparse(u, v);
    BilinearForm<typename decltype(u)::SolutionType, decltype(space), decltype(space),
      Math::Matrix<Real>>
      dense(u, v);
    auto kernel = [](const Point&, const Point&) { return 2.0; };
    auto integral = Integral(Potential(kernel, u), v);
    Math::Vector<Real> moments(3);
    moments << 0.5, 1.0, 0.5;
    const Math::Matrix<Real> expected = 2 * moments * moments.transpose();
    sparse += integral;
    dense += integral;
    for (int repeat = 0; repeat < 2; ++repeat)
    {
      sparse.assemble();
      dense.assemble();
      EXPECT_LT((Math::Matrix<Real>(sparse.getOperator()) - expected).norm(), 1e-12);
      EXPECT_LT((dense.getOperator() - expected).norm(), 1e-12) << dense.getOperator();
    }
    auto excluded = integral;
    excluded.setTrialAttributes({999});
    sparse.getGlobalIntegrators().clear();
    sparse += excluded;
    dense.getGlobalIntegrators().clear();
    dense += excluded;
    sparse.assemble();
    dense.assemble();
    EXPECT_EQ(sparse.getOperator().norm(), 0);
    EXPECT_EQ(dense.getOperator().norm(), 0);
    auto excludedTest = integral;
    excludedTest.setTestAttributes({999});
    sparse.getGlobalIntegrators().clear();
    sparse += excludedTest;
    dense.getGlobalIntegrators().clear();
    dense += excludedTest;
    sparse.assemble();
    dense.assemble();
    EXPECT_EQ(sparse.getOperator().norm(), 0);
    EXPECT_EQ(dense.getOperator().norm(), 0);
    Problem problem(u, v);
    problem = integral;
    problem.assemble();
    EXPECT_LT(
      (Math::Matrix<Real>(problem.getLinearSystem().getOperator()) - expected).norm(),
      1e-12);
  }

  /**
   * @brief Checks empty component ranges and scalar-cache publication on rebinding.
   *
   * For @f$V_h=P_1\otimes\mathbb K^0=\{0\}@f$, the mass, weighted mass
   * and directional-derivative forms have no local entries. Binding must not
   * evaluate coefficients or divide by the component count. Point cells are
   * embedded in three dimensions; all positive-dimensional cell families are
   * also included, over real and complex fields.
   *
   * ## State transitions
   *
   * A local integrator is rebound after changes between zero, two and three
   * components and quadrature orders zero and four. No global linear system
   * is retained across space changes. Every nonempty entry must equal a fresh
   * integrator exactly. Logical entity indices and coefficient-call counts
   * check the binding contract without coordinate comparisons or tolerances.
   */
  TEST(Rodin_Variational_P1QuadratureRule, EmptyComponentsAndRebinding)
  {
    const auto checkScalar = []<class Scalar>() {
      using Range = Math::SpatialVector<Scalar>;
      for (const auto geometry :
        {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
          Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
          Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
      {
        SCOPED_TRACE(int(geometry));
        SCOPED_TRACE((std::is_same_v<Scalar, Complex> ? "Complex" : "Real"));
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        auto mesh = [&]() {
          if (dimension == 0)
          {
            auto builder = LocalMesh::Builder();
            builder.initialize(3).nodes(1).vertex({0, 0, 0});
            return builder.finalize();
          }
          const Array<size_t> shape = Array<size_t>::Constant(dimension, 2);
          return LocalMesh::UniformGrid(geometry, shape);
        }();
        for (size_t d = 1; d <= dimension; ++d)
        {
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
        P1<Range> space(mesh, 0);
        TrialFunction u(space);
        TestFunction v(space);
        size_t coefficientCalls = 0;
        const auto coefficientCallable = [&](const Point&) -> Scalar {
          ++coefficientCalls;
          return Scalar(2);
        };
        const auto coefficient = [&]() {
          if constexpr (std::is_same_v<Scalar, Real>)
            return RealFunction(coefficientCallable);
          else
            return ComplexFunction(coefficientCallable);
        }();
        const size_t ambient = mesh.getSpaceDimension();
        const VectorFunction direction(ambient, [&](const Point&) {
          ++coefficientCalls;
          Math::SpatialVector<Real> value(ambient);
          for (size_t axis = 0; axis < ambient; ++axis)
            value[axis] = Real(axis + 1);
          return value;
        });
        const auto check = [&](auto makeForm, bool hasCoefficient) {
          auto reused = makeForm();
          for (const size_t components :
            {size_t(0), size_t(2), size_t(0), size_t(3), size_t(0)})
          {
            SCOPED_TRACE(components);
            space = P1<Range>(mesh, components);
            ASSERT_EQ(space.getSize(), mesh.getVertexCount() * components);
            for (const size_t order : {size_t(0), size_t(4), size_t(0)})
            {
              SCOPED_TRACE(order);
              reused.setOrder(order);
              for (auto cell = mesh.getCell(); cell; ++cell)
              {
                SCOPED_TRACE(cell->getIndex());
                const size_t expectedCalls = components == 0 || !hasCoefficient
                  ? 0
                  : QF::QuadratureFormula::get(order, geometry).getSize();
                coefficientCalls = 0;
                ASSERT_NO_THROW(reused.setPolytope(*cell));
                EXPECT_EQ(coefficientCalls, expectedCalls);
                EXPECT_EQ(reused.getPolytope().getDimension(), cell->getDimension());
                EXPECT_EQ(reused.getPolytope().getIndex(), cell->getIndex());
                auto fresh = makeForm();
                fresh.setOrder(order);
                coefficientCalls = 0;
                ASSERT_NO_THROW(fresh.setPolytope(*cell));
                EXPECT_EQ(coefficientCalls, expectedCalls);
                for (size_t te = 0; te < v.getDOFs(*cell); ++te)
                {
                  for (size_t tr = 0; tr < u.getDOFs(*cell); ++tr)
                  {
                    const auto actual = reused.integrate(tr, te);
                    EXPECT_EQ(actual, fresh.integrate(tr, te));
                  }
                }
              }
            }
          }
        };
        check([&]() { return Integral(u, v); }, false);
        check([&]() { return Integral(coefficient * u, v); }, true);
        check([&]() { return Integral(Jacobian(u) * direction, v); }, true);
      }
    };
    checkScalar.template operator()<Real>();
    checkScalar.template operator()<Complex>();
  }

  /// @brief Verifies linear form scalar linear coefficient uses multi point quadrature for variational P1 quadrature rule by checking tolerance-based numerical results, exact expected values, form assembly.
  TEST(Rodin_Variational_P1QuadratureRule,
    LinearForm_ScalarLinearCoefficient_UsesMultiPointQuadrature)
  {
    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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
    // The coefficient is linear, but built from a lambda, so it reports no
    // degree and order inference sees only the basis. The integrand's own
    // degree is stated instead, which is what setOrder is for.
    auto integ = Integral(c, v);
    integ.setOrder(2);
    lf = integ;
    lf.assemble();

    const auto& b = lf.getVector();
    ASSERT_EQ(b.size(), 3);

    // Exact values on the reference triangle:
    // \int(x+y)φ0 = 1/12, \int(x+y)φ1 = 1/8, \int(x+y)φ2 = 1/8.
    EXPECT_NEAR(b(0), 1.0 / 12.0, 1e-12);
    EXPECT_NEAR(b(1), 1.0 / 8.0, 1e-12);
    EXPECT_NEAR(b(2), 1.0 / 8.0, 1e-12);
    EXPECT_NEAR(b.sum(), 1.0 / 3.0, 1e-12);
  }

  /// @brief Verifies bilinear form linear coefficient mass uses multi point quadrature for variational P1 quadrature rule by checking tolerance-based numerical results, exact expected values, form assembly.
  TEST(Rodin_Variational_P1QuadratureRule,
    BilinearForm_LinearCoefficientMass_UsesMultiPointQuadrature)
  {
    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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
    // The coefficient is linear, but built from a lambda, so it reports no
    // degree and order inference sees only the basis. The integrand's own
    // degree is stated instead, which is what setOrder is for.
    auto integ = Integral(Dot(Mult(c, u), v));
    integ.setOrder(3);
    bf = integ;
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

  /// @brief Verifies mixed spaces grad grad assembles for variational P1 quadrature rule by checking exact expected values, form assembly.
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

  /// @brief Verifies mixed spaces vector mass assembles for variational P1 quadrature rule by checking exact expected values, form assembly.
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

  /// @brief Verifies mixed spaces vector mass no coeff assembles for variational P1 quadrature rule by checking exact expected values, form assembly.
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

  /// @brief Verifies vector mass matches grid function linear form for variational P1 quadrature rule by checking tolerance-based numerical results, exact expected values, form assembly.
  TEST(Rodin_Variational_P1QuadratureRule, VectorMass_MatchesGridFunctionLinearForm)
  {
    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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

  /// @brief Checks weighted vector mass against the grid-function linear form.
  TEST(Rodin_Variational_P1QuadratureRule,
    ScalarWeightedVectorMass_MatchesGridFunctionLinearForm)
  {
    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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

    RealFunction coeff(
      [](const Geometry::Point& p) { return 1.0 + p.x() + 2.0 * p.y(); });

    BilinearForm mass(u, v);
    // The coefficient is linear, but built from a lambda, so it reports no
    // degree and order inference sees only the basis. The integrand's own
    // degree is stated instead, which is what setOrder is for.
    auto massInteg = Integral(Dot(Mult(coeff, u), v));
    massInteg.setOrder(3);
    mass = massInteg;
    mass.assemble();

    LinearForm load(v);
    // The same order as the bilinear form above, so both are integrated
    // with one rule and the identity between them is exact.
    auto loadInteg = Integral(Dot(coeff * gf, v));
    loadInteg.setOrder(3);
    load = loadInteg;
    load.assemble();

    const auto residual = mass.getOperator() * x - load.getVector();
    EXPECT_NEAR(residual.norm(), 0.0, 1e-14);
  }

  /// @brief Verifies outer scalar weighted vector mass matches grid function linear form for variational P1 quadrature rule by checking tolerance-based numerical results, exact expected values, form assembly.
  TEST(Rodin_Variational_P1QuadratureRule,
    OuterScalarWeightedVectorMass_MatchesGridFunctionLinearForm)
  {
    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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

    RealFunction coeff(
      [](const Geometry::Point& p) { return 1.0 + p.x() + 2.0 * p.y(); });

    BilinearForm mass(u, v);
    // The coefficient is linear, but built from a lambda, so it reports no
    // degree and order inference sees only the basis. The integrand's own
    // degree is stated instead, which is what setOrder is for.
    auto massInteg = Integral(coeff * Dot(u, v));
    massInteg.setOrder(3);
    mass = massInteg;
    mass.assemble();

    LinearForm load(v);
    // The same order as the bilinear form above, so both are integrated
    // with one rule and the identity between them is exact.
    auto loadInteg = Integral(coeff * Dot(gf, v));
    loadInteg.setOrder(3);
    load = loadInteg;
    load.assemble();

    const auto residual = mass.getOperator() * x - load.getVector();
    EXPECT_NEAR(residual.norm(), 0.0, 1e-14);
  }

  /// @brief Verifies matrix weighted vector mass matches grid function linear form for variational P1 quadrature rule by checking tolerance-based numerical results, exact expected values, form assembly.
  TEST(Rodin_Variational_P1QuadratureRule,
    MatrixWeightedVectorMass_MatchesGridFunctionLinearForm)
  {
    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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

  /// @brief Verifies vector potential constant kernel uses separate trial and test basis values for variational P1 quadrature rule by checking tolerance-based numerical results, exact expected values, form assembly.
  TEST(Rodin_Variational_P1QuadratureRule,
    VectorPotentialConstantKernel_UsesSeparateTrialAndTestBasisValues)
  {
    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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

    auto kernel = [](Math::SpatialMatrix<Real>& out, const Point&, const Point&) {
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
            const Eigen::Index row = static_cast<Eigen::Index>(testVertex + testComp * 3);
            const Eigen::Index col =
              static_cast<Eigen::Index>(trialVertex + trialComp * 3);
            EXPECT_NEAR(A(row, col),
              scalarA(static_cast<Eigen::Index>(testVertex),
                static_cast<Eigen::Index>(trialVertex)) *
                k[testComp][trialComp],
              1e-12);
          }
        }
      }
    }
  }

  /// @brief Verifies a complex potential kernel is not conjugated during contraction or assembly.
  TEST(Rodin_Variational_P1QuadratureRule, ComplexVectorPotentialPreservesKernel)
  {
    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
                  .initialize(2)
                  .nodes(3)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({0, 1})
                  .polytope(Polytope::Type::Triangle, {{0, 1, 2}})
                  .finalize();

    P1<Real> scalarFES(mesh);
    TrialFunction scalarU(scalarFES);
    TestFunction scalarV(scalarFES);
    DenseProblem scalarProblem(scalarU, scalarV);
    scalarProblem = Integral(
      Potential([](const Point&, const Point&) { return 1.0; }, scalarU), scalarV);
    scalarProblem.assemble();
    const auto& scalarA = scalarProblem.getLinearSystem().getOperator();

    P1<Math::SpatialVector<Complex>> fes(mesh, 2);
    TrialFunction u(fes);
    TestFunction v(fes);
    const Complex entries[2][2] = {
      {Complex(1, 2), Complex(-3, 1)}, {Complex(2, -4), Complex(5, 3)}};
    auto kernel = [&](Math::SpatialMatrix<Complex>& out, const Point&, const Point&) {
      out.resize(2, 2);
      for (size_t i = 0; i < 2; ++i)
      {
        for (size_t j = 0; j < 2; ++j)
          out(i, j) = entries[i][j];
      }
    };

    DenseProblem problem(u, v);
    problem = Integral(Potential(kernel, u), v);
    problem.assemble();
    const auto& A = problem.getLinearSystem().getOperator();

    for (size_t testVertex = 0; testVertex < 3; ++testVertex)
    {
      for (size_t trialVertex = 0; trialVertex < 3; ++trialVertex)
      {
        for (size_t testComp = 0; testComp < 2; ++testComp)
        {
          for (size_t trialComp = 0; trialComp < 2; ++trialComp)
          {
            const Eigen::Index row = static_cast<Eigen::Index>(testVertex + testComp * 3);
            const Eigen::Index col =
              static_cast<Eigen::Index>(trialVertex + trialComp * 3);
            const Complex expected = scalarA(static_cast<Eigen::Index>(testVertex),
                                       static_cast<Eigen::Index>(trialVertex)) *
              entries[testComp][trialComp];
            EXPECT_NEAR(std::abs(A(row, col) - expected), 0.0, 1e-12);
          }
        }
      }
    }
  }

  /// @brief Verifies divergence pressure coupling assembles for variational P1 quadrature rule by checking exact expected values, form assembly.
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

  /// @brief Verifies pressure divergence assembles for variational P1 quadrature rule by checking exact expected values, form assembly.
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

  /// @brief Verifies divergence couplings transform gradients on a sheared tetrahedron.
  TEST(Rodin_Variational_P1QuadratureRule,
    DivergenceCouplingsTransformGradientsOnShearedTetrahedron)
  {
    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
                  .initialize(3)
                  .nodes(4)
                  .vertex({0, 0, 0})
                  .vertex({2, 0, 0})
                  .vertex({1, 3, 0})
                  .vertex({0, 1, 4})
                  .polytope(Polytope::Type::Tetrahedron, {{0, 1, 2, 3}})
                  .finalize();

    P1<Math::SpatialVector<Real>> velocityFES(mesh, 3);
    P1<Real> pressureFES(mesh);
    GridFunction velocity(velocityFES);
    velocity.getData() << 0.0, 2.0, 7.0, 2.0, 0.0, -2.0, 8.0, 19.0, 0.0, 4.0, 2.0, -8.0;
    GridFunction pressure(pressureFES);
    pressure.getData().setOnes();

    TrialFunction u(velocityFES);
    TestFunction q(pressureFES);
    BilinearForm divergencePressure(u, q);
    divergencePressure = Integral(Div(u), q);
    divergencePressure.assemble();

    TrialFunction p(pressureFES);
    TestFunction v(velocityFES);
    BilinearForm pressureDivergence(p, v);
    pressureDivergence = Integral(p, Div(v));
    pressureDivergence.assemble();

    // The affine velocity has gradient
    // [ 1  2  0; -1  3  4; 2  0 -2 ], hence divergence 2.
    // The tetrahedron has volume 4, so both actions equal 8.
    const Real trialDivergence =
      (divergencePressure.getOperator() * velocity.getData()).dot(pressure.getData());
    const Real testDivergence =
      (pressureDivergence.getOperator() * pressure.getData()).dot(velocity.getData());
    EXPECT_NEAR(trialDivergence, 8.0, 1e-12);
    EXPECT_NEAR(testDivergence, 8.0, 1e-12);
  }
}
