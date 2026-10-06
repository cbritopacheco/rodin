/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Matrix Jacobian and nonlocal quadrature regressions. */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Assembly.h"
#include "Rodin/Geometry/ParametricTransformation.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  using G = Polytope::Type;
  constexpr size_t quadratureOrder = 6;
  // Entry comparisons use the same quadrature formula, so only rounding differs.
  constexpr Real tolerance = 1e-10;

  template <class T>
  concept ScalarIntegral = requires(const T& value) { Integral(value); };

  Mesh<> makeMesh(G geometry)
  {
    const size_t d = Polytope::Traits(geometry).getDimension();
    auto mesh = d == 1 ? Mesh<>::UniformGrid(geometry, {2})
      : d == 2         ? Mesh<>::UniformGrid(geometry, {2, 2})
                       : Mesh<>::UniformGrid(geometry, {2, 2, 2});
    for (size_t from = 1; from <= d; ++from)
      for (size_t to = 0; to < from; ++to)
        mesh.getConnectivity().compute(from, to);
    return mesh;
  }

  template <class Integrand>
  void checkEntries(Integrand integrand, const Polytope& cell)
  {
    auto rule = Integral(integrand);
    rule.setOrder(quadratureOrder);
    rule.setPolytope(cell);
    const size_t ntr = integrand.getLHS().getDOFs(cell);
    const size_t nte = integrand.getRHS().getDOFs(cell);
    const std::vector<size_t> trials = {0, ntr / 3, ntr / 2, ntr - 1};
    const std::vector<size_t> tests = {0, nte / 3, nte / 2, nte - 1};
    using Scalar = typename FormLanguage::Traits<Integrand>::ScalarType;
    Math::Matrix<Scalar> expected = Math::Matrix<Scalar>::Zero(4, 4);
    const auto& formula =
      QF::PolytopeQuadratureFormula::get(quadratureOrder, cell.getGeometry());
    const auto& quadrature = cell.getQuadrature(formula);
    // Evaluate the expression directly, independently of specialized tabulation
    // and local-matrix assembly. Rectangular shapes expose component permutations.
    for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
    {
      const auto& point = quadrature.getPoint(qp);
      const IntegrationPoint ip(point, &formula, qp);
      integrand.setIntegrationPoint(ip);
      const Real weight = formula.getWeight(qp) * point.getDistortion();
      for (size_t i = 0; i < tests.size(); ++i)
        for (size_t j = 0; j < trials.size(); ++j)
          expected(i, j) += weight *
            Math::dot(integrand.getLHS().getBasis(trials[j]),
              integrand.getRHS().getBasis(tests[i]));
    }
    auto copied = rule;
    auto moved = std::move(copied);
    moved.setPolytope(cell);
    for (size_t i = 0; i < tests.size(); ++i)
      for (size_t j = 0; j < trials.size(); ++j)
      {
        EXPECT_LE(std::abs(rule.integrate(trials[j], tests[i]) - expected(i, j)),
          tolerance * (1 + std::abs(expected(i, j))));
        EXPECT_LE(std::abs(moved.integrate(trials[j], tests[i]) - expected(i, j)),
          tolerance * (1 + std::abs(expected(i, j))));
      }
  }

  template <class TrialFES, class TestFES>
  void checkJacobians(const TrialFES& trialFES, const TestFES& testFES)
  {
    using Scalar = typename FormLanguage::Traits<TrialFES>::ScalarType;
    const size_t dimension = trialFES.getMesh().getSpaceDimension();
    Math::SpatialTensor<Scalar> first(
      trialFES.getRows(), trialFES.getColumns(), dimension);
    auto second = first;
    for (size_t i = 0; i < first.size(); ++i)
    {
      first[i] = Scalar(1 + Real(i) / 10);
      second[i] = Scalar(-2 + Real(i) / 7);
      if constexpr (std::is_same_v<Scalar, Complex>)
      {
        first[i] += Complex(0, Real(i + 1) / 11);
        second[i] += Complex(0, -Real(i + 2) / 13);
      }
    }
    TensorFunction A(first), B(second);
    VectorFunction velocity(dimension, [dimension](const Point& p) {
      Math::SpatialVector<Real> value(dimension);
      for (size_t k = 0; k < dimension; ++k)
        value(k) = 1 + Real(k) + p.getPhysicalCoordinates()(k);
      return value;
    });
    RealFunction weight([](const Point& p) { return 2 + p.x(); });
    TrialFunction u(trialFES);
    TestFunction v(testFES);
    // Matrix Jacobians are rank three: vector-only handlers must not match.
    static_assert(
      !FormLanguage::IsSpecialized<
        QuadratureRule<decltype(Dot(weight * Jacobian(u), Jacobian(v)))>>::Value);
    static_assert(
      !FormLanguage::IsSpecialized<
        QuadratureRule<decltype(Dot(Dot(A, Jacobian(u)), Dot(B, Jacobian(v))))>>::Value);
    static_assert(!FormLanguage::IsSpecialized<
                  QuadratureRule<decltype(Dot(Jacobian(u) * velocity, v))>>::Value);
    for (auto cell = trialFES.getMesh().getCell(); cell; ++cell)
    {
      checkEntries(Dot(Jacobian(u), Jacobian(v)), *cell);
      checkEntries(Dot(weight * Jacobian(u), Jacobian(v)), *cell);
      checkEntries(Dot(Dot(A, Jacobian(u)), Dot(B, Jacobian(v))), *cell);
      checkEntries(Dot(Dot(A, Jacobian(u)), Dot(A, Jacobian(v))), *cell);
      checkEntries(Dot(Jacobian(u) * velocity, v), *cell);
    }
  }

  template <class Scalar>
  void checkFamilies(const Mesh<>& mesh, size_t rows, size_t cols)
  {
    using Matrix = Math::SpatialMatrix<Scalar>;
    P0<Matrix> p0(mesh, rows, cols);
    P0g<Matrix, Mesh<>> p0g(mesh, rows, cols);
    P1<Matrix> p1(mesh, rows, cols);
    H1<1, Matrix> h1(std::integral_constant<size_t, 1>{}, mesh, rows, cols);
    H1<2, Matrix> h2(std::integral_constant<size_t, 2>{}, mesh, rows, cols);
    H1<3, Matrix> h3(std::integral_constant<size_t, 3>{}, mesh, rows, cols);
    checkJacobians(p0, p0);
    checkJacobians(p0g, p0g);
    checkJacobians(p1, p1);
    checkJacobians(h1, h1);
    checkJacobians(h2, h2);
    checkJacobians(h3, h3);
    checkJacobians(h1, h2);
    checkJacobians(h3, h1);
  }
}

TEST(MatrixQuadrature, IntegralRequiresScalarRange)
{
  auto mesh = makeMesh(G::Triangle);
  auto check = [&]<class FES>(const FES& fes) {
    TestFunction v(fes);
    static_assert(!ScalarIntegral<decltype(v)>);
    TrialFunction u(fes);
    static_assert(requires { Integral(u, v); });
    static_assert(ScalarIntegral<decltype(Dot(u, v))>);
  };
  P1 vectorP1(mesh, 2);
  P1 matrixP1(mesh, 2, 3);
  P0<Math::SpatialVector<Real>> vectorP0(mesh, 2);
  P0 matrixP0(mesh, 2, 3);
  P0g<Math::SpatialVector<Real>, Mesh<>> vectorP0g(mesh, 2);
  P0g matrixP0g(mesh, 2, 3);
  H1 vectorH1(std::integral_constant<size_t, 2>{}, mesh, 2);
  H1 matrixH1(std::integral_constant<size_t, 2>{}, mesh, 2, 3);
  check(vectorP1);
  check(matrixP1);
  check(vectorP0);
  check(matrixP0);
  check(vectorP0g);
  check(matrixP0g);
  check(vectorH1);
  check(matrixH1);
  P1 scalar(mesh);
  TestFunction v(scalar);
  static_assert(ScalarIntegral<decltype(v)>);
}

namespace
{
  template <class Scalar>
  void checkPotential(const Mesh<>& mesh, size_t rows, size_t cols)
  {
    using Matrix = Math::SpatialMatrix<Scalar>;
    P1<Matrix> fes(mesh, rows, cols);
    P1<Scalar> scalar(mesh);
    TrialFunction u(fes);
    TrialFunction su(scalar);
    TestFunction v(fes);
    TestFunction sv(scalar);
    const auto scalarKernel = [](const Point&, const Point&) { return Scalar(1); };
    static_assert(std::is_same_v<
      typename FormLanguage::Traits<decltype(Potential(scalarKernel, u))>::ScalarType,
      Scalar>);
    auto scalarRule = Integral(Potential(scalarKernel, su), sv);
    auto matrixRule = Integral(Potential(scalarKernel, u), v);
    Matrix coupling(rows, rows);
    for (size_t i = 0; i < rows; ++i)
      for (size_t j = 0; j < rows; ++j)
        coupling(i, j) = Scalar(1 + 2 * i + 3 * j);
    const auto rowKernel = [&](
                             Matrix& out, const Point&, const Point&) { out = coupling; };
    auto rowRule = Integral(Potential(rowKernel, u), v);
    Math::SpatialTensor<Scalar, 4> tensor(rows, cols, rows, cols);
    for (size_t i = 0; i < rows; ++i)
      for (size_t j = 0; j < cols; ++j)
        for (size_t k = 0; k < rows; ++k)
          for (size_t l = 0; l < cols; ++l)
          {
            tensor(i, j, k, l) = Scalar(1 + i + 2 * j + 3 * k + 4 * l);
            if constexpr (std::is_same_v<Scalar, Complex>)
              tensor(i, j, k, l) += Complex(0, 2 + i + j + k + l);
          }
    const auto tensorKernel = [&](Math::SpatialTensor<Scalar, 4>& out, const Point&,
                                const Point&) { out = tensor; };
    auto tensorRule = Integral(Potential(tensorKernel, u), v);
    // Test coincident and distinct cell pairs, then rebind to a coincident pair.
    const auto first = *mesh.getCell();
    Index lastIndex = first.getIndex();
    for (auto cell = mesh.getCell(); cell; ++cell)
      lastIndex = cell->getIndex();
    const auto last = *mesh.getCell(lastIndex);
    const size_t components = rows * cols;
    for (const auto& test : {first, last, first})
    {
      scalarRule.setPolytope(first, test);
      matrixRule.setPolytope(first, test);
      rowRule.setPolytope(first, test);
      tensorRule.setPolytope(first, test);
      const size_t ntr =
        scalar.getFiniteElement(first.getDimension(), first.getIndex()).getCount();
      const size_t nte =
        scalar.getFiniteElement(test.getDimension(), test.getIndex()).getCount();
      for (size_t tr = 0; tr < ntr * components; ++tr)
        for (size_t te = 0; te < nte * components; ++te)
        {
          const auto block = scalarRule.integrate(tr / components, te / components);
          const size_t tc = te % components, rc = tr % components;
          EXPECT_LE(
            std::abs(matrixRule.integrate(tr, te) - (tc == rc ? block : Scalar(0))),
            tolerance);
          EXPECT_LE(std::abs(rowRule.integrate(tr, te) -
                      (tc % cols == rc % cols ? coupling(tc / cols, rc / cols) * block
                                              : Scalar(0))),
            tolerance);
          EXPECT_LE(std::abs(tensorRule.integrate(tr, te) -
                      tensor(tc / cols, tc % cols, rc / cols, rc % cols) * block),
            tolerance);
        }
    }
    // Verify global nonlocal assembly uses the space's actual DOF maps.
    // Independent reference: the constant kernel separates the two integrals.
    // Non-triangle panels use the existing centroid rule; coincident affine
    // triangles integrate their linear bases exactly with the collapsed rule.
    Math::Vector<Scalar> moments = Math::Vector<Scalar>::Zero(scalar.getSize());
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const QF::Centroid formula(cell->getGeometry());
      const Point point(*cell, formula.getPoint(0));
      const auto& fe = scalar.getFiniteElement(cell->getDimension(), cell->getIndex());
      const auto& dofs = scalar.getDOFs(cell->getDimension(), cell->getIndex());
      for (size_t local = 0; local < fe.getCount(); ++local)
        moments(dofs(local)) += formula.getWeight(0) * point.getDistortion() *
          fe.getBasis(local)(formula.getPoint(0));
    }
    for (size_t action = 0; action < 3; ++action)
    {
      SCOPED_TRACE(::testing::Message() << "kernel action " << action);
      DenseProblem problem(u, v);
      if (action == 0)
        problem = Integral(Potential(scalarKernel, u), v);
      else if (action == 1)
        problem = Integral(Potential(rowKernel, u), v);
      else
        problem = Integral(Potential(tensorKernel, u), v);
      problem.assemble();
      const auto& actual = problem.getLinearSystem().getOperator();
      EXPECT_GT(actual.norm(), 0);
      for (size_t te = 0; te < fes.getSize(); ++te)
        for (size_t tr = 0; tr < fes.getSize(); ++tr)
        {
          const size_t tc = te % components, rc = tr % components;
          const Scalar coefficient = action == 0 ? Scalar(tc == rc)
            : action == 1
            ? (tc % cols == rc % cols ? coupling(tc / cols, rc / cols) : Scalar(0))
            : tensor(tc / cols, tc % cols, rc / cols, rc % cols);
          const Scalar expected =
            Math::dot(moments(te / components), moments(tr / components)) * coefficient;
          EXPECT_LE(
            std::abs(actual(te, tr) - expected), tolerance * (1 + std::abs(expected)))
            << "entry " << te << ", " << tr << ": " << actual(te, tr) << " vs "
            << expected;
        }
    }
  }
}

TEST(MatrixQuadrature, MatrixPotentialAllGeometriesAndShapes)
{
  for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
         G::Hexahedron, G::Wedge, G::Pyramid})
  {
    auto mesh = makeMesh(geometry);
    for (size_t rows = 1; rows <= 3; ++rows)
      for (size_t cols = 1; cols <= 3; ++cols)
      {
        SCOPED_TRACE(
          ::testing::Message() << int(geometry) << ": " << rows << "x" << cols);
        checkPotential<Real>(mesh, rows, cols);
      }
  }
}

TEST(MatrixQuadrature, ComplexMatrixPotential)
{
  const auto mesh = makeMesh(G::Triangle);
  checkPotential<Complex>(mesh, 2, 3);
}

TEST(MatrixQuadrature, JacobianFormsAllSpacesGeometriesAndShapes)
{
  for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
         G::Hexahedron, G::Wedge, G::Pyramid})
  {
    const auto mesh = makeMesh(geometry);
    for (size_t rows = 1; rows <= 3; ++rows)
      for (size_t cols = 1; cols <= 3; ++cols)
      {
        SCOPED_TRACE(
          ::testing::Message() << int(geometry) << ": " << rows << "x" << cols);
        checkFamilies<Real>(mesh, rows, cols);
      }
  }
}

TEST(MatrixQuadrature, ComplexJacobianForms)
{
  for (auto geometry : {G::Triangle, G::Tetrahedron})
  {
    const auto mesh = makeMesh(geometry);
    checkFamilies<Complex>(mesh, 2, 3);
  }
}

TEST(MatrixQuadrature, CurvedEmbeddedAndPointJacobians)
{
  auto builder = Mesh<>::Builder();
  builder.initialize(3).nodes(3);
  builder.vertex({0, 0, 0}).vertex({1, 0, 1}).vertex({0, 1, 1});
  auto mesh = builder.polytope(G::Triangle, {0, 1, 2}).finalize();
  mesh.getConnectivity().compute(2, 1);
  checkFamilies<Real>(mesh, 2, 3);
  RealH1Element<2> geometryElement(G::Triangle);
  PointCloud nodes(3, geometryElement.getCount());
  for (size_t a = 0; a < geometryElement.getCount(); ++a)
  {
    const auto& r = geometryElement.getNode(a);
    nodes(0, a) = r[0] + 0.1 * r[0] * r[1];
    nodes(1, a) = r[1] + 0.2 * r[0] * r[1];
    nodes(2, a) = r[0] + r[1] + 0.3 * r[0] * r[1];
  }
  mesh.setPolytopeTransformation(
    {2, 0}, new ParametricTransformation<RealH1Element<2>>(nodes, geometryElement));
  checkFamilies<Real>(mesh, 2, 3);
  auto pointBuilder = Mesh<>::Builder();
  pointBuilder.initialize(3).nodes(1).vertex({0, 0, 0});
  const auto pointMesh = pointBuilder.finalize();
  checkFamilies<Real>(pointMesh, 2, 3);
}
