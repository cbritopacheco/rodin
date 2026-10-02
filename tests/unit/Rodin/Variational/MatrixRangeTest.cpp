/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Matrix finite element, tensor operator, and backend regressions. */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Assembly.h"
#include "Rodin/Geometry/ParametricTransformation.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  std::vector<std::pair<size_t, size_t>> matrixShapes()
  {
    std::vector<std::pair<size_t, size_t>> shapes;
    for (size_t rows = 1; rows <= 3; ++rows)
      for (size_t cols = 1; cols <= 3; ++cols)
        shapes.emplace_back(rows, cols);
    return shapes;
  }

  using Matrix = Math::SpatialMatrix<Real>;
  using G = Polytope::Type;

  template <class Element>
  void checkElement(G geometry, size_t rows, size_t cols)
  {
    Element fe(geometry, rows, cols);
    const auto& scalar = fe.getScalarElement();
    ASSERT_EQ(fe.getCount(), scalar.getCount() * rows * cols);
    EXPECT_EQ(fe.getOrder(), scalar.getOrder());
    for (size_t a = 0; a < fe.getCount(); ++a)
      for (size_t b = 0; b < fe.getCount(); ++b)
        EXPECT_NEAR(fe.getLinearForm(a)(fe.getBasis(b)), a == b ? 1.0 : 0.0, 1e-10);
  }

  template <class FES, class ScalarFES>
  void checkField(const FES& fes, const ScalarFES& scalar)
  {
    const size_t rows = fes.getRows(), cols = fes.getColumns(), components = rows * cols;
    EXPECT_EQ(fes.getSize(), scalar.getSize() * components);
    Matrix value(rows, cols);
    for (size_t r = 0; r < rows; ++r)
      for (size_t c = 0; c < cols; ++c)
        value(r, c) = 1 + 10 * r + c;
    GridFunction field(fes);
    field = value;
    const FES copied(fes);
    const FES moved(std::move(FES(copied)));
    FES assigned(copied);
    assigned = fes;
    EXPECT_EQ(moved.getSize(), fes.getSize());
    EXPECT_EQ(field.getRows(), rows);
    EXPECT_EQ(field.getColumns(), cols);
    TrialFunction trial(fes);
    TestFunction test(fes);
    TrialFunction scalarTrial(scalar);
    TestFunction scalarTest(scalar);
    auto mass = Integral(trial, test);
    auto scalarMass = Integral(scalarTrial, scalarTest);
    for (auto it = fes.getMesh().getCell(); it; ++it)
    {
      const auto& poly = *it;
      const auto& fe = fes.getFiniteElement(poly.getDimension(), poly.getIndex());
      Geometry::Point point(poly, Polytope::Traits(poly.getGeometry()).getCentroid());
      const auto evaluated = field.getValue(point);
      ASSERT_EQ(evaluated.rows(), rows);
      ASSERT_EQ(evaluated.cols(), cols);
      for (size_t r = 0; r < rows; ++r)
        for (size_t c = 0; c < cols; ++c)
          EXPECT_NEAR(evaluated(r, c), value(r, c), 1e-10);
      const auto& dofs = fes.getDOFs(poly.getDimension(), poly.getIndex());
      EXPECT_TRUE((moved.getDOFs(poly.getDimension(), poly.getIndex()) == dofs).all());
      EXPECT_TRUE((assigned.getDOFs(poly.getDimension(), poly.getIndex()) == dofs).all());
      const auto& scalarDOFs = scalar.getDOFs(poly.getDimension(), poly.getIndex());
      ASSERT_EQ(dofs.size(), fe.getCount());
      for (size_t a = 0; a < fe.getCount(); ++a)
      {
        EXPECT_EQ(dofs[a], scalarDOFs[a / components] * components + a % components);
        EXPECT_EQ(fes.getGlobalIndex({poly.getDimension(), poly.getIndex()}, a), dofs[a]);
      }
      mass.setPolytope(poly);
      scalarMass.setPolytope(poly);
      for (size_t a = 0; a < fe.getCount(); ++a)
        for (size_t b = 0; b < fe.getCount(); ++b)
          EXPECT_NEAR(mass.integrate(a, b),
            a % components == b % components
              ? scalarMass.integrate(a / components, b / components)
              : 0.0,
            1e-10);
      auto load = Integral(MatrixFunction(value), test);
      load.setPolytope(poly);
      // Summing the load entries tests the partition of unity for each component.
      for (size_t c = 0; c < components; ++c)
      {
        Real sum = 0;
        for (size_t a = c; a < fe.getCount(); a += components)
          sum += load.integrate(a);
        EXPECT_NEAR(sum, value(c / cols, c % cols) * poly.getMeasure(), 1e-9);
      }
    }
    BilinearForm globalMass(trial, test);
    BilinearForm globalScalarMass(scalarTrial, scalarTest);
    globalMass = Integral(trial, test);
    globalScalarMass = Integral(scalarTrial, scalarTest);
    globalMass.assemble();
    globalScalarMass.assemble();
    Math::SparseMatrix<Real> sequential;
    Assembly::Sequential<Math::SparseMatrix<Real>, decltype(globalMass)> assembler;
    assembler.execute(sequential,
      {fes, fes, globalMass.getLocalIntegrators(), globalMass.getGlobalIntegrators()});
    EXPECT_LE((sequential - globalMass.getOperator()).norm(), 1e-10);
    using DenseForm =
      BilinearForm<typename decltype(trial)::SolutionType, FES, FES, Math::Matrix<Real>>;
    DenseForm dense(trial, test);
    dense = Integral(trial, test);
    dense.assemble();
    EXPECT_LE((dense.getOperator() - Math::Matrix<Real>(sequential)).norm(), 1e-10);
    const auto& actual = globalMass.getOperator();
    const auto& expected = globalScalarMass.getOperator();
    for (size_t a = 0; a < fes.getSize(); ++a)
      for (size_t b = 0; b < fes.getSize(); ++b)
        EXPECT_NEAR(actual.coeff(a, b),
          a % components == b % components
            ? expected.coeff(a / components, b / components)
            : 0.0,
          1e-10);
  }
}

TEST(MatrixRange, NodalDualityAllGeometriesAndShapes)
{
  for (auto geometry : {G::Point, G::Segment, G::Triangle, G::Quadrilateral,
         G::Tetrahedron, G::Pyramid, G::Wedge, G::Hexahedron})
    for (size_t rows = 1; rows <= 3; ++rows)
      for (size_t cols = 1; cols <= 3; ++cols)
      {
        checkElement<P0Element<Matrix>>(geometry, rows, cols);
        checkElement<P0gElement<Matrix>>(geometry, rows, cols);
        checkElement<P1Element<Matrix>>(geometry, rows, cols);
        checkElement<H1Element<1, Matrix>>(geometry, rows, cols);
        checkElement<H1Element<2, Matrix>>(geometry, rows, cols);
        checkElement<H1Element<3, Matrix>>(geometry, rows, cols);
      }
}

TEST(MatrixRange, ProjectionAndAssemblyAllSpaces)
{
  for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
         G::Pyramid, G::Wedge, G::Hexahedron})
  {
    const bool threeD = Polytope::Traits(geometry).getDimension() == 3;
    Mesh mesh = geometry == G::Segment ? LocalMesh::UniformGrid(geometry, {2})
      : threeD                         ? LocalMesh::UniformGrid(geometry, {2, 2, 2})
                                       : LocalMesh::UniformGrid(geometry, {2, 2});
    for (size_t d = 1; d <= mesh.getDimension(); ++d)
      for (size_t lower = 0; lower < d; ++lower)
        mesh.getConnectivity().compute(d, lower);
    for (const auto& shape : matrixShapes())
    {
      const auto [rows, cols] = shape;
      P0 scalarP0(mesh);
      P0 matrixP0(mesh, rows, cols);
      P0g scalarP0g(mesh);
      P0g matrixP0g(mesh, rows, cols);
      P1 scalarP1(mesh);
      P1 matrixP1(mesh, rows, cols);
      H1 scalarH1(std::integral_constant<size_t, 1>{}, mesh);
      H1 matrixH1(std::integral_constant<size_t, 1>{}, mesh, rows, cols);
      H1 scalarH2(std::integral_constant<size_t, 2>{}, mesh);
      H1 matrixH2(std::integral_constant<size_t, 2>{}, mesh, rows, cols);
      H1 scalarH3(std::integral_constant<size_t, 3>{}, mesh);
      H1 matrixH3(std::integral_constant<size_t, 3>{}, mesh, rows, cols);
      checkField(matrixP0, scalarP0);
      checkField(matrixP0g, scalarP0g);
      checkField(matrixP1, scalarP1);
      checkField(matrixH1, scalarH1);
      checkField(matrixH2, scalarH2);
      checkField(matrixH3, scalarH3);
    }
  }
}

TEST(MatrixRange, RejectsInvalidDimensions)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  EXPECT_THROW(P1(mesh, 0, 3), std::exception);
  EXPECT_THROW(P0(mesh, 3, 4), std::exception);
  EXPECT_THROW(P0g(mesh, 4, 3), std::exception);
  EXPECT_THROW((H1(std::integral_constant<size_t, 3>{}, mesh, 3, 0)), std::exception);
}

TEST(MatrixRange, ComponentAndTransposeForms)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  P1 fes(mesh, 2, 3);
  TrialFunction trial(fes);
  TestFunction test(fes);
  auto componentMass = Integral(Component(trial, 1, 2), Component(test, 1, 2));
  auto transposeMass = Integral(Transpose(trial), Transpose(test));
  auto mass = Integral(trial, test);
  for (auto it = mesh.getCell(); it; ++it)
  {
    componentMass.setPolytope(*it);
    transposeMass.setPolytope(*it);
    mass.setPolytope(*it);
    const size_t count = fes.getFiniteElement(2, it->getIndex()).getCount();
    for (size_t a = 0; a < count; ++a)
      for (size_t b = 0; b < count; ++b)
      {
        EXPECT_NEAR(transposeMass.integrate(a, b), mass.integrate(a, b), 1e-12);
        EXPECT_NEAR(componentMass.integrate(a, b),
          a % 6 == 5 && b % 6 == 5 ? mass.integrate(a, b) : 0.0, 1e-12);
      }
  }
}

TEST(MatrixRange, CubicProjection)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {3, 3});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  H1 fes(std::integral_constant<size_t, 3>{}, mesh, 3, 3);
  GridFunction field(fes);
  auto polynomial = [](const Geometry::Point& p) {
    Matrix value(3, 3);
    for (size_t r = 0; r < 3; ++r)
      for (size_t c = 0; c < 3; ++c)
        value(r, c) = (r + 1) * p.x() * p.x() * p.y() + (c + 1) * p.y() + 10 * r + c;
    return value;
  };
  field.project(polynomial);
  for (auto it = mesh.getCell(); it; ++it)
  {
    Geometry::Point p(*it, Math::SpatialPoint{0.2, 0.3});
    const auto actual = field.getValue(p), expected = polynomial(p);
    for (size_t r = 0; r < 3; ++r)
      for (size_t c = 0; c < 3; ++c)
        EXPECT_NEAR(actual(r, c), expected(r, c), 1e-10);
    EXPECT_NEAR(Component(field, 1, 2).getValue(p), expected(1, 2), 1e-10);
  }
}

TEST(MatrixRange, ComplexFrobeniusAndSpaces)
{
  Math::SpatialMatrix<Complex> a(2, 3), b(2, 3);
  a.setZero();
  b.setZero();
  a(1, 2) = Complex(2, 3);
  b(1, 2) = Complex(4, 5);
  const Complex expected = a(1, 2) * std::conj(b(1, 2));
  const Math::Matrix<Complex> eigenA = a.getData().topLeftCorner(2, 3);
  const Math::Matrix<Complex> eigenB = b.getData().topLeftCorner(2, 3);
  EXPECT_EQ(Math::dot(a, b), expected);
  EXPECT_EQ(Math::dot(eigenA, b), expected);
  EXPECT_EQ(Math::dot(a, eigenB), expected);
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  MatrixP1<LocalMesh> realSpace(mesh, 2, 3);
  P1<Math::SpatialMatrix<Complex>, LocalMesh> complexSpace(mesh, 2, 3);
  H1<3, Math::SpatialMatrix<Complex>, LocalMesh> higherSpace(
    std::integral_constant<size_t, 3>{}, mesh, 2, 3);
  EXPECT_EQ(complexSpace.getSize(), realSpace.getSize());
  GridFunction field(higherSpace);
  field = a;
  auto cell = mesh.getCell();
  Geometry::Point p(*cell, Math::SpatialPoint{0.2, 0.3});
  const auto actual = field.getValue(p);
  EXPECT_NEAR(std::abs(actual(1, 2) - a(1, 2)), 0, 1e-10);
  auto checkMass = [&](const auto& fes) {
    TrialFunction trial(fes);
    TestFunction test(fes);
    auto mass = Integral(trial, test);
    auto load = Integral(MatrixFunction(a), test);
    mass.setPolytope(*cell);
    load.setPolytope(*cell);
    const size_t count = fes.getFiniteElement(2, cell->getIndex()).getCount();
    for (size_t b = 0; b < count; ++b)
    {
      Complex product = 0;
      for (size_t i = 0; i < count; ++i)
        product += mass.integrate(i, b) * a((i % 6) / 3, i % 3);
      EXPECT_NEAR(std::abs(product - load.integrate(b)), 0, 1e-10);
    }
  };
  checkMass(complexSpace);
  checkMass(higherSpace);
  using ComplexMatrix = Math::SpatialMatrix<Complex>;
  checkMass(P0<ComplexMatrix>(mesh, 2, 3));
  checkMass(P0g<ComplexMatrix, LocalMesh>(mesh, 2, 3));
  checkMass(H1<1, ComplexMatrix>(std::integral_constant<size_t, 1>{}, mesh, 2, 3));
  checkMass(H1<2, ComplexMatrix>(std::integral_constant<size_t, 2>{}, mesh, 2, 3));
}

TEST(MatrixRange, EvaluationCacheSurvivesAddressReuse)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {2, 2});
  using FES = P0<Matrix, LocalMesh>;
  std::optional<FES> fes;
  std::optional<GridFunction<FES, Math::Vector<Real>>> field;
  auto cell = mesh.getCell();
  ++cell;
  Geometry::Point p(*cell, Math::SpatialPoint{0.2, 0.3});
  for (size_t cols : {3u, 2u})
  {
    field.reset();
    fes.emplace(mesh, 3, cols);
    field.emplace(*fes);
    Matrix value(3, cols);
    value.setZero();
    value(2, cols - 1) = 7;
    *field = value;
    const auto actual = field->getValue(p);
    EXPECT_EQ(actual.rows(), 3);
    EXPECT_EQ(actual.cols(), cols);
    EXPECT_EQ(actual(2, cols - 1), 7);
  }
}

TEST(MatrixRange, MatrixDirichletValues)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {2, 2});
  mesh.getConnectivity().compute(1, 2);
  P1 fes(mesh, 3, 3);
  TrialFunction trial(fes);
  Matrix value(3, 3);
  for (size_t r = 0; r < 3; ++r)
    for (size_t c = 0; c < 3; ++c)
      value(r, c) = 1 + 10 * r + c;
  DirichletBC boundary(trial, MatrixFunction(value));
  boundary.assemble();
  ASSERT_EQ(
    std::get<typename decltype(boundary)::ValueDOFs>(boundary.getDOFs()).size(), 4 * 9);
  for (const auto& [global, prescribed] :
    std::get<typename decltype(boundary)::ValueDOFs>(boundary.getDOFs()))
    EXPECT_EQ(prescribed, value((global % 9) / 3, global % 3));
}

TEST(MatrixRange, ScalarAndMatrixWeightedMass)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  auto check = [&](const auto& fes) {
    TrialFunction trial(fes);
    TestFunction test(fes);
    Matrix coefficient(2, 2);
    coefficient(0, 0) = 2;
    coefficient(0, 1) = 3;
    coefficient(1, 0) = 5;
    coefficient(1, 1) = 7;
    auto plain = Integral(trial, test);
    auto scaled = Integral(RealFunction(4) * trial, test);
    auto weighted = Integral(MatrixFunction(coefficient) * trial, test);
    for (auto it = mesh.getCell(); it; ++it)
    {
      plain.setPolytope(*it);
      scaled.setPolytope(*it);
      weighted.setPolytope(*it);
      const size_t count = fes.getFiniteElement(2, it->getIndex()).getCount();
      for (size_t a = 0; a < count; ++a)
        for (size_t b = 0; b < count; ++b)
        {
          EXPECT_NEAR(scaled.integrate(a, b), 4 * plain.integrate(a, b), 1e-10);
          const size_t ca = a % 6, cb = b % 6;
          const Real expected = ca % 3 == cb % 3
            ? coefficient(cb / 3, ca / 3) * plain.integrate((a / 6) * 6, (b / 6) * 6)
            : 0;
          EXPECT_NEAR(weighted.integrate(a, b), expected, 1e-10);
        }
    }
  };
  P1 p1(mesh, 2, 3);
  check(p1);
  H1 h1(std::integral_constant<size_t, 1>{}, mesh, 2, 3);
  check(h1);
  H1 h2(std::integral_constant<size_t, 2>{}, mesh, 2, 3);
  check(h2);
  H1 h3(std::integral_constant<size_t, 3>{}, mesh, 2, 3);
  check(h3);
}

TEST(MatrixRange, TensorArithmeticAndContraction)
{
  Math::SpatialTensor<Real> a(2, 3, 2);
  a.setZero();
  a(1, 2, 1) = 5;
  EXPECT_EQ(a.size(), 12);
  EXPECT_EQ(a[11], 5);
  EXPECT_EQ(a.norm(), 5);
  EXPECT_EQ(Math::dot(a, a), 25);
  Math::SpatialTensor<Complex> b(2, 3, 2);
  b.setZero();
  b(1, 2, 1) = Complex(2, 3);
  EXPECT_EQ(Math::dot(a, b), 5.0 * std::conj(Complex(2, 3)));
  EXPECT_EQ(Math::dot(b, a), Complex(10, 15));
  EXPECT_NEAR((b + b - 2 * b).norm(), 0, 1e-12);
  EXPECT_EQ((a / 5)(1, 2, 1), 1);
  Math::SpatialVector<Real> direction{2, 3};
  const auto contraction = a * direction;
  EXPECT_EQ(contraction.rows(), 2);
  EXPECT_EQ(contraction.cols(), 3);
  EXPECT_EQ(contraction(1, 2), 15);
  Math::SpatialTensor<Real, 4> fourth(2, 3, 2, 3);
  fourth.setConstant(1);
  EXPECT_EQ(fourth.size(), 36);
  EXPECT_EQ(Math::dot(fourth, fourth), 36);
  EXPECT_THROW((Math::SpatialTensor<Real>(4, 2, 2)), std::exception);
}

TEST(MatrixRange, MatrixDifferentialOperatorsAllSpacesAndGeometries)
{
  constexpr Real step = 1e-6;
  constexpr Real tolerance = 2e-7;
  for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
         G::Pyramid, G::Wedge, G::Hexahedron})
  {
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    Mesh mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {2})
      : dimension == 2         ? LocalMesh::UniformGrid(geometry, {2, 2})
                               : LocalMesh::UniformGrid(geometry, {2, 2, 2});
    for (size_t d = 1; d <= dimension; ++d)
      for (size_t lower = 0; lower < d; ++lower)
        mesh.getConnectivity().compute(d, lower);
    auto check = [&](const auto& fes, bool constant) {
      GridFunction field(fes);
      auto affine = [&](const Geometry::Point& point) {
        Matrix value(fes.getRows(), fes.getColumns());
        for (size_t row = 0; row < value.rows(); ++row)
          for (size_t col = 0; col < value.cols(); ++col)
          {
            value(row, col) = 1 + row + col;
            if (!constant)
              for (size_t k = 0; k < dimension; ++k)
                value(row, col) +=
                  (row + 1) * (col + 2) * (k + 1) * point.getPhysicalCoordinates()(k);
          }
        return value;
      };
      field.project(affine);
      auto gradient = Grad(field);
      auto jacobian = Jacobian(field);
      TrialFunction trial(fes);
      TestFunction test(fes);
      auto shapeGradient = Grad(trial);
      auto stiffness = Integral(Grad(trial), Grad(test));
      auto cell = mesh.getCell();
      Geometry::Point point(*cell, Polytope::Traits(geometry).getCentroid());
      const IntegrationPoint ip(point);
      shapeGradient.setIntegrationPoint(ip);
      const auto evaluated = gradient.getValue(ip);
      EXPECT_NEAR((jacobian.getValue(ip) - evaluated).norm(), 0, tolerance);
      EXPECT_NEAR(Frobenius(gradient).getValue(ip), evaluated.norm(), tolerance);
      for (size_t k = 0; k < dimension; ++k)
      {
        const auto derivative = Derivative(k, field).getValue(ip);
        auto partial = Derivative(k, trial);
        partial.setIntegrationPoint(ip);
        auto physicalPlus = point.getPhysicalCoordinates(), physicalMinus = physicalPlus;
        physicalPlus(k) += step;
        physicalMinus(k) -= step;
        Math::SpatialPoint plus, minus;
        cell->getTransformation().inverse(plus, physicalPlus);
        cell->getTransformation().inverse(minus, physicalMinus);
        const auto finiteDifference = (field.getValue(Geometry::Point(*cell, plus)) -
                                        field.getValue(Geometry::Point(*cell, minus))) /
          (2 * step);
        for (size_t row = 0; row < fes.getRows(); ++row)
          for (size_t col = 0; col < fes.getColumns(); ++col)
          {
            const Real expected = constant ? 0 : (row + 1) * (col + 2) * (k + 1);
            EXPECT_NEAR(evaluated(row, col, k), expected, tolerance);
            EXPECT_NEAR(derivative(row, col), expected, tolerance);
            EXPECT_NEAR(finiteDifference(row, col), expected, tolerance);
          }
        EXPECT_NEAR(
          partial.getBasis(0)(0, 0), shapeGradient.getBasis(0)(0, 0, k), tolerance);
      }
      if (fes.getColumns() == dimension)
      {
        const auto divergence = Div(field).getValue(ip);
        auto basisDivergence = Div(test);
        basisDivergence.setIntegrationPoint(ip);
        for (size_t row = 0; row < fes.getRows(); ++row)
        {
          Real expected = 0;
          for (size_t k = 0; k < dimension; ++k)
            expected += evaluated(row, k, k);
          EXPECT_NEAR(divergence(row), expected, tolerance);
        }
      }
      stiffness.setPolytope(*cell);
      const auto count = fes.getFiniteElement(dimension, cell->getIndex()).getCount();
      Real energy = 0;
      const auto& dofs = fes.getDOFs(dimension, cell->getIndex());
      for (size_t a = 0; a < count; ++a)
        for (size_t b = 0; b < count; ++b)
          energy += field[dofs[a]] * stiffness.integrate(a, b) * field[dofs[b]];
      EXPECT_NEAR(energy, evaluated.squaredNorm() * cell->getMeasure(), 1e-6);
    };
    for (const auto shape : matrixShapes())
    {
      const auto [rows, cols] = shape;
      P0 p0(mesh, rows, cols);
      check(p0, true);
      P0g p0g(mesh, rows, cols);
      check(p0g, true);
      P1 p1(mesh, rows, cols);
      check(p1, false);
      H1 h1(std::integral_constant<size_t, 1>{}, mesh, rows, cols);
      check(h1, false);
      H1 h2(std::integral_constant<size_t, 2>{}, mesh, rows, cols);
      check(h2, false);
      H1 h3(std::integral_constant<size_t, 3>{}, mesh, rows, cols);
      check(h3, false);
    }
  }
}

TEST(MatrixRange, TensorCoefficientCouplesAllMatrixEntries)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  Math::SpatialTensor<Real, 4> material(2, 3, 2, 3);
  for (size_t a = 0; a < material.size(); ++a)
    material[a] = 1 + a;
  auto check = [&](const auto& fes) {
    TrialFunction trial(fes);
    TestFunction test(fes);
    auto weighted = Integral(TensorFunction(material) * trial, test);
    auto mass = Integral(trial, test);
    auto cell = mesh.getCell();
    weighted.setPolytope(*cell);
    mass.setPolytope(*cell);
    const size_t count = fes.getFiniteElement(2, cell->getIndex()).getCount();
    for (size_t a = 0; a < count; ++a)
      for (size_t b = 0; b < count; ++b)
      {
        const size_t ca = a % 6, cb = b % 6;
        EXPECT_NEAR(weighted.integrate(a, b),
          material(cb / 3, cb % 3, ca / 3, ca % 3) * mass.integrate(a / 6 * 6, b / 6 * 6),
          1e-10);
      }
  };
  P0 p0(mesh, 2, 3);
  check(p0);
  P0g p0g(mesh, 2, 3);
  check(p0g);
  P1 p1(mesh, 2, 3);
  check(p1);
  H1 h1(std::integral_constant<size_t, 1>{}, mesh, 2, 3);
  check(h1);
  H1 h2(std::integral_constant<size_t, 2>{}, mesh, 2, 3);
  check(h2);
  H1 h3(std::integral_constant<size_t, 3>{}, mesh, 2, 3);
  check(h3);
}

TEST(MatrixRange, HigherOrdersPreserveComponentDuality)
{
  auto check = []<size_t K>(std::integral_constant<size_t, K>, G geometry) {
    H1Element<K, Matrix> fe(geometry, 2, 3);
    for (size_t a = 0; a < fe.getCount(); ++a)
    {
      EXPECT_NEAR(fe.getLinearForm(a)(fe.getBasis(a)), 1, 1e-8);
      EXPECT_NEAR(fe.getLinearForm(a)(fe.getBasis((a + 1) % fe.getCount())), 0, 1e-8);
    }
  };
  for (auto geometry : {G::Point, G::Segment, G::Triangle, G::Quadrilateral,
         G::Tetrahedron, G::Pyramid, G::Wedge, G::Hexahedron})
  {
    check(std::integral_constant<size_t, 4>{}, geometry);
    check(std::integral_constant<size_t, 5>{}, geometry);
  }
}

TEST(MatrixRange, ComplexMatrixDifferentials)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  using ComplexMatrix = Math::SpatialMatrix<Complex>;
  auto check = [&](const auto& fes) {
    GridFunction field(fes);
    field = [](const Geometry::Point& p) {
      ComplexMatrix value(2, 2);
      for (size_t r = 0; r < 2; ++r)
        for (size_t c = 0; c < 2; ++c)
          value(r, c) = Complex((1 + r) * p.x(), (1 + c) * p.y());
      return value;
    };
    auto cell = mesh.getCell();
    Geometry::Point point(*cell, Polytope::Traits(G::Triangle).getCentroid());
    const auto gradient = Grad(field).getValue(point);
    for (size_t r = 0; r < 2; ++r)
      for (size_t c = 0; c < 2; ++c)
      {
        EXPECT_NEAR(std::abs(gradient(r, c, 0) - Complex(1 + r, 0)), 0, 1e-10);
        EXPECT_NEAR(std::abs(gradient(r, c, 1) - Complex(0, 1 + c)), 0, 1e-10);
      }
    EXPECT_NEAR(gradient.squaredNorm(), 20, 1e-10);
    const auto value = field.getValue(point);
    EXPECT_NEAR(value.squaredNorm(), Math::dot(value, value).real(), 1e-10);
    EXPECT_NEAR(value.norm() * value.norm(), value.squaredNorm(), 1e-10);
    EXPECT_NEAR(std::abs(Component(field, 1, 0).getValue(point) - value(1, 0)), 0, 1e-10);
  };
  check(P1<ComplexMatrix>(mesh, 2, 2));
  check(H1<1, ComplexMatrix>(std::integral_constant<size_t, 1>{}, mesh, 2, 2));
  check(H1<2, ComplexMatrix>(std::integral_constant<size_t, 2>{}, mesh, 2, 2));
  check(H1<3, ComplexMatrix>(std::integral_constant<size_t, 3>{}, mesh, 2, 2));
}

TEST(MatrixRange, BoundaryGradientsMatchScalarTraceSemantics)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {3, 3});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(1, 0);
  H1 fes(std::integral_constant<size_t, 3>{}, mesh, 2, 3);
  H1 scalar(std::integral_constant<size_t, 3>{}, mesh);
  GridFunction field(fes);
  field = [](const Geometry::Point& p) {
    Matrix value(2, 3);
    value.setConstant(p.x() + 2 * p.y());
    return value;
  };
  TrialFunction u(fes);
  TrialFunction su(scalar);
  TestFunction v(fes);
  TestFunction sv(scalar);
  auto matrixForm = BoundaryIntegral(Grad(u), Grad(v));
  auto scalarForm = BoundaryIntegral(Grad(su), Grad(sv));
  for (auto face = mesh.getBoundary(); face; ++face)
  {
    Geometry::Point point(*face, Polytope::Traits(G::Segment).getCentroid());
    const auto gradient = Grad(field).getValue(point);
    EXPECT_NEAR(gradient(1, 2, 0), 1, 1e-9);
    EXPECT_NEAR(gradient(1, 2, 1), 2, 1e-9);
    matrixForm.setPolytope(*face);
    scalarForm.setPolytope(*face);
    const size_t count = fes.getFiniteElement(1, face->getIndex()).getCount();
    for (size_t a = 0; a < count; ++a)
      for (size_t b = 0; b < count; ++b)
        EXPECT_NEAR(matrixForm.integrate(a, b),
          a % 6 == b % 6 ? scalarForm.integrate(a / 6, b / 6) : 0, 1e-9);
  }
}

TEST(MatrixRange, CurvedAndEmbeddedChainRule)
{
  for (const size_t spaceDimension : {2, 3})
  {
    auto builder = LocalMesh::Builder();
    builder.initialize(spaceDimension).nodes(3);
    if (spaceDimension == 2)
      builder.vertex({0, 0}).vertex({1, 0}).vertex({0, 1});
    else
      builder.vertex({0, 0, 0}).vertex({1, 0, 0}).vertex({0, 1, 0});
    Mesh mesh = builder.polytope(G::Triangle, {0, 1, 2}).finalize();
    RealH1Element<2> geometryElement(G::Triangle);
    Geometry::PointCloud nodes(spaceDimension, geometryElement.getCount());
    for (size_t a = 0; a < geometryElement.getCount(); ++a)
    {
      const auto& r = geometryElement.getNode(a);
      nodes(0, a) = r[0] + 0.1 * r[0] * r[1];
      nodes(1, a) = r[1] + 0.2 * r[0] * r[1];
      if (spaceDimension == 3)
        nodes(2, a) = 0.3 * r[0] * r[1];
    }
    mesh.setPolytopeTransformation({2, 0},
      new Geometry::ParametricTransformation<RealH1Element<2>>(nodes, geometryElement));
    mesh.getConnectivity().compute(2, 1);
    auto check = [&](const auto& fes) {
      GridFunction field(fes);
      field = [](const Geometry::Point& p) {
        Matrix value(2, 3);
        for (size_t i = 0; i < 2; ++i)
          for (size_t j = 0; j < 3; ++j)
            value(i, j) = (i + 1) * p.x() + (j + 1) * p.y();
        return value;
      };
      auto cell = mesh.getCell();
      auto reference = Polytope::Traits(G::Triangle).getCentroid();
      Geometry::Point point(*cell, reference);
      const auto gradient = Grad(field).getValue(point);
      const auto& inverse = point.getJacobianInverse();
      Math::SpatialTensor<Real> expected(2, 3, spaceDimension);
      expected.setZero();
      for (size_t l = 0; l < 2; ++l)
      {
        auto plus = reference, minus = reference;
        plus(l) += 1e-6;
        minus(l) -= 1e-6;
        const auto difference = (field.getValue(Geometry::Point(*cell, plus)) -
                                  field.getValue(Geometry::Point(*cell, minus))) /
          2e-6;
        for (size_t i = 0; i < 2; ++i)
          for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < spaceDimension; ++k)
              expected(i, j, k) += difference(i, j) * inverse(l, k);
      }
      EXPECT_LE((gradient - expected).norm(), 1e-7);
    };
    check(P0(mesh, 2, 3));
    check(P0g(mesh, 2, 3));
    check(P1(mesh, 2, 3));
    check(H1(std::integral_constant<size_t, 1>{}, mesh, 2, 3));
    check(H1(std::integral_constant<size_t, 2>{}, mesh, 2, 3));
    check(H1(std::integral_constant<size_t, 3>{}, mesh, 2, 3));
  }
}

TEST(MatrixRange, MixedSpaceMassUsesFrobeniusProduct)
{
  for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
         G::Pyramid, G::Wedge, G::Hexahedron})
  {
    const size_t D = Polytope::Traits(geometry).getDimension();
    Mesh mesh = D == 1 ? LocalMesh::UniformGrid(geometry, {2})
      : D == 2         ? LocalMesh::UniformGrid(geometry, {2, 2})
                       : LocalMesh::UniformGrid(geometry, {2, 2, 2});
    for (size_t d = 1; d <= D; ++d)
      for (size_t l = 0; l < d; ++l)
        mesh.getConnectivity().compute(d, l);
    auto check = [&](const auto& trialSpace, const auto& testSpace) {
      TrialFunction u(trialSpace);
      TestFunction v(testSpace);
      TrialFunction su(trialSpace.getScalarSpace());
      TestFunction sv(testSpace.getScalarSpace());
      auto matrixMass = Integral(u, v);
      auto scalarMass = Integral(su, sv);
      auto cell = mesh.getCell();
      matrixMass.setPolytope(*cell);
      scalarMass.setPolytope(*cell);
      const size_t n = trialSpace.getFiniteElement(D, cell->getIndex()).getCount();
      const size_t m = testSpace.getFiniteElement(D, cell->getIndex()).getCount();
      for (size_t a = 0; a < n; ++a)
        for (size_t b = 0; b < m; ++b)
          EXPECT_NEAR(matrixMass.integrate(a, b),
            a % 6 == b % 6 ? scalarMass.integrate(a / 6, b / 6) : 0, 1e-9);
    };
    P0 p0(mesh, 2, 3);
    P0g p0g(mesh, 2, 3);
    P1 p1(mesh, 2, 3);
    H1 h2(std::integral_constant<size_t, 2>{}, mesh, 2, 3);
    H1 h3(std::integral_constant<size_t, 3>{}, mesh, 2, 3);
    check(p0, p1);
    check(p1, p0);
    check(p0g, h2);
    check(h2, p0g);
    check(p1, h2);
    check(h2, p1);
    check(h2, h3);
    check(h3, h2);
  }
}

TEST(MatrixRange, SubMeshInclusionAndRestriction)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {3, 3});
  mesh.getConnectivity().compute(2, 1);
  for (size_t a = 0; a < mesh.getCellCount(); ++a)
    mesh.setAttribute({2, a}, a < mesh.getCellCount() / 2 ? 1 : 2);
  auto sub = mesh.keep(1);
  sub.getConnectivity().compute(2, 1);
  auto check = [&](const auto& parentSpace, const auto& childSpace) {
    GridFunction parent(parentSpace), child(childSpace);
    auto value = [](const Geometry::Point& p) {
      Matrix v(2, 3);
      if constexpr (std::is_same_v<
                      typename std::decay_t<decltype(parentSpace)>::ElementType,
                      P0gElement<Matrix>>)
        v.setConstant(1);
      else
        v.setConstant(1 + p.x() + 2 * p.y());
      return v;
    };
    parent = value;
    child = value;
    for (auto cell = sub.getCell(); cell; ++cell)
    {
      Geometry::Point point(*cell, Polytope::Traits(G::Triangle).getCentroid());
      const auto ancestor = mesh.inclusion(point);
      ASSERT_TRUE(ancestor);
      EXPECT_LE((parent.getValue(point) - child.getValue(*ancestor)).norm(), 1e-9);
      EXPECT_LE(
        (Grad(parent).getValue(point) - Grad(child).getValue(*ancestor)).norm(), 1e-9);
    }
  };
  check(P0(mesh, 2, 3), P0(sub, 2, 3));
  check(P0g(mesh, 2, 3), P0g(sub, 2, 3));
  check(P1(mesh, 2, 3), P1(sub, 2, 3));
  check(H1(std::integral_constant<size_t, 1>{}, mesh, 2, 3),
    H1(std::integral_constant<size_t, 1>{}, sub, 2, 3));
  check(H1(std::integral_constant<size_t, 2>{}, mesh, 2, 3),
    H1(std::integral_constant<size_t, 2>{}, sub, 2, 3));
  check(H1(std::integral_constant<size_t, 3>{}, mesh, 2, 3),
    H1(std::integral_constant<size_t, 3>{}, sub, 2, 3));
}

TEST(MatrixRange, InteriorDerivativeRequiresTraceDomain)
{
  Mesh mesh = LocalMesh::UniformGrid(G::Triangle, {3, 3});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  for (auto cell = mesh.getCell(); cell; ++cell)
  {
    Geometry::Point point(*cell, Polytope::Traits(G::Triangle).getCentroid());
    mesh.setAttribute(cell.key(), point.x() < 1 ? 1 : 2);
  }
  P1 fes(mesh, 2, 2);
  GridFunction field(fes);
  field = [](const Geometry::Point& p) {
    Matrix value(2, 2);
    value.setConstant(std::abs(p.x() - 1));
    return value;
  };
  size_t checked = 0;
  for (auto face = mesh.getFace(); face; ++face)
  {
    if (mesh.isBoundary(face->getIndex()))
      continue;
    Geometry::Point point(*face, Polytope::Traits(G::Segment).getCentroid());
    if (std::abs(point.x() - 1) > 1e-12)
      continue;
    EXPECT_THROW(Grad(field).getValue(point), std::exception);
    const auto left = Grad(field).traceOf(1).getValue(point);
    const auto right = Grad(field).traceOf(2).getValue(point);
    EXPECT_NEAR(left(0, 0, 0), -1, 1e-9);
    EXPECT_NEAR(right(0, 0, 0), 1, 1e-9);
    ++checked;
  }
  EXPECT_GT(checked, 0);
}

TEST(MatrixRange, ZeroDimensionalFields)
{
  Mesh mesh = LocalMesh::Builder().initialize(1).nodes(1).vertex({0.2}).finalize();
  auto check = [&](const auto& fes) {
    GridFunction field(fes);
    Matrix value(fes.getRows(), fes.getColumns());
    for (size_t i = 0; i < value.rows(); ++i)
      for (size_t j = 0; j < value.cols(); ++j)
        value(i, j) = 1 + 7 * i + j;
    field = value;
    auto vertex = mesh.getVertex();
    Geometry::Point point(
      *vertex, Polytope::Traits(G::Point).getCentroid(), vertex->getCoordinates());
    EXPECT_LE((field.getValue(point) - value).norm(), 1e-12);
    EXPECT_NEAR(Grad(field).getValue(point).norm(), 0, 1e-12);
  };
  for (auto [r, c] : matrixShapes())
  {
    check(P0(mesh, r, c));
    check(P0g(mesh, r, c));
    check(P1(mesh, r, c));
    check(H1(std::integral_constant<size_t, 1>{}, mesh, r, c));
    check(H1(std::integral_constant<size_t, 2>{}, mesh, r, c));
    check(H1(std::integral_constant<size_t, 3>{}, mesh, r, c));
  }
}

TEST(MatrixRange, StressVelocityCouplingPreservesPhysicalIndices)
{
  auto check = []<size_t K>(std::integral_constant<size_t, K> order, G geometry) {
    const size_t D = Polytope::Traits(geometry).getDimension();
    Mesh mesh = D == 2 ? LocalMesh::UniformGrid(geometry, {2, 2})
                       : LocalMesh::UniformGrid(geometry, {2, 2, 2});
    for (size_t d = 1; d <= D; ++d)
      for (size_t l = 0; l < d; ++l)
        mesh.getConnectivity().compute(d, l);
    H1 stressSpace(order, mesh, 2, D);
    H1 velocitySpace(order, mesh, 2);
    H1 scalarSpace(order, mesh);
    TrialFunction stress(stressSpace);
    TestFunction stressTest(stressSpace);
    TrialFunction velocity(velocitySpace);
    TestFunction velocityTest(velocitySpace);
    TrialFunction scalar(scalarSpace);
    TestFunction scalarTest(scalarSpace);
    SCOPED_TRACE(std::to_string(K) + ":" + std::to_string(static_cast<int>(geometry)));
    auto divergence = Integral(Div(stress), velocityTest);
    auto jacobian = Integral(Jacobian(velocity), stressTest);
    if (geometry == G::Pyramid)
    {
      divergence.setOrder(12);
      jacobian.setOrder(12);
    }
    auto cell = mesh.getCell();
    divergence.setPolytope(*cell);
    jacobian.setPolytope(*cell);
    const size_t matrixCount =
      stressSpace.getFiniteElement(D, cell->getIndex()).getCount();
    const size_t vectorCount =
      velocitySpace.getFiniteElement(D, cell->getIndex()).getCount();
    for (size_t direction = 0; direction < D; ++direction)
    {
      auto derivative = Integral(Component(Grad(scalar), direction), scalarTest);
      if (geometry == G::Pyramid)
        derivative.setOrder(12);
      derivative.setPolytope(*cell);
      for (size_t a = 0; a < matrixCount; ++a)
        for (size_t b = 0; b < vectorCount; ++b)
        {
          const size_t row = (a % (2 * D)) / D, column = a % D;
          if (column == direction)
            EXPECT_NEAR(divergence.integrate(a, b),
              row == b % 2 ? derivative.integrate(a / (2 * D), b / 2) : 0, 1e-9);
          const size_t testRow = (a % (2 * D)) / D;
          if (column == direction)
            EXPECT_NEAR(jacobian.integrate(b, a),
              testRow == b % 2 ? derivative.integrate(b / 2, a / (2 * D)) : 0, 1e-9);
        }
    }
    BilinearForm assembled(stress, velocityTest);
    assembled = Integral(Div(stress), velocityTest);
    assembled.assemble();
    EXPECT_EQ(assembled.getOperator().rows(), velocitySpace.getSize());
    EXPECT_EQ(assembled.getOperator().cols(), stressSpace.getSize());
  };
  for (auto geometry :
    {G::Triangle, G::Quadrilateral, G::Tetrahedron, G::Pyramid, G::Wedge, G::Hexahedron})
  {
    check(std::integral_constant<size_t, 1>{}, geometry);
    check(std::integral_constant<size_t, 2>{}, geometry);
    check(std::integral_constant<size_t, 3>{}, geometry);
  }
}
