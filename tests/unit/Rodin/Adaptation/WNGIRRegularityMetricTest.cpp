/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include "Rodin/Adaptation/WNGIR/Distribution.h"
#include "Rodin/Geometry.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Adaptation_WNGIRRegularityMetric, DeviatoricCurrentStrainKernelAndEnergyP1P2P3)
{
  const auto check = []<size_t Order, size_t Dimension>() {
    auto mesh = [&] {
      if constexpr (Dimension == 2)
        return LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
      else
        return LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    }();
    for (size_t from = 1; from <= Dimension; ++from)
      for (size_t to = 0; to <= Dimension; ++to)
        if (from != to)
          mesh.getConnectivity().compute(from, to);
    auto fes = [&] {
      if constexpr (Order == 1)
        return P1<Math::SpatialVector<Real>, LocalMesh>(mesh, Dimension);
      else
        return H1(std::integral_constant<size_t, Order>{}, mesh, Dimension);
    }();
    TrialFunction trial(fes);
    TestFunction test(fes);
    GridFunction current(fes), position(fes), field(fes);
    current = VectorFunction(Dimension, [](const Point& p) {
      Math::SpatialVector<Real> value(Dimension);
      value.setZero();
      value(0) = Real(0.2) * p.getCoordinates()(0) * p.getCoordinates()(0);
      value(1) = Real(-0.1) * p.getCoordinates()(1);
      return value;
    });
    position = VectorFunction(Dimension,
      [](const Point& p) { return Math::SpatialVector<Real>(p.getCoordinates()); });
    position += current;
    constexpr Real coefficient = Real(0.7);
    const auto weight = Adaptation::wngirCurrentVolumeWeight(current, Dimension);
    auto integral = Integral(
      coefficient * weight * Adaptation::wngirCurrentStrain(trial, current, Dimension),
      Adaptation::wngirCurrentStrain(test, current, Dimension));
    integral.setOrder(2 * Order);
    auto traceIntegral = Integral((coefficient / Real(Dimension)) * weight *
        Trace(Adaptation::wngirCurrentStrain(trial, current, Dimension)),
      Trace(Adaptation::wngirCurrentStrain(test, current, Dimension)));
    traceIntegral.setOrder(2 * Order);
    BilinearForm form(trial, test);
    form = integral - traceIntegral;
    form.assemble();
    BilinearForm tabulated(trial, test);
    tabulated =
      Adaptation::WNGIRDistribution(trial, test, current, coefficient, 2 * Order);
    tabulated.assemble();
    EXPECT_LT((tabulated.getOperator() - form.getOperator()).norm(),
      Real(1e-12) * form.getOperator().norm());
    Math::Matrix<Real> matrix(form.getOperator());
    EXPECT_LT((matrix - matrix.transpose()).norm(), Real(1e-12) * matrix.norm());
    Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eigen(matrix);
    EXPECT_GE(eigen.eigenvalues().minCoeff(), Real(-1e-11) * matrix.norm());
    field = VectorFunction(Dimension, [](const Point&) {
      Math::SpatialVector<Real> value(Dimension);
      value.setConstant(Real(1));
      return value;
    });
    EXPECT_LT((matrix * field.getData()).norm(),
      Real(1e-10) * matrix.norm() * field.getData().norm());
    field = position;
    EXPECT_LT((matrix * field.getData()).norm(),
      Real(1e-10) * matrix.norm() * field.getData().norm());
    Math::Matrix<Real> rotation = Math::Matrix<Real>::Zero(Dimension, Dimension);
    rotation(0, 1) = Real(-1);
    rotation(1, 0) = Real(1);
    field = MatrixFunction(rotation) * position;
    EXPECT_LT((matrix * field.getData()).norm(),
      Real(1e-10) * matrix.norm() * field.getData().norm());
    field = VectorFunction(Dimension, [](const Point& p) {
      Math::SpatialVector<Real> value(Dimension);
      for (size_t i = 0; i < Dimension; ++i)
        value(i) = Real(i + 1) * p.getCoordinates()(i) * p.getCoordinates()(i);
      return value;
    });
    auto jc = Jacobian(current), jf = Jacobian(field);
    Real norm = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const auto& qf = QF::PolytopeQuadratureFormula::get(2 * Order, cell->getGeometry());
      const auto& quadrature = cell->getQuadrature(qf);
      for (size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& point = quadrature.getPoint(q);
        const IntegrationPoint ip(point, &qf, q);
        Adaptation::CellDeformation d(Dimension);
        d.setDisplacementGradient(jc.getValue(ip));
        const Math::SpatialMatrix<Real> L(
          jf.getValue(ip) * d.getDeformationGradient().inverse());
        const Math::SpatialMatrix<Real> strain(Real(0.5) * (L + L.transpose()));
        const Real w = qf.getWeight(q) * point.getDistortion() * d.getJacobian();
        norm +=
          w * (strain.squaredNorm() - strain.trace() * strain.trace() / Real(Dimension));
      }
    }
    const Real expected = coefficient * norm;
    const Real actual = field.getData().dot(matrix * field.getData());
    EXPECT_GT(actual, Real(0));
    EXPECT_NEAR(actual, expected, Real(1e-10) * std::max(Real(1), std::abs(expected)));
  };
  check.template operator()<1, 2>();
  check.template operator()<2, 2>();
  check.template operator()<1, 3>();
  check.template operator()<2, 3>();
  check.template operator()<3, 2>();
  check.template operator()<3, 3>();
}

TEST(Rodin_Adaptation_WNGIRRegularityMetric, StableQuadratureAndCompleteConformalKernel)
{
  const auto check = []<size_t Order, size_t Dimension>(bool curved) {
    auto mesh = [&] {
      if constexpr (Dimension == 2)
        return LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
      else
        return LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    }();
    for (size_t from = 1; from <= Dimension; ++from)
      for (size_t to = 0; to <= Dimension; ++to)
        if (from != to)
          mesh.getConnectivity().compute(from, to);
    auto fes = [&] {
      if constexpr (Order == 1)
        return P1<Math::SpatialVector<Real>, LocalMesh>(mesh, Dimension);
      else
        return H1(std::integral_constant<size_t, Order>{}, mesh, Dimension);
    }();
    TrialFunction trial(fes);
    TestFunction test(fes);
    GridFunction current(fes), position(fes), mode(fes);
    current = VectorFunction(Dimension, [=](const Point& p) {
      Math::SpatialVector<Real> value(Dimension);
      const auto& x = p.getCoordinates();
      value.setZero();
      value(0) = Real(0.2) * (curved ? x(0) * x(0) : x(0));
      value(1) = Real(-0.1) * x(1);
      return value;
    });
    position = VectorFunction(Dimension,
      [](const Point& p) { return Math::SpatialVector<Real>(p.getCoordinates()); });
    position += current;
    const auto assemble = [&](size_t order) {
      BilinearForm form(trial, test);
      form = Adaptation::WNGIRDistribution(trial, test, current, Real(1), order);
      form.assemble();
      Math::Matrix<Real> matrix(form.getOperator());
      return Math::Matrix<Real>((Real(0.5) * (matrix + matrix.transpose())).eval());
    };
    const auto actual = assemble(2 * Order), reference = assemble(12);
    if (!curved || Order == 1)
      EXPECT_LT((actual - reference).norm(), Real(1e-11) * reference.norm());
    const size_t count = Dimension * (Dimension + 1) / 2 + 1;
    Math::Matrix<Real> modes(actual.rows(), count);
    size_t column = 0;
    for (size_t axis = 0; axis < Dimension; ++axis)
    {
      mode = VectorFunction(Dimension, [=](const Point&) {
        Math::SpatialVector<Real> value = Math::SpatialVector<Real>::Zero(Dimension);
        value(axis) = Real(1);
        return value;
      });
      modes.col(column++) = mode.getData();
    }
    for (size_t a = 0; a < Dimension; ++a)
      for (size_t b = a + 1; b < Dimension; ++b)
      {
        Math::Matrix<Real> rotation = Math::Matrix<Real>::Zero(Dimension, Dimension);
        rotation(a, b) = Real(-1);
        rotation(b, a) = Real(1);
        mode = MatrixFunction(rotation) * position;
        modes.col(column++) = mode.getData();
      }
    modes.col(column) = position.getData();
    Eigen::HouseholderQR<Math::Matrix<Real>> qr(modes);
    const Math::Matrix<Real> orthogonal =
      qr.householderQ() * Math::Matrix<Real>::Identity(actual.rows(), actual.rows());
    const Math::Matrix<Real> kernel = orthogonal.leftCols(count);
    EXPECT_LT((actual * kernel).norm(), Real(1e-10) * actual.norm());
    Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> spectrum(actual);
    ASSERT_EQ(spectrum.info(), Eigen::Success);
    const size_t kernelCount =
      Order == 1 ? count : (Dimension == 2 ? 2 * (Order + 1) : 10);
    const Real tolerance = Real(1e-10) * actual.norm();
    EXPECT_EQ((spectrum.eigenvalues().array().abs() < tolerance).count(), kernelCount);
    EXPECT_GE(spectrum.eigenvalues().minCoeff(), -tolerance);
    EXPECT_GT(spectrum.eigenvalues()(kernelCount), tolerance);
    const Math::Matrix<Real> complement =
      spectrum.eigenvectors().rightCols(actual.rows() - kernelCount);
    const Math::Matrix<Real> reduced = complement.transpose() * actual * complement;
    const Math::Matrix<Real> referenceReduced =
      complement.transpose() * reference * complement;
    Eigen::GeneralizedSelfAdjointEigenSolver<Math::Matrix<Real>> comparison(
      reduced, referenceReduced);
    ASSERT_EQ(comparison.info(), Eigen::Success);
    EXPECT_GT(comparison.eigenvalues().minCoeff(), Real(0.95));
    EXPECT_LT(comparison.eigenvalues().maxCoeff(), Real(1.05));
    std::cout << "quadrature P" << Order << " d=" << Dimension << " curved=" << curved
              << " dofs=" << actual.rows() << " kernel=" << kernelCount << " ratio=["
              << comparison.eigenvalues().minCoeff() << ","
              << comparison.eigenvalues().maxCoeff() << "]\n";
    const auto geometry =
      Dimension == 2 ? Polytope::Type::Triangle : Polytope::Type::Tetrahedron;
    const auto& rule = QF::PolytopeQuadratureFormula::get(2 * Order, geometry);
    for (size_t q = 0; q < rule.getSize(); ++q)
      EXPECT_GT(rule.getWeight(q), Real(0));
    if constexpr (Order == 2)
    {
      // Positive centroid weights alone do not prevent hourglass modes.
      const auto underintegrated = assemble(1);
      Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eigen(underintegrated);
      EXPECT_GT((eigen.eigenvalues().array().abs() < Real(1e-10) * underintegrated.norm())
                  .count(),
        kernelCount);
    }
  };
  for (const bool curved : {false, true})
  {
    check.template operator()<1, 2>(curved);
    check.template operator()<1, 3>(curved);
    if (!curved)
    {
      check.template operator()<2, 2>(curved);
      check.template operator()<2, 3>(curved);
      check.template operator()<3, 2>(curved);
      check.template operator()<3, 3>(curved);
    }
  }
}
