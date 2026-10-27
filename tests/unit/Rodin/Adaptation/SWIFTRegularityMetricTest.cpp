/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <Eigen/Eigenvalues>
#include "Rodin/Adaptation/SWIFT/Distribution.h"
#include "Rodin/Geometry.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Adaptation_SWIFTRegularityMetric, CenteredVarianceKernelAndNativeFormsP1P2P3)
{
  const auto check = []<size_t Order, size_t Dimension>(bool deformed) {
    auto mesh = [&] {
      if constexpr (Dimension == 1)
        return LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
      else if constexpr (Dimension == 2)
        return LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
      else
        return LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    }();
    for (size_t from = 1; from <= Dimension; ++from)
    {
      for (size_t to = 0; to <= Dimension; ++to)
      {
        if (from != to)
          mesh.getConnectivity().compute(from, to);
      }
    }
    auto fes = [&] {
      if constexpr (Order == 1)
        return P1<Math::SpatialVector<Real>, LocalMesh>(mesh, Dimension);
      else
        return H1(std::integral_constant<size_t, Order>{}, mesh, Dimension);
    }();
    TrialFunction trial(fes);
    TestFunction test(fes);
    GridFunction current(fes), position(fes), field(fes);
    current = VectorFunction(Dimension, [=](const Point& point) {
      Math::SpatialVector<Real> value = Math::SpatialVector<Real>::Zero(Dimension);
      for (size_t axis = 0; axis < Dimension; ++axis)
      {
        value(axis) = deformed
          ? Real(0.1) * point.getCoordinates()(axis) * point.getCoordinates()(axis)
          : Real(0);
      }
      return value;
    });
    position = VectorFunction(Dimension, [](const Point& point) {
      return Math::SpatialVector<Real>(point.getCoordinates());
    });
    position += current;
    const size_t order = 2 * Order;
    const auto weight =
      RealFunction([gradient = Jacobian(current)](const auto& point) -> Real {
        Adaptation::CellDeformation deformation(Dimension);
        deformation.setDisplacementGradient(gradient.getValue(point));
        return deformation.getJacobian();
      });
    const auto trialGradient =
      Jacobian(trial) * Adaptation::SWIFT::CurrentInverse(current, Dimension);
    const auto testGradient =
      Jacobian(test) * Adaptation::SWIFT::CurrentInverse(current, Dimension);
    const auto trialStrain = Real(0.5) * (trialGradient + Transpose(trialGradient));
    const auto testStrain = Real(0.5) * (testGradient + Transpose(testGradient));
    for (const auto& [dev, div] : {std::pair<Real, Real>{0, 0}, {0, 1}, {1, 0},
           {Real(0.1), 10}, {10, Real(0.1)}, {1, 1}})
    {
      SCOPED_TRACE(testing::Message()
        << "P" << Order << " d=" << Dimension << " deformed=" << deformed
        << " dev=" << dev << " div=" << div);
      BilinearForm core(trial, test), native(trial, test);
      const Adaptation::SWIFT::Distribution distribution(
        trial, test, current, dev, div, order);
      core = distribution;
      core.assemble();
      auto strainIntegral = Integral(dev * weight * trialStrain, testStrain);
      strainIntegral.setOrder(order);
      auto traceIntegral = Integral(
        (div - dev) / Real(Dimension) * weight * Trace(trialStrain), Trace(testStrain));
      traceIntegral.setOrder(order);
      native = strainIntegral + traceIntegral;
      native.assemble();
      EXPECT_LT((core.getOperator() - native.getOperator()).norm(),
        Real(1e-11) * std::max(Real(1), native.getOperator().norm()));
      const auto couplings = distribution.getCentering();
      const Math::Matrix<Real> matrix =
        Math::Matrix<Real>(core.getOperator()) - couplings * couplings.transpose();
      const Real tolerance = Real(1e-10) * std::max(Real(1), matrix.norm());
      EXPECT_LT((matrix - matrix.transpose()).norm(), tolerance);
      Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> spectrum(matrix);
      ASSERT_EQ(spectrum.info(), Eigen::Success);
      EXPECT_GE(spectrum.eigenvalues().minCoeff(), -tolerance);
      if (dev > 0 && div > 0)
      {
        EXPECT_EQ((spectrum.eigenvalues().array().abs() < tolerance).count(),
          Dimension * (Dimension + 1));
      }

      // Every global affine current-coordinate motion has constant strain.
      for (size_t row = 0; row < Dimension; ++row)
      {
        field = VectorFunction(Dimension, [=](const Point&) {
          Math::SpatialVector<Real> value = Math::SpatialVector<Real>::Zero(Dimension);
          value(row) = 1;
          return value;
        });
        EXPECT_LT((matrix * field.getData()).norm(), tolerance * field.getData().norm());
        for (size_t column = 0; column < Dimension; ++column)
        {
          Math::Matrix<Real> affine = Math::Matrix<Real>::Zero(Dimension, Dimension);
          affine(row, column) = 1;
          field = MatrixFunction(affine) * position;
          EXPECT_LT(
            (matrix * field.getData()).norm(), tolerance * field.getData().norm());
        }
      }
      field = VectorFunction(Dimension, [](const Point& point) {
        Math::SpatialVector<Real> value(Dimension);
        for (size_t axis = 0; axis < Dimension; ++axis)
        {
          value(axis) =
            Real(axis + 1) * point.getCoordinates()(axis) * point.getCoordinates()(axis);
        }
        return value;
      });
      auto currentJacobian = Jacobian(current);
      auto fieldJacobian = Jacobian(field);
      Real volume = 0, secondDev = 0, firstDiv = 0, secondDiv = 0;
      Math::SpatialMatrix<Real> firstDev(Dimension, Dimension);
      firstDev.setZero();
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto& rule = QF::QuadratureFormula::get(order, cell->getGeometry());
        const auto& quadrature = cell->getQuadrature(rule);
        for (size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& point = quadrature.getPoint(q);
          const IntegrationPoint ip(point, &rule, q);
          Adaptation::CellDeformation deformation(Dimension);
          deformation.setDisplacementGradient(currentJacobian.getValue(ip));
          const Math::SpatialMatrix<Real> gradient(
            fieldJacobian.getValue(ip) * deformation.getDeformationGradient().inverse());
          Math::SpatialMatrix<Real> strain(Real(0.5) * (gradient + gradient.transpose()));
          const Real divergence = strain.trace();
          for (size_t axis = 0; axis < Dimension; ++axis)
            strain(axis, axis) -= divergence / Real(Dimension);
          const Real w =
            rule.getWeight(q) * point.getDistortion() * deformation.getJacobian();
          volume += w;
          firstDev += w * strain;
          secondDev += w * strain.squaredNorm();
          firstDiv += w * divergence;
          secondDiv += w * divergence * divergence;
        }
      }
      const Real expected = dev * (secondDev - firstDev.squaredNorm() / volume) +
        div / Real(Dimension) * (secondDiv - firstDiv * firstDiv / volume);
      EXPECT_NEAR(field.getData().dot(matrix * field.getData()), expected, tolerance);
    }
  };
  for (const bool deformed : {false, true})
  {
    check.template operator()<1, 1>(deformed);
    check.template operator()<2, 1>(deformed);
    check.template operator()<3, 1>(deformed);
    check.template operator()<1, 2>(deformed);
    check.template operator()<2, 2>(deformed);
    check.template operator()<3, 2>(deformed);
    check.template operator()<1, 3>(deformed);
    check.template operator()<2, 3>(deformed);
    check.template operator()<3, 3>(deformed);
  }
}

TEST(Rodin_Adaptation_SWIFTRegularityMetric,
  BackgroundScaleMultipliesTheCompleteCenteredForm)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  for (size_t from = 1; from <= 2; ++from)
  {
    for (size_t to = 0; to <= 2; ++to)
    {
      if (from != to)
        mesh.getConnectivity().compute(from, to);
    }
  }
  P1<Math::SpatialVector<Real>, LocalMesh> fes(mesh, 2);
  TrialFunction trial(fes);
  TestFunction test(fes);
  GridFunction current(fes);
  current = VectorFunction(Real(0), Real(0));
  Math::Matrix<Real> reference;
  for (const Real h : {Real(1), Real(0.5), Real(0.25)})
  {
    BilinearForm core(trial, test);
    const Adaptation::SWIFT::Distribution distribution(trial, test, current, h, 2 * h, 2);
    core = distribution;
    core.assemble();
    const auto couplings = distribution.getCentering();
    const Math::Matrix<Real> matrix =
      Math::Matrix<Real>(core.getOperator()) - couplings * couplings.transpose();
    if (h == 1)
      reference = matrix;
    EXPECT_LT((matrix - h * reference).norm(), Real(1e-12) * reference.norm());
  }
}

/// @brief Curved reference geometry does not change the zero energy of constant current strain.
TEST(Rodin_Adaptation_SWIFTRegularityMetric, CurvedGeometryWithLowerOrderDisplacement)
{
  const auto check = []<size_t Degree>(
                       Polytope::Type geometry, Real curvature, Real motion) {
    auto mesh = [&] {
      if (geometry == Polytope::Type::Triangle)
      {
        return LocalMesh::Builder()
          .initialize(2)
          .nodes(3)
          .vertex({0, 0})
          .vertex({1 + curvature, 0})
          .vertex({0, 1})
          .polytope(geometry, {0, 1, 2})
          .finalize();
      }
      return LocalMesh::Builder()
        .initialize(2)
        .nodes(4)
        .vertex({0, 0})
        .vertex({1 + curvature, 0})
        .vertex({1 + curvature, 1})
        .vertex({0, 1})
        .polytope(geometry, {0, 1, 2, 3})
        .finalize();
    }();
    for (size_t from = 1; from <= 2; ++from)
    {
      for (size_t to = 0; to <= 2; ++to)
      {
        if (from != to)
          mesh.getConnectivity().compute(from, to);
      }
    }
    RealH1Element<2> element(geometry);
    PointCloud nodes(2, element.getCount());
    for (size_t local = 0; local < element.getCount(); ++local)
    {
      const auto r = element.getNode(local);
      nodes(0, local) = r[0] + curvature * r[0] * r[0];
      nodes(1, local) = r[1];
    }
    mesh.setPolytopeTransformation({2, 0},
      new ParametricTransformation<RealH1Element<2>>(
        std::move(nodes), std::move(element)));
    auto space = [&] {
      if constexpr (Degree == 1)
        return P1<Math::SpatialVector<Real>>(mesh, 2);
      else
      {
        return H1<Degree, Math::SpatialVector<Real>>(
          std::integral_constant<size_t, Degree>{}, mesh, 2);
      }
    }();
    TrialFunction trial(space);
    TestFunction test(space);
    GridFunction current(space), field(space);
    current = VectorFunction(size_t(2),
      [motion](const Point& p) { return Math::SpatialVector<Real>{0, motion * p.y()}; });
    field = VectorFunction(
      size_t(2), [](const Point& p) { return Math::SpatialVector<Real>{0, p.y()}; });
    const Adaptation::SWIFT::Distribution distribution(trial, test, current, 1, 1, 12);
    BilinearForm core(trial, test);
    core = distribution;
    core.assemble();
    const auto coupling = distribution.getCentering();
    Math::Matrix<Real> matrix = core.getOperator();
    const auto coefficients = field.getData().eval();
    const Real volume =
      geometry == Polytope::Type::Triangle ? Real(0.5) + curvature / 3 : 1 + curvature;
    EXPECT_NEAR(coefficients.dot(matrix * coefficients), volume / (1 + motion), 3e-12);
    matrix -= coupling * coupling.transpose();
    EXPECT_NEAR(coefficients.dot(matrix * coefficients), 0, 5e-12);
    const Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eigenvalues(matrix);
    EXPECT_GT(eigenvalues.eigenvalues().minCoeff(), -1e-11);
  };
  for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral})
  {
    for (const Real curvature : {Real(-0.4), Real(1)})
    {
      for (const Real motion : {Real(0), Real(0.2)})
      {
        check.template operator()<1>(geometry, curvature, motion);
        check.template operator()<2>(geometry, curvature, motion);
      }
    }
  }
}
