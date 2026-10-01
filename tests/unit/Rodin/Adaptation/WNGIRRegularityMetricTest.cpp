/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include "Rodin/Adaptation/WNGIRRegularityMetric.h"
#include "Rodin/Geometry.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Adaptation_WNGIRRegularityMetric, CurrentStrainKernelAndEnergyP1P2)
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
    const auto weight = Adaptation::Detail::wngirCurrentVolumeWeight(current, Dimension);
    auto integral = Integral(coefficient * weight *
        Adaptation::Detail::wngirCurrentStrain(trial, current, Dimension),
      Adaptation::Detail::wngirCurrentStrain(test, current, Dimension));
    integral.setOrder(2 * Order);
    BilinearForm form(trial, test);
    form = integral;
    form.assemble();
    const auto modes = Adaptation::Detail::wngirCurrentStrainCouplings(
      test, current, Dimension, coefficient, 2 * Order);
    ASSERT_EQ(modes.size(), 1u);
    Math::Matrix<Real> matrix(form.getOperator());
    matrix -= modes[0] * modes[0].transpose();
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
    Real norm = 0, meanTrace = 0, measure = 0;
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
        norm += w * strain.squaredNorm();
        meanTrace += w * strain.trace();
        measure += w;
      }
    }
    const Real expected =
      coefficient * (norm - meanTrace * meanTrace / (Real(Dimension) * measure));
    const Real actual = field.getData().dot(matrix * field.getData());
    EXPECT_GT(actual, Real(0));
    EXPECT_NEAR(actual, expected, Real(1e-10) * std::max(Real(1), std::abs(expected)));
    const auto empty = Adaptation::Detail::wngirCurrentStrainCouplings(
      test, current, Dimension, Real(0), 2 * Order);
    EXPECT_TRUE(empty.empty());
    EXPECT_THROW(Adaptation::Detail::wngirCurrentStrainCouplings(
                   test, current, Dimension, Real(-1), 2 * Order),
      Alert::Exception);
  };
  check.template operator()<1, 2>();
  check.template operator()<2, 2>();
  check.template operator()<1, 3>();
  check.template operator()<2, 3>();
}

TEST(Rodin_Adaptation_WNGIRRegularityMetric, SignedInertiaMatchesDense)
{
  Math::Matrix<Real> dense = Math::Matrix<Real>::Identity(5, 5);
  std::vector<Math::Vector<Real>> modes;
  auto mode = Math::Vector<Real>::Zero(5).eval();
  mode(0) = Real(2);
  modes.push_back(mode);
  Math::SparseMatrix<Real> sparse = dense.sparseView();
  EXPECT_EQ(Adaptation::Detail::wngirRegularityInertia(sparse, modes, {Real(-1)}), 1);
  EXPECT_EQ(Adaptation::Detail::wngirRegularityInertia(sparse, modes, {Real(1)}), 0);
  mode(1) = Real(1);
  modes.push_back(mode);
  const std::vector<Real> weights{Real(-1), Real(2)};
  dense -= modes[0] * modes[0].transpose();
  dense += Real(2) * modes[1] * modes[1].transpose();
  Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eigen(dense);
  EXPECT_EQ(Adaptation::Detail::wngirRegularityInertia(sparse, modes, weights),
    (eigen.eigenvalues().array() < Real(0)).count());
}
