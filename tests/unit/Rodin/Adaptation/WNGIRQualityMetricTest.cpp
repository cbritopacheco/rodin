/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <Eigen/Eigenvalues>
#include "Rodin/Adaptation/WNGIRQualityMetric.h"
#include "Rodin/Assembly.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Adaptation;

TEST(Rodin_Adaptation_WNGIRQualityMetric, FullShapeCurvatureP1P2In2D3D)
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
    GridFunction current(fes), direction(fes);
    current = VectorFunction(Dimension, [](const Point& p) {
      Math::SpatialVector<Real> v(Dimension);
      v.setZero();
      v(0) = Real(0.4) * p.getCoordinates()(0);
      v(1) = Real(-0.3) * p.getCoordinates()(1);
      return v;
    });
    direction = VectorFunction(Dimension, [](const Point& p) {
      Math::SpatialVector<Real> v(Dimension);
      v = Real(0.2) * p.getCoordinates();
      v(0) += Real(0.1) * p.getCoordinates()(1) * p.getCoordinates()(1);
      return v;
    });
    TrialFunction trial(fes);
    TestFunction test(fes);
    WNGIRParameters parameters;
    parameters.h = Real(0.5);
    parameters.kappaBulk = Real(0.7);
    parameters.kappaC = Real(1.3);
    BilinearForm metric(trial, test);
    metric = Detail::WNGIRQualityMetric(trial, test, current, parameters);
    metric.assemble();
    const auto gradient = Jacobian(current), increment = Jacobian(direction);
    const auto energy = [&](Real t) {
      Real result = 0;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(2 * Order, cell->getGeometry());
        const auto& quadrature = cell->getQuadrature(qf);
        for (size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& point = quadrature.getPoint(q);
          const IntegrationPoint ip(point, &qf, q);
          CellDeformation deformation(Dimension);
          deformation.setDisplacementGradient(Math::SpatialMatrix<Real>(
            gradient.getValue(ip) + t * increment.getValue(ip)));
          result += qf.getWeight(q) * point.getDistortion() * Real(Dimension) / Real(4) *
            (deformation.getRelativeDistortion() - Real(1));
        }
      }
      return parameters.h * parameters.kappaBulk * parameters.kappaC * result;
    };
    const Real actual =
      direction.getData().dot(metric.getOperator() * direction.getData());
    constexpr Real eps = Real(1e-4);
    const Real expected =
      (energy(eps) - Real(2) * energy(0) + energy(-eps)) / (eps * eps);
    EXPECT_NEAR(actual, expected, Real(2e-6) * std::max(Real(1), std::abs(expected)));
    const Math::Matrix<Real> dense(metric.getOperator());
    EXPECT_LT((dense - dense.transpose()).norm(), Real(1e-12));
    Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eigen(dense);
    ASSERT_EQ(eigen.info(), Eigen::Success);
    EXPECT_LT(eigen.eigenvalues().minCoeff(), Real(-1e-7));
    parameters.kappaC *= Real(2);
    metric = Detail::WNGIRQualityMetric(trial, test, current, parameters);
    metric.assemble();
    EXPECT_NEAR(direction.getData().dot(metric.getOperator() * direction.getData()),
      Real(2) * actual, Real(1e-12));
  };
  check.template operator()<1, 2>();
  check.template operator()<2, 2>();
  check.template operator()<1, 3>();
  check.template operator()<2, 3>();
}

TEST(Rodin_Adaptation_WNGIRQualityMetric, IdentityIsDeviatoricStrainWithoutVolumePenalty)
{
  for (size_t dimension : {2u, 3u})
  {
    CellDeformation deformation(dimension);
    for (size_t mode = 0; mode < 3; ++mode)
    {
      Math::SpatialMatrix<Real> G =
        Math::SpatialMatrix<Real>::Identity(dimension, dimension);
      if (mode == 1)
      {
        G.setZero();
        G(0, 1) = Real(1);
        G(1, 0) = Real(-1);
      }
      else if (mode == 2)
      {
        G(0, 0) = Real(0.3);
        G(1, 1) = Real(-0.2);
        G(0, 1) = Real(0.7);
      }
      Math::SpatialMatrix<Real> strain(Real(0.5) * (G + G.transpose()));
      const Real mean = strain.trace() / Real(dimension);
      for (size_t i = 0; i < dimension; ++i)
        strain(i, i) -= mean;
      const Real actual =
        Real(dimension) / Real(4) * deformation.getRelativeDistortionSecondAction(G, G);
      EXPECT_NEAR(actual, strain.squaredNorm(), Real(1e-12));
      if (mode < 2)
        EXPECT_NEAR(actual, Real(0), Real(1e-12));
    }
  }
}
