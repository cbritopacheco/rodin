/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "../../../QuadratureReference.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief The tangential field (0,y,z) on F(r,s)=(r,s,s) has divergence one.
  template <size_t Degree>
  void checkEmbeddedDivergence()
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(3)
                  .nodes(3)
                  .vertex({0, 0, 0})
                  .vertex({1, 0, 0})
                  .vertex({0, 1, 1})
                  .polytope(Polytope::Type::Triangle, {0, 1, 2})
                  .finalize();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    auto vectors = [&] {
      if constexpr (Degree == 0)
        return P1<Math::SpatialVector<Real>>(mesh, 3);
      else
      {
        return H1<Degree, Math::SpatialVector<Real>>(
          std::integral_constant<size_t, Degree>{}, mesh, 3);
      }
    }();
    auto scalars = [&] {
      if constexpr (Degree == 0)
        return P1<Real>(mesh);
      else
        return H1<Degree, Real>(std::integral_constant<size_t, Degree>{}, mesh);
    }();
    GridFunction field(vectors);
    field = VectorFunction(size_t(3), [](const Point& point) {
      return Math::SpatialVector<Real>{0, point.y(), point.z()};
    });
    TrialFunction u(vectors);
    TrialFunction w(scalars);
    TestFunction v(vectors);
    TestFunction q(scalars);
    const auto cell = mesh.getCell();
    const auto& formula = QF::QuadratureFormula::get(5, cell->getGeometry());
    const auto& quadrature = cell->getQuadrature(formula);
    auto derivative = Div(field);
    auto basis = Div(u);
    basis.setPolytope(*cell);
    for (size_t sample = 0; sample < quadrature.getSize(); ++sample)
    {
      const auto& point = quadrature.getPoint(sample);
      const IntegrationPoint ip(quadrature, sample);
      EXPECT_NEAR(derivative.getValue(point), 1, 1e-12);
      EXPECT_NEAR(derivative.getValue(ip), 1, 1e-12);
      basis.setIntegrationPoint(ip);
      Real reconstructed = 0;
      for (size_t local = 0; local < u.getDOFs(*cell); ++local)
      {
        reconstructed +=
          field.getData()[vectors.getGlobalIndex({2, cell->getIndex()}, local)] *
          basis.getBasis(local);
      }
      EXPECT_NEAR(reconstructed, 1, 1e-12);
    }
    const auto check = [&](auto integral) {
      integral.setOrder(5);
      Rodin::Tests::QuadratureReference<decltype(integral)> reference(integral);
      reference.assemble(*cell, 5);
      integral.setPolytope(*cell);
      const auto& expected = reference.getOperator();
      for (Eigen::Index row = 0; row < expected.rows(); ++row)
      {
        for (Eigen::Index col = 0; col < expected.cols(); ++col)
          EXPECT_NEAR(integral.integrate(col, row), expected(row, col), 1e-12);
      }
    };
    check(Integral(Div(u), q));
    check(Integral(w, Div(v)));
  }
}

/// @brief Native P1 fields, shapes and both mixed-integral orientations use ambient axes.
TEST(EmbeddedDiv, P1)
{
  checkEmbeddedDivergence<0>();
}

/// @brief Linear H1 fields, shapes and both mixed-integral orientations use ambient axes.
TEST(EmbeddedDiv, H1Linear)
{
  checkEmbeddedDivergence<1>();
}

/// @brief Quadratic H1 evaluation retains the ambient divergence at every sample.
TEST(EmbeddedDiv, H1Quadratic)
{
  checkEmbeddedDivergence<2>();
}
