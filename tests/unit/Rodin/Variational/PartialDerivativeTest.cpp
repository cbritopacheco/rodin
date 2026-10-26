/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief Scalar partial derivatives use the physical gradient binding and order.
TEST(PartialDerivative, ScalarShapeBinding)
{
  auto mesh = LocalMesh::Builder()
                .initialize(2)
                .nodes(3)
                .vertex({0, 0})
                .vertex({1, 0})
                .vertex({0, 1})
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .finalize();
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  const auto cell = mesh.getCell();
  const Point point(*cell, Math::SpatialPoint{0.2, 0.3});
  const IntegrationPoint ip(point);
  const auto check = [&]<class Space>(const Space& space) {
    TrialFunction u(space);
    TestFunction v(space);
    auto dx = Dx(u);
    auto dy = Dy(u);
    auto gradient = Grad(u);
    dx.setPolytope(*cell).setIntegrationPoint(ip);
    dy.setPolytope(*cell).setIntegrationPoint(ip);
    gradient.setPolytope(*cell).setIntegrationPoint(ip);
    EXPECT_EQ(dx.getOrder(*cell), gradient.getOrder(*cell));
    auto copied = dx;
    auto moved = std::move(dy);
    copied.setPolytope(*cell).setIntegrationPoint(ip);
    moved.setPolytope(*cell).setIntegrationPoint(ip);
    for (size_t local = 0; local < u.getDOFs(*cell); ++local)
    {
      EXPECT_EQ(copied.getBasis(local), gradient.getBasis(local)[0]);
      EXPECT_EQ(moved.getBasis(local), gradient.getBasis(local)[1]);
    }
    auto invalid = Derivative(size_t(2), u);
    invalid.setPolytope(*cell).setIntegrationPoint(ip);
    EXPECT_THROW(invalid.getBasis(0), Alert::Exception);
    auto integral = Integral(Dx(u), v);
    integral.setPolytope(*cell);
    EXPECT_TRUE(std::isfinite(std::abs(integral.integrate(0, 0))));
  };
  P1<Real> real(mesh);
  P1<Complex> complex(mesh);
  H1<2, Real> quadratic(std::integral_constant<size_t, 2>{}, mesh);
  check(real);
  check(complex);
  check(quadratic);
}
