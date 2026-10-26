/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/// @file @brief Checked derivative dispatch on mesh ancestry and excluded cells.
#include <gtest/gtest.h>
#include <array>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <size_t Degree>
  void checkDerivativeDomains()
  {
    auto parent = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    parent.getConnectivity().compute(2, 1);
    parent.getConnectivity().compute(1, 0);
    const LocalMesh unrelated(parent);
    SubMesh<Context::Local>::Builder builder;
    builder.initialize(parent).include(2, 0);
    for (const Index edge : parent.getConnectivity().getIncidence(2, 1).at(0))
      builder.include(1, edge);
    auto submesh = builder.finalize();
    submesh.getConnectivity().compute(2, 1);
    submesh.getConnectivity().compute(1, 0);
    const auto scalarSpace = [&] {
      if constexpr (Degree == 1)
        return P1<Real>(submesh);
      else
        return H1<Degree, Real>(std::integral_constant<size_t, Degree>{}, submesh);
    }();
    const auto vectorSpace = [&] {
      if constexpr (Degree == 1)
        return P1<Math::SpatialVector<Real>>(submesh, 2);
      else
      {
        return H1<Degree, Math::SpatialVector<Real>>(
          std::integral_constant<size_t, Degree>{}, submesh, 2);
      }
    }();
    GridFunction scalar(scalarSpace);
    GridFunction vector(vectorSpace);
    scalar = RealFunction([](const Point& p) { return p.x() + 2 * p.y(); });
    vector = VectorFunction(size_t(2),
      [](const Point& p) { return Math::SpatialVector<Real>{p.x(), 2 * p.y()}; });
    const auto gradient = Grad(scalar);
    const auto derivative = Derivative(size_t(0), scalar);
    const auto jacobian = Jacobian(vector);
    const auto divergence = Div(vector);
    const auto check = [&](const auto& location, bool included) {
      if (included)
      {
        EXPECT_NEAR(gradient.getValue(location)[0], 1, 2e-12);
        EXPECT_NEAR(gradient.getValue(location)[1], 2, 2e-12);
        EXPECT_NEAR(derivative.getValue(location), 1, 2e-12);
        EXPECT_NEAR(jacobian.getValue(location)(1, 1), 2, 2e-12);
        EXPECT_NEAR(divergence.getValue(location), 3, 2e-12);
      }
      else
      {
        EXPECT_THROW(gradient.getValue(location), Alert::Exception);
        EXPECT_THROW(derivative.getValue(location), Alert::Exception);
        EXPECT_THROW(jacobian.getValue(location), Alert::Exception);
        EXPECT_THROW(divergence.getValue(location), Alert::Exception);
      }
    };
    for (const auto* mesh : std::array<const LocalMesh*, 2>{&parent, &unrelated})
    {
      for (size_t cell = 0; cell < mesh->getCellCount(); ++cell)
      {
        const auto entity = mesh->getCell(cell);
        const auto& formula = QF::QuadratureFormula::get(3, entity->getGeometry());
        const auto& quadrature = entity->getQuadrature(formula);
        const Point& point = quadrature.getPoint(0);
        const bool included = mesh == &parent && cell == 0;
        check(point, included);
        check(IntegrationPoint(point), included);
        check(IntegrationPoint(quadrature, 0), included);
      }
    }
    // A non-submesh field must reject a wholly unrelated mesh too.
    P1<Real> wholeSpace(parent);
    GridFunction whole(wholeSpace);
    whole = RealFunction(1);
    const Point foreign(*unrelated.getCell(), Math::SpatialPoint{0.2, 0.3});
    EXPECT_THROW(Grad(whole).getValue(foreign), Alert::Exception);
    EXPECT_THROW(Grad(whole).getValue(IntegrationPoint(foreign)), Alert::Exception);
  }
}

TEST(DerivativeDomain, P1)
{
  checkDerivativeDomains<1>();
}
TEST(DerivativeDomain, H1)
{
  checkDerivativeDomains<2>();
}
