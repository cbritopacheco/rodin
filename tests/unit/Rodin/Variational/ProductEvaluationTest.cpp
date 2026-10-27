/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <memory>
#include <type_traits>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief Coefficients on the right stay live when evaluation has no quadrature metadata.
  TEST(ProductEvaluation, RightCoefficientsWithoutQuadratureMetadataRemainLive)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = *mesh.getCell();
    const Math::SpatialPoint reference{0, 0};
    const Point point(cell, reference);
    const IntegrationPoint ip(point);
    const auto check = [&]<class Range>(const P1<Range>& space) {
      size_t calls = 0;
      Real scale = 2;
      auto coefficient = RealFunction([&](const Point& location) {
        ++calls;
        return scale + location(0);
      });
      static_assert(decltype(coefficient)::Cache::Enabled);
      TestFunction v(space);
      auto expression = v * coefficient;
      expression.setIntegrationPoint(ip);
      EXPECT_EQ(calls, 0u);
      const auto verify = [](const auto& value, Real expected) {
        if constexpr (std::is_same_v<Range, Real>)
          EXPECT_EQ(value, expected);
        else
        {
          ASSERT_EQ(value.size(), 2u);
          EXPECT_EQ(value[0], expected);
          EXPECT_EQ(value[1], 0);
        }
      };
      verify(expression.getBasis(0), 2);
      EXPECT_EQ(calls, 1u);
      scale = 5;
      verify(expression.getBasis(0), 5);
      EXPECT_EQ(calls, 2u);
      auto copied = expression;
      auto moved = std::move(expression);
      copied.setIntegrationPoint(ip);
      moved.setIntegrationPoint(ip);
      EXPECT_EQ(calls, 2u);
      scale = 7;
      verify(copied.getBasis(0), 7);
      verify(moved.getBasis(0), 7);
      EXPECT_EQ(calls, 4u);
    };
    P1<Real> scalars(mesh);
    P1<Math::SpatialVector<Real>> vectors(mesh, 2);
    check(scalars);
    check(vectors);
  }

  /// @brief Scaled integrators expose their current cell and cloned bindings remain independent.
  TEST(ProductEvaluation, ScaledIntegratorCopiesRebindIndependently)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    auto cells = mesh.getCell();
    const auto first = *cells;
    ++cells;
    const auto second = *cells;
    const auto check = [&](auto scaled) {
      scaled.setPolytope(first);
      EXPECT_EQ(&scaled.getPolytope(), &first);
      std::unique_ptr<decltype(scaled)> copied(scaled.copy());
      copied->setPolytope(second);
      EXPECT_EQ(&copied->getPolytope(), &second);
      EXPECT_EQ(&scaled.getPolytope(), &first);
      EXPECT_NE(copied->getPolytope().getIndex(), scaled.getPolytope().getIndex());
    };
    check(Real(2) * Integral(v));
    check(Real(2) * Integral(u, v));
  }
}
