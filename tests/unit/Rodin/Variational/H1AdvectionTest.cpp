/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <memory>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief Mixed-order vector transport integrates an exact polynomial on a scaled interval.
  TEST(H1Advection, MixedOrdersAndComponentsOnScaledIntervals)
  {
    constexpr Real tolerance = 1e-10;
    for (size_t nodes : {2, 5})
    {
      SCOPED_TRACE(nodes);
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {nodes});
      for (auto vertex = mesh.getVertex(); vertex; ++vertex)
      {
        auto coordinate = vertex->getCoordinates();
        coordinate(0) = 3 * coordinate(0) / static_cast<Real>(nodes - 1) - 2;
        mesh.setVertexCoordinates(vertex->getIndex(), coordinate);
      }
      H1<2, Math::SpatialVector<Real>, LocalMesh> trialSpace(
        std::integral_constant<size_t, 2>{}, mesh, 2);
      H1<3, Math::SpatialVector<Real>, LocalMesh> testSpace(
        std::integral_constant<size_t, 3>{}, mesh, 2);
      const auto square = RealFunction([](const Point& p) {
        const Real s = (p(0) + 2) / 3;
        return s * s;
      });
      const auto firstTest =
        RealFunction([](const Point& p) { return 1 + (p(0) + 2) / 3; });
      const auto secondTest =
        RealFunction([](const Point& p) { return 2 - (p(0) + 2) / 3; });
      VectorFunction beta(size_t{1}, [](const Point& p) {
        Math::SpatialVector<Real> value(1);
        value[0] = 1 + 2 * (p(0) + 2) / 3;
        return value;
      });
      TrialFunction u(trialSpace);
      TestFunction v(testSpace);
      GridFunction state(trialSpace);
      GridFunction weight(testSpace);
      state.project(VectorFunction{square, 3 * square});
      weight.project(VectorFunction{firstTest, secondTest});
      BilinearForm transport(u, v);
      auto integral = Integral(Dot(Jacobian(u) * beta, v));
      static_assert(decltype(integral)::Specialized);
      integral.setOrder(6);
      transport = integral;
      transport.assemble();
      ASSERT_EQ(transport.getOperator().cols(), state.getData().size());
      ASSERT_EQ(transport.getOperator().rows(), weight.getData().size());
      // In s=(x+2)/3, the integrand times dx is 14s+24s^2-8s^3.
      EXPECT_NEAR(
        (transport.getOperator() * state.getData()).dot(weight.getData()), 13, tolerance);
      state.project(VectorFunction{square, RealFunction(0)});
      weight.project(VectorFunction{RealFunction(0), secondTest});
      EXPECT_NEAR(
        (transport.getOperator() * state.getData()).dot(weight.getData()), 0, tolerance);
    }
  }

  /// @brief Moving a bound transport rule preserves its cached matrix and supports rebinding.
  TEST(H1Advection, BoundRuleMoveCopyAndRebindPreservePolynomialAction)
  {
    constexpr Real tolerance = 1e-10;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      auto coordinate = vertex->getCoordinates();
      coordinate(0) *= Real(0.5);
      mesh.setVertexCoordinates(vertex->getIndex(), coordinate);
    }
    H1<2, Math::SpatialVector<Real>, LocalMesh> space(
      std::integral_constant<size_t, 2>{}, mesh, 2);
    TrialFunction u(space);
    TestFunction v(space);
    GridFunction state(space);
    GridFunction weight(space);
    const auto square = RealFunction([](const Point& p) { return p(0) * p(0); });
    state.project(VectorFunction{square, 3 * square});
    weight.project(VectorFunction{RealFunction([](const Point& p) { return 1 + p(0); }),
      RealFunction([](const Point& p) { return 2 - p(0); })});
    VectorFunction beta(size_t{1}, [](const Point& p) {
      Math::SpatialVector<Real> value(1);
      value[0] = 1 + 2 * p(0);
      return value;
    });
    auto rule = Integral(Dot(Jacobian(u) * beta, v));
    static_assert(decltype(rule)::Specialized);
    rule.setOrder(6);
    auto cell = mesh.getCell();
    rule.setPolytope(*cell);
    EXPECT_EQ(&rule.getPolytope(), &*cell);
    auto moved = std::move(rule);
    EXPECT_EQ(&moved.getPolytope(), &*cell);
    const auto action = [&](auto& integrator, const auto& polytope) {
      const auto dofs = space.getDOFs(polytope.getDimension(), polytope.getIndex());
      Real value = 0;
      for (Index i = 0; i < static_cast<Index>(dofs.size()); ++i)
      {
        for (Index j = 0; j < static_cast<Index>(dofs.size()); ++j)
          value += state[dofs[i]] * integrator.integrate(i, j) * weight[dofs[j]];
      }
      return value;
    };
    // The antiderivative is 7s^2+8s^3-2s^4, hence 42/16 on [0,1/2].
    EXPECT_NEAR(action(moved, *cell), Real(42) / 16, tolerance);
    std::unique_ptr<decltype(moved)> copied(moved.copy());
    ++cell;
    copied->setPolytope(*cell);
    EXPECT_EQ(&copied->getPolytope(), &*cell);
    EXPECT_NEAR(action(*copied, *cell), Real(166) / 16, tolerance);
    moved.setPolytope(*cell);
    EXPECT_NEAR(action(moved, *cell), Real(166) / 16, tolerance);
  }
}
