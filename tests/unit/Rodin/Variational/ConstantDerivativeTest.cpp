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

namespace Rodin::Tests::Unit
{
  /// @brief Broken P0 gradients report zero order over both scalar fields.
  TEST(ConstantDerivative, P0GradientOrder)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = mesh.getCell();
    const Point point(*cell, Math::SpatialPoint{0.2, 0.3});
    const auto check = [&]<class Scalar>() {
      P0<Scalar> space(mesh);
      GridFunction field(space);
      field.getData().setConstant(Scalar(2));
      auto gradient = Grad(field);
      EXPECT_EQ(gradient.getOrder(*cell), Optional<size_t>(0));
      EXPECT_EQ(gradient.getValue(point).norm(), 0);
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
  }

  /// @brief Constant vector derivatives retain their leaves in assembled expressions.
  TEST(ConstantDerivative, GlobalConstantShapeLeaves)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = mesh.getCell();
    const auto check = [&]<class Scalar>() {
      P0g<Math::SpatialVector<Scalar>, LocalMesh> vectors(mesh, 2);
      P0g<Scalar, LocalMesh> scalars(mesh);
      TrialFunction u(vectors);
      TestFunction v(vectors);
      TestFunction q(scalars);
      auto divergence = Integral(Div(u), q);
      auto jacobian = Integral(Jacobian(u), Jacobian(v));
      divergence.setPolytope(*cell);
      jacobian.setPolytope(*cell);
      for (size_t trial = 0; trial < 2; ++trial)
      {
        EXPECT_EQ(divergence.integrate(trial, 0), Scalar(0));
        for (size_t test = 0; test < 2; ++test)
          EXPECT_EQ(jacobian.integrate(trial, test), Scalar(0));
      }
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
  }
}
