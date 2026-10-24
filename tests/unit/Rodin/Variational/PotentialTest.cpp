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
  /// @brief Scalar evaluation integrates a regular kernel without uninitialized state.
  TEST(Potential, ScalarEvaluation)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = mesh.getCell();
    const Point point(*cell, Math::SpatialPoint{0.2, 0.3});
    auto kernel = [](const Point&, const Point&) { return Real(2); };
    RealFunction constant(Real(3));
    auto potential = Potential(kernel, constant);
    EXPECT_NEAR(potential.getValue(point), 6, 1e-13);
    RealFunction affine([](const Point& source) { return source.x() + source.y(); });
    auto affinePotential = Potential(kernel, affine);
    affinePotential.setOrder(3);
    EXPECT_NEAR(affinePotential.getValue(point), 2, 1e-13);
    auto copied = affinePotential;
    auto moved = std::move(affinePotential);
    EXPECT_NEAR(copied.getValue(point), 2, 1e-13);
    EXPECT_NEAR(moved.getValue(point), 2, 1e-13);
    ComplexFunction complex(Complex(1, 2));
    EXPECT_NEAR(
      std::abs(Potential(kernel, complex).getValue(point) - Complex(2, 4)), 0, 1e-13);
  }

  /// @brief Cloned global integrators own temporary stateful kernels.
  TEST(Potential, TemporaryKernelOwnership)
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(2)
                  .nodes(3)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({0, 1})
                  .polytope(Polytope::Type::Triangle, {0, 1, 2})
                  .finalize();
    P1<Real> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    struct Kernel
    {
      Real scale;
      Real operator()([[maybe_unused]] const Point& target,
        [[maybe_unused]] const Point& source) const
      {
        return scale;
      }
    };
    auto integral = Integral(Potential(Kernel{2}, u), v);
    auto copied = integral;
    const auto cell = mesh.getCell();
    integral.setPolytope(*cell, *cell);
    copied.setPolytope(*cell, *cell);
    EXPECT_NEAR(integral.integrate(0, 0), Real(1) / 18, 1e-13);
    EXPECT_NEAR(copied.integrate(0, 0), Real(1) / 18, 1e-13);
    Kernel kernel{2};
    auto expression = Potential(kernel, u);
    EXPECT_NE(&expression.getKernel(), &kernel);
    kernel.scale = 5;
    EXPECT_EQ(expression.getKernel().scale, 2);
  }

  /// @brief Unimplemented vector evaluation raises an Alert in optimized builds too.
  TEST(Potential, UnsupportedFunctionRange)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = mesh.getCell();
    const Point point(*cell, Math::SpatialPoint{0.2, 0.3});
    auto kernel = [](const Point&, const Point&) { return Real(1); };
    VectorFunction vector{1., 2.};
    auto potential = Potential(kernel, vector);
    EXPECT_THROW(potential.getValue(point), Alert::Exception);
  }
}
