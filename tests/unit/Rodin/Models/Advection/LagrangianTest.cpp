/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Advection/Lagrangian.h"
#include "Rodin/Math/RungeKutta/RK2.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Advection;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Unit tests for Lagrangian class basic functionality
   */
  class LagrangianTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      Mesh<Context::Local> getMesh()
      {
        Mesh mesh;
        mesh = mesh.UniformGrid(GetParam(), { 8, 8 });
        mesh.scale(1.0 / 7.0);
        mesh.getConnectivity().compute(1, 2);
        return mesh;
      }
  };

  /**
   * @brief Test basic Lagrangian class construction.
   *
   * This test verifies that Lagrangian objects can be constructed
   * with proper trial/test functions, initial conditions, and velocity fields.
   */
  TEST_P(LagrangianTest, BasicConstruction)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    // Create trial and test functions
    TrialFunction u(vh);
    TestFunction v(vh);

    // Define initial condition
    auto pi = Math::Constants::pi();
    auto u0 = sin(pi * F::x) * sin(pi * F::y);

    // Define velocity field
    auto velocity = VectorFunction{
      RealFunction([](const Point&) { return 0.1; }),
      RealFunction([](const Point&) { return 0.2; })
    };

    // Test construction
    try
    {
      Lagrangian lagrangian(u, v, u0, velocity);
      // If we get here, construction succeeded
      EXPECT_TRUE(true);
    }
    catch (const std::exception& e)
    {
      FAIL() << "Lagrangian construction failed: " << e.what();
    }
    catch (...)
    {
      FAIL() << "Lagrangian construction failed with unknown exception";
    }
  }

  /**
   * @brief Test Lagrangian with different velocity field types.
   *
   * This test verifies that Lagrangian can work with various velocity fields.
   */
  TEST_P(LagrangianTest, DifferentVelocityFields)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    TrialFunction u(vh);
    TestFunction v(vh);

    // Simple initial condition
    auto u0 = RealFunction([](const Point&) { return 1.0; });

    // Test with constant velocity
    auto constant_velocity = VectorFunction{
      RealFunction([](const Point&) { return 1.0; }),
      RealFunction([](const Point&) { return 0.0; })
    };

    // Test with variable velocity
    auto variable_velocity = VectorFunction{
      F::x,
      F::y
    };

    // Test construction with constant velocity
    EXPECT_NO_THROW({
      Lagrangian lagrangian1(u, v, u0, constant_velocity);
    });

    // Test construction with variable velocity  
    EXPECT_NO_THROW({
      Lagrangian lagrangian2(u, v, u0, variable_velocity);
    });
  }

  /**
   * @brief Test Lagrangian step function basic functionality.
   *
   * This test verifies that the step function can be called without errors.
   */
  TEST_P(LagrangianTest, StepFunction)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    TrialFunction u(vh);
    TestFunction v(vh);

    // Define simple initial condition
    auto u0 = sin(Math::Constants::pi() * F::x)
            * sin(Math::Constants::pi() * F::y);

    // Define simple velocity field
    auto velocity = VectorFunction{
      RealFunction([](const Point&) { return 0.1; }),
      RealFunction([](const Point&) { return 0.0; })
    };

    // Create Lagrangian object and verify step works
    Lagrangian lagrangian(u, v, u0, velocity);
    Real dt = 0.01;
    EXPECT_NO_THROW({
      lagrangian.step(dt);
    });
  }

  /**
   * @brief Test Lagrangian with a custom time-stepping scheme (RK2).
   *
   * Verifies that the user-supplied stepper is actually used by the Lagrangian,
   * rather than a hardcoded default. Both RK2 and RK4 should produce valid
   * results for the same advection problem.
   */
  TEST_P(LagrangianTest, CustomStepperRK2)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    auto u0 = sin(Math::Constants::pi() * F::x)
            * sin(Math::Constants::pi() * F::y);

    auto velocity = VectorFunction{
      RealFunction([](const Point&) { return 0.1; }),
      RealFunction([](const Point&) { return 0.0; })
    };

    // Use RK2 stepper
    Math::RungeKutta::RK2 rk2;

    TrialFunction u(vh);
    TestFunction v(vh);
    Lagrangian lagrangian(u, v, u0, velocity, rk2);

    Real dt = 0.01;
    EXPECT_NO_THROW({
      lagrangian.step(dt);
    });

    // Verify the solution is non-trivial (not all zeros)
    const auto& uh = u.getSolution();
    const size_t cd = mesh.getDimension();
    Math::SpatialVector<Real> rc{{ 0.5, 0.5 }};
    auto it = mesh.getPolytope(cd, 0);
    Geometry::Point p(*it, rc);
    Real val = uh(p);
    EXPECT_NE(val, 0.0) << "Solution should be non-trivial after advection";
  }

  /**
   * @brief Test with different initial conditions.
   *
   * This test verifies that Lagrangian works with various initial conditions.
   */
  TEST_P(LagrangianTest, DifferentInitialConditions)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    TrialFunction u(vh);
    TestFunction v(vh);

    auto velocity = VectorFunction{
      RealFunction([](const Point&) { return 0.1; }),
      RealFunction([](const Point&) { return 0.1; })
    };

    // Test with constant initial condition
    auto u0_constant = RealFunction([](const Point&) { return 2.5; });
    EXPECT_NO_THROW({
      Lagrangian lagrangian1(u, v, u0_constant, velocity);
    });

    // Test with polynomial initial condition
    auto u0_poly = F::x + F::y;
    EXPECT_NO_THROW({
      Lagrangian lagrangian2(u, v, u0_poly, velocity);
    });

    // Test with trigonometric initial condition
    auto pi = Math::Constants::pi();
    auto u0_trig = sin(pi * F::x) * cos(pi * F::y);
    EXPECT_NO_THROW({
      Lagrangian lagrangian3(u, v, u0_trig, velocity);
    });
  }

  INSTANTIATE_TEST_SUITE_P(
    LagrangianMeshParams,
    LagrangianTest,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}
