/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Variational/Flow.h"

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"

#include "Rodin/Models/Advection/Lagrangian.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Models::Advection;

/**
 * @brief Manufactured solution tests for Lagrangian advection.
 *
 * These tests verify the mathematical setup and evaluate known analytical 
 * solutions to the advection equation ∂u/∂t + v·∇u = 0.
 * 
 * Note: The current Lagrangian implementation has API issues, so these tests
 * focus on mathematical validation rather than full solver integration.
 */
namespace Rodin::Tests::Manufactured::AdvectionLagrangian
{
  template <size_t M>
  class ManufacturedAdvectionTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      Mesh<Context::Local> getMesh()
      {
        Mesh mesh;
        mesh = mesh.UniformGrid(GetParam(), { M, M });
        mesh.scale(1.0 / (M - 1));
        mesh.getConnectivity().compute(1, 2);
        return mesh;
      }
  };

  using ManufacturedAdvectionTest16x16 = ManufacturedAdvectionTest<16>;
  using ManufacturedAdvectionTest32x32 = ManufacturedAdvectionTest<32>;

  /**
   * @brief Test manufactured solution setup: Translation with constant velocity.
   *
   * Exact solution: u(x,y,t) = sin(π(x-vx*t)) * sin(π(y-vy*t))
   * where vx = 0.5, vy = 0.3 are constant velocities.
   * This test verifies that the basic mathematical setup is correct.
   */
  TEST_P(ManufacturedAdvectionTest16x16, ConstantVelocityTranslation_Setup)
  {
    const auto pi = Math::Constants::pi();
    const auto& mesh = this->getMesh();
    P1 vh(mesh);

    // constant velocity
    const Real vx = 0.5, vy = 0.3;
    auto velocity = VectorFunction{
      RealFunction([vx](const Point&) { return vx; }),
      RealFunction([vy](const Point&) { return vy; })
    };

    // u0(x,y) = sin(pi x) sin(pi y)
    auto u0 = sin(pi * F::x) * sin(pi * F::y);

    // exact at time t: u(x,y,t) = sin(pi(x - vx t)) sin(pi(y - vy t))
    const Real t_final = 0;
    auto u_exact = sin(pi * (F::x - vx * t_final)) * sin(pi * (F::y - vy * t_final));

    // trial / test
    TrialFunction u(vh);
    TestFunction v(vh);

    // construct Lagrangian object
    Models::Advection::Lagrangian lagrangian(u, v, u0, velocity);
    lagrangian.step(1);

    u.getSolution().save("u0.gf");
    mesh.save("u0.mesh");
    std::exit(1);

    // choose a concrete cell and a reference point; then use PHYSICAL coords
    auto polytope = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::SpatialVector<Real> rc{{0.5, 0.5}};                  // reference coords
    Point p(*polytope, rc);
    const auto pc = p.getPhysicalCoordinates();              // physical coords
    const Real x = pc[0], y = pc[1];

    // evaluate fields
    const Real u0_val = u0(p);
    const Real uex_val = u_exact(p);
    const auto vval = velocity(p);

    // sanity
    EXPECT_TRUE(std::isfinite(u0_val));
    EXPECT_TRUE(std::isfinite(uex_val));
    EXPECT_LE(std::abs(u0_val), 1.0);
    EXPECT_LE(std::abs(uex_val), 1.0);

    // velocity is constant
    EXPECT_NEAR(vval[0], vx, 1e-12);
    EXPECT_NEAR(vval[1], vy, 1e-12);

    // periodic wrap on [0,1]^2
    auto wrap01 = [](Real z) {
      // works for any real z
      z -= std::floor(z);
      return z;
    };

    const Real x0 = wrap01(x - vx * t_final);
    const Real y0 = wrap01(y - vy * t_final);
    const Real expected = std::sin(pi * x0) * std::sin(pi * y0);

    // advection identity under periodic BCs
    EXPECT_NEAR(uex_val, expected, 1e-12);
  }

  // /**
  //  * @brief Test manufactured solution setup: Rotation about origin.
  //  *
  //  * Velocity field: v = ω*(-y, x) for rotation with angular velocity ω.
  //  * This tests the setup for solid body rotation problems.
  //  */
  // TEST_P(ManufacturedAdvectionTest16x16, RigidBodyRotation_Setup)
  // {
  //   auto pi = Math::Constants::pi();
  //   Mesh mesh = this->getMesh();
  //   P1 vh(mesh);

  //   // Angular velocity
  //   Real omega = 1.0;
  //   
  //   // Rotational velocity field: v = ω*(-y, x)
  //   auto velocity = VectorFunction{
  //     -omega * F::y,
  //     omega * F::x
  //   };

  //   // Initial condition: Gaussian-like profile centered at (0.75, 0.5)
  //   auto u0 = exp(-10.0 * ((F::x - 0.75) * (F::x - 0.75) + (F::y - 0.5) * (F::y - 0.5)));

  //   // Create trial and test functions
  //   TrialFunction u(vh);
  //   TestFunction v(vh);

  //   // Test that Lagrangian can be constructed
  //   EXPECT_NO_THROW({
  //     Lagrangian lagrangian(u, v, u0, velocity);
  //   });

  //   // Test rotational velocity field properties
  //   std::vector<std::pair<Real, Real>> test_points = {
  //     {0.5, 0.5}, {0.3, 0.7}, {0.8, 0.2}
  //   };

  //   for (const auto& [x_test, y_test] : test_points)
  //   {
  //     Math::SpatialVector<Real> coords{{ x_test, y_test }};
  //     auto polytope = mesh.getPolytope(2, 0);
  //     Geometry::Point testPoint(*polytope, coords);
  //     
  //     auto vel_value = velocity(testPoint);
  //     Real expected_vx = -omega * y_test;
  //     Real expected_vy = omega * x_test;
  //     
  //     EXPECT_NEAR(vel_value[0], expected_vx, 1e-3);
  //     EXPECT_NEAR(vel_value[1], expected_vy, 1e-3);
  //     
  //     // Test that velocity is tangent to circles centered at origin
  //     Real vel_magnitude = std::sqrt(vel_value[0]*vel_value[0] + vel_value[1]*vel_value[1]);
  //     Real radius = std::sqrt(x_test*x_test + y_test*y_test);
  //     Real expected_magnitude = omega * radius;
  //     EXPECT_NEAR(vel_magnitude, expected_magnitude, 1e-10);
  //   }
  // }

  // /**
  //  * @brief Test manufactured solution setup: Linear shear flow.
  //  *
  //  * Velocity field: v = (y, 0) - linear shear in x-direction
  //  * Initial condition: u₀(x,y) = sin(2πx)
  //  * Exact solution: u(x,y,t) = sin(2π(x - yt))
  //  */
  // TEST_P(ManufacturedAdvectionTest32x32, LinearShearFlow_Setup)
  // {
  //   auto pi = Math::Constants::pi();
  //   Mesh mesh = this->getMesh();
  //   P1 vh(mesh);

  //   // Linear shear velocity field: v = (y, 0)
  //   auto velocity = VectorFunction{
  //     F::y,
  //     Zero()
  //   };

  //   // Initial condition: u₀(x,y) = sin(2π*x)
  //   auto u0 = sin(2.0 * pi * F::x);

  //   // Exact solution at time t: u(x,y,t) = sin(2π*(x - y*t))
  //   Real t_final = 0.1;
  //   auto u_exact = sin(2.0 * pi * (F::x - F::y * t_final));

  //   // Create trial and test functions
  //   TrialFunction u(vh);
  //   TestFunction v(vh);

  //   // Test that Lagrangian can be constructed
  //   EXPECT_NO_THROW({
  //     Lagrangian lagrangian(u, v, u0, velocity);
  //   });

  //   // Validate exact solution at test points
  //   std::vector<std::pair<Real, Real>> test_points = {
  //     {0.25, 0.25}, {0.5, 0.5}, {0.75, 0.25}
  //   };

  //   for (const auto& [x_test, y_test] : test_points)
  //   {
  //     Math::SpatialVector<Real> coords{{ x_test, y_test }};
  //     auto polytope = mesh.getPolytope(2, 0);
  //     Geometry::Point testPoint(*polytope, coords);
  //     
  //     Real initial_value = u0(testPoint);
  //     Real exact_value = u_exact(testPoint);
  //     auto vel_value = velocity(testPoint);
  //     
  //     EXPECT_TRUE(std::isfinite(initial_value));
  //     EXPECT_TRUE(std::isfinite(exact_value));
  //     EXPECT_TRUE(std::abs(initial_value) <= 1.0); // sin function is bounded
  //     EXPECT_TRUE(std::abs(exact_value) <= 1.0);   // sin function is bounded
  //     
  //     // Velocity should be (y, 0)
  //     EXPECT_NEAR(vel_value[0], y_test, 1e-10);
  //     EXPECT_NEAR(vel_value[1], 0.0, 1e-10);
  //     
  //     // Test shear property: u(x,y,t) should equal u₀(x - y*t, y)
  //     Real x_sheared = x_test - y_test * t_final;
  //     Real expected_exact = std::sin(2.0 * pi * x_sheared);
  //     EXPECT_NEAR(exact_value, expected_exact, 1e-10);
  //   }
  // }

  // /**
  //  * @brief Test velocity field conservation properties.
  //  *
  //  * This test verifies properties that should be preserved by the advection equation,
  //  * such as mass conservation in divergence-free velocity fields.
  //  */
  // TEST_P(ManufacturedAdvectionTest32x32, VelocityFieldProperties)
  // {
  //   Mesh mesh = this->getMesh();
  //   P1 vh(mesh);

  //   // Test divergence-free rotational field
  //   Real omega = 0.5;
  //   auto rotation_velocity = VectorFunction{
  //     -omega * F::y,
  //     omega * F::x
  //   };

  //   // Test constant velocity field (also divergence-free)
  //   auto constant_velocity = VectorFunction{
  //     RealFunction([](const Point&) { return 1.0; }),
  //     RealFunction([](const Point&) { return 0.5; })
  //   };

  //   // Test velocity field evaluations
  //   std::vector<std::pair<Real, Real>> test_points = {
  //     {0.2, 0.3}, {0.6, 0.8}, {0.9, 0.1}
  //   };

  //   for (const auto& [x_test, y_test] : test_points)
  //   {
  //     Math::SpatialVector<Real> coords{{ x_test, y_test }};
  //     auto polytope = mesh.getPolytope(2, 0);
  //     Geometry::Point testPoint(*polytope, coords);
  //     
  //     // Test rotational velocity
  //     auto rot_vel = rotation_velocity(testPoint);
  //     EXPECT_NEAR(rot_vel[0], -omega * y_test, 1e-10);
  //     EXPECT_NEAR(rot_vel[1], omega * x_test, 1e-10);
  //     
  //     // Test constant velocity
  //     auto const_vel = constant_velocity(testPoint);
  //     EXPECT_NEAR(const_vel[0], 1.0, 1e-10);
  //     EXPECT_NEAR(const_vel[1], 0.5, 1e-10);
  //     
  //     // Both fields should be finite
  //     EXPECT_TRUE(std::isfinite(rot_vel[0]) && std::isfinite(rot_vel[1]));
  //     EXPECT_TRUE(std::isfinite(const_vel[0]) && std::isfinite(const_vel[1]));
  //   }

  //   // Test construction with different velocity fields
  //   TrialFunction u(vh);
  //   TestFunction v(vh);
  //   auto simple_initial = RealFunction([](const Point&) { return 1.0; });

  //   EXPECT_NO_THROW({
  //     Lagrangian lagrangian1(u, v, simple_initial, rotation_velocity);
  //     Lagrangian lagrangian2(u, v, simple_initial, constant_velocity);
  //   });
  // }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    ManufacturedAdvectionTest16x16,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

  INSTANTIATE_TEST_SUITE_P(
    MeshParams32x32,
    ManufacturedAdvectionTest32x32,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}
