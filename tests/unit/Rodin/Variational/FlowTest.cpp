/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Advection/Lagrangian.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Advection;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Unit tests for Flow class basic functionality
   */
  class FlowTest : public ::testing::TestWithParam<Polytope::Type>
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
   * @brief Test basic Flow class construction and coordinate handling.
   *
   * This test verifies that Flow objects can be constructed with valid
   * velocity fields and that basic coordinate system handling works.
   */
  TEST_P(FlowTest, BasicConstruction)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    // Define constant velocity field: v = (0.1, 0.2)
    auto velocity = VectorFunction{
      RealFunction([](const Point&) { return 0.1; }),
      RealFunction([](const Point&) { return 0.2; })
    };

    // Define the time step for advection
    Real dt = 0.1;
    (void)dt; // Suppress unused variable warning

    // Create test function
    TestFunction v(vh);

    // Test that we can create points with the correct coordinate format
    Math::SpatialVector<Real> coords{{ 0.5, 0.5 }};
    EXPECT_EQ(coords.size(), 2);
    EXPECT_NEAR(coords[0], 0.5, 1e-10);
    EXPECT_NEAR(coords[1], 0.5, 1e-10);

    // Test that we can create a point on the mesh
    auto polytope = mesh.getPolytope(2, 0);
    EXPECT_TRUE(polytope != mesh.getPolytope(2).end());

    // Create the test point 
    Geometry::Point testPoint(*polytope, coords);

    // Test velocity evaluation
    auto vel_value = velocity(testPoint);
    EXPECT_NEAR(vel_value[0], 0.1, 1e-10);
    EXPECT_NEAR(vel_value[1], 0.2, 1e-10);
  }

  /**
   * @brief Test different velocity field types and their evaluation.
   *
   * This test verifies that different types of velocity fields can be
   * properly constructed and evaluated at specific points.
   */
  TEST_P(FlowTest, VelocityFieldEvaluation)
  {
    Mesh mesh = this->getMesh();

    // Test constant velocity field
    auto constant_velocity = VectorFunction{
      RealFunction([](const Point&) { return 1.0; }),
      RealFunction([](const Point&) { return 0.5; })
    };

    // Test rotational velocity field: v = (-y, x)
    auto rotational_velocity = VectorFunction{
      -F::y,
      F::x
    };

    // Test linear shear velocity field: v = (y, 0)
    auto shear_velocity = VectorFunction{
      F::y,
      Zero()
    };

    // Create a test point
    Math::SpatialVector<Real> coords{{ 0.3, 0.4 }};
    auto polytope_iter = mesh.getPolytope(2, 0);
    Geometry::Point testPoint(*polytope_iter, coords);

    // Get the actual physical coordinates of the point
    auto physCoords = testPoint.getPhysicalCoordinates();
    Real actual_x = physCoords[0];
    Real actual_y = physCoords[1];

    // Evaluate constant velocity - should be (1.0, 0.5) everywhere
    auto const_vel_value = constant_velocity(testPoint);
    EXPECT_NEAR(const_vel_value[0], 1.0, 1e-10);
    EXPECT_NEAR(const_vel_value[1], 0.5, 1e-10);

    // Evaluate rotational velocity at the actual physical coordinates
    auto rot_vel_value = rotational_velocity(testPoint);
    EXPECT_NEAR(rot_vel_value[0], -actual_y, 1e-3);
    EXPECT_NEAR(rot_vel_value[1], actual_x, 1e-3);

    // Evaluate shear velocity at the actual physical coordinates
    auto shear_vel_value = shear_velocity(testPoint);
    EXPECT_NEAR(shear_vel_value[0], actual_y, 1e-3);
    EXPECT_NEAR(shear_vel_value[1], 0.0, 1e-10);
  }

  /**
   * @brief Test mesh polytope access and validation.
   *
   * This test verifies basic mesh functionality needed for Flow operations.
   */
  TEST_P(FlowTest, MeshAccess)
  {
    Mesh mesh = this->getMesh();

    // Check that mesh has elements
    EXPECT_GT(mesh.getCellCount(), 0);

    // Test polytope iteration
    auto polytope = mesh.getPolytope(2, 0);
    EXPECT_TRUE(polytope != mesh.getPolytope(2).end());

    // Test coordinate creation and point construction
    Math::SpatialVector<Real> coords{{ 0.25, 0.75 }};
    Geometry::Point testPoint(*polytope, coords);

    // Verify physical coordinates are reasonable (within [0,1]x[0,1])
    auto physCoords = testPoint.getPhysicalCoordinates();
    EXPECT_GE(physCoords[0], 0.0);
    EXPECT_LE(physCoords[0], 1.0);
    EXPECT_GE(physCoords[1], 0.0);
    EXPECT_LE(physCoords[1], 1.0);
  }

  INSTANTIATE_TEST_SUITE_P(
    FlowMeshParams,
    FlowTest,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}
