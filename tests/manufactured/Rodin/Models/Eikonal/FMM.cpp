/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <gtest/gtest.h>

#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/VectorFunction.h"
#include "Rodin/Models/Eikonal/FMM.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Manufactured::Eikonal
{
  /**
   * @brief Manufactured solution tests for Fast Marching Method (FMM)
   * 
   * These tests verify the FMM implementation against known analytical solutions
   * for various geometric configurations and speed functions.
   */
  class FMMManufacturedTest : public ::testing::Test
  {
    protected:
      static constexpr Real TOLERANCE = 1e-2;  // Relaxed tolerance for manufactured solutions
      static constexpr Real TIGHT_TOLERANCE = 1e-3;  // Tighter tolerance for simple cases
  };

  // Test 1: Point source with constant speed - 2D Euclidean distance
  TEST_F(FMMManufacturedTest, PointSource_ConstantSpeed_2D_EuclideanDistance)
  {
    const size_t n = 33;
    const Real h = 2.0 / n;
    const Real tol = h;
    const Real exclude_r = 2.0 * h;

    // Create fine 2D triangular mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.scale(2.0 / n);  // Scale to [0, 2] x [0, 2] 
    mesh.displace(VectorFunction{ -1.0, -1.0 });  // Center at origin: [-1, 1] x [-1, 1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed = 1
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    // Set source at origin (0, 0)
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance = coord.norm();
      if (distance < 0.05)
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    Models::Eikonal::FMM fmm(u, speed);
    fmm.seed(interface).solve();

    // Verify against analytical solution: u(x,y) = sqrt(x^2 + y^2)
    Real max_error = 0.0;
    Real total_error = 0.0;
    Index count = 0;

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto x = mesh.getVertexCoordinates(it->getIndex());
      if (x.norm() < exclude_r)
        continue; // skip nodes near source
      Real err = std::abs(u[it->getIndex()] - x.norm());
      max_error = std::max(max_error, err);
      total_error += err; ++count;
    }

    ASSERT_GT(count, 0);
    EXPECT_LT(max_error, tol);
    EXPECT_LT(total_error / count, 0.5 * tol);
  }

  // Test 2: Point source with variable speed - 2D radial speed function
  TEST_F(FMMManufacturedTest, PointSource_RadialSpeed_2D)
  {
    // Create 2D mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 25, 25 });
    mesh.scale(1.0 / 24.0);  // Scale to [0, 1] x [0, 1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Radial speed function: s(r) = 1 + r
    auto speed = [](const Geometry::Point& p) -> Real 
    {
      Real r = (p.getCoordinates() - Math::SpatialVector<Real>{{0.5, 0.5}}).norm();
      return 1.0 + r;
    };

    // Set source at center (0.5, 0.5)
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance = (coord - Math::SpatialVector<Real>{{0.5, 0.5}}).norm();
      if (distance < 0.03)
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    Models::Eikonal::FMM fmm(u, speed);
    fmm.seed(interface).solve();

    // For radial speed s(r) = 1 + r, the analytical solution is:
    // u(r) = ln(1 + r) (for small r approximation)
    // We'll check that the solution is monotonic and has reasonable values

    Real center_val = std::numeric_limits<Real>::infinity();
    Real edge_val = 0.0;

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real r = (coord - Math::SpatialVector<Real>{{0.5, 0.5}}).norm();
      Real computed = u[it->getIndex()];

      if (r < 0.05)
        center_val = std::min(center_val, computed);
      if (r > 0.4)
        edge_val = std::max(edge_val, computed);
    }

    EXPECT_LT(center_val, TIGHT_TOLERANCE)
      << "Center should have near-zero distance";

    EXPECT_GT(edge_val, 0.2)
      << "Edge should have reasonable distance";

    EXPECT_LT(center_val, edge_val)
      << "Distance should increase with radius";
  }

  // Test 3: Surface mesh test - Box surface with point source
  TEST_F(FMMManufacturedTest, BoxSurface_PointSource)
  {
    // Create surface mesh
    Mesh mesh;
    mesh = mesh.Box(Polytope::Type::Triangle, { 16, 16, 16 });
    mesh.scale(2.0 / 16.0);  // Scale to [0, 2]^3
    mesh.displace(VectorFunction{ -1.0, -1.0, -1.0 });  // Center at origin: [-1, 1]^3
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 2);
    mesh.getConnectivity().compute(0, 0);


    EXPECT_EQ(mesh.getDimension(), 2)
      << "Surface mesh should be 2D";

    EXPECT_EQ(mesh.getSpaceDimension(), 3)
      << "Surface mesh should be embedded in 3D";

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set source near one face center (e.g., x = 1 face center)
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      const Real distance = (coord - Math::SpatialVector<Real>{{ -1.0, -1.0, -1.0 }}).norm();
      if (distance < 0.07)
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    fmm.seed(interface).solve();

    // Verify solution properties on surface
    Real min_dist = std::numeric_limits<Real>::infinity();
    Real max_dist = 0.0;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val))
        << "Solution should not contain NaN";

      EXPECT_GE(val, 0.0)
        << "Distance should be non-negative";

      if (std::isfinite(val))
      {
        min_dist = std::min(min_dist, val);
        max_dist = std::max(max_dist, val);
      }
    }

    EXPECT_LT(min_dist, TIGHT_TOLERANCE)
      << "Minimum distance should be near zero at source";

    EXPECT_GT(max_dist, 1.0)
      << "Maximum distance should span reasonable range on surface";
  }

  // Test 4: Sphere mapping test - Map cube to sphere and test geodesic distances
  TEST_F(FMMManufacturedTest, SphereMappingTest_GeodesicDistance)
  {
    // Create cube mesh
    Mesh mesh;
    mesh = mesh.Box(Polytope::Type::Triangle, { 32, 32, 32 });
    mesh.scale(2.0 / 32.0);  // Scale to [0, 2]^3
    mesh.displace(VectorFunction{ -1.0, -1.0, -1.0 });  // Center at origin: [-1, 1]^3

    // Extract surface
    mesh.getConnectivity().compute(2, 0);

    // Map surface vertices to sphere
    Math::SpatialVector<Real> coord;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      coord = mesh.getVertexCoordinates(it->getIndex());
      Real norm = coord.norm();
      if (norm > 1e-12)  // Avoid division by zero
        coord = coord / norm;  // Project to unit sphere
      mesh.setVertexCoordinates(it->getIndex(), coord);
    }
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed on sphere surface
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set source at one point on the sphere (e.g., (1, 0, 0))
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance = (coord - Math::SpatialVector<Real>{{1.0, 0.0, 0.0}}).norm();
      if (distance < 0.05)
        interface.push_back(it->getIndex());
    }

    EXPECT_FALSE(interface.empty());

    fmm.seed(interface).solve();

    // For a unit sphere, the geodesic distance between two points is the arc length
    // For points (1,0,0) and (-1,0,0), the geodesic distance is π
    Real max_expected_distance = M_PI;
    Real min_dist = std::numeric_limits<Real>::infinity();
    Real max_dist = 0.0;
    Real opposite_point_distance = std::numeric_limits<Real>::infinity();
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real val = u[it->getIndex()];
      if (std::isfinite(val))
      {
        min_dist = std::min(min_dist, val);
        max_dist = std::max(max_dist, val);
        // Check point approximately opposite to source (1,0,0) -> (-1,0,0)
        Real dist_to_opposite = (coord - Math::SpatialVector<Real>{{-1.0, 0.0, 0.0}}).norm();
        if (dist_to_opposite < 0.05)
          opposite_point_distance = std::min(opposite_point_distance, val);
      }
    }

    EXPECT_LT(min_dist, TIGHT_TOLERANCE)
      << "Minimum distance should be near zero at source";

    EXPECT_LT(max_dist, max_expected_distance + TOLERANCE)
      << "Maximum distance should not exceed sphere circumference";

    // Check that opposite point has distance close to π (if it exists in midsection)
    if (std::isfinite(opposite_point_distance))
    {
      Real expected_opposite_dist = M_PI;
      Real error = std::abs(opposite_point_distance - expected_opposite_dist);
      EXPECT_LT(error, 1.0 / 16)
        << "Opposite point should have geodesic distance ≈ π";
    }
  }

  // Test 5: 3D volumetric manufactured test - Constant speed in cube
  TEST_F(FMMManufacturedTest, Volume3D_ConstantSpeed_CubeCenter)
  {
    // Create 3D tetrahedral mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, { 16, 16, 16 });
    mesh.scale(2.0 / 15.0);  // Scale to [0, 2]^3
    mesh.displace(VectorFunction{ -1.0, -1.0, -1.0 });  // Center at origin: [-1, 1]^3
    mesh.getConnectivity().compute(3, 0);
    mesh.getConnectivity().compute(2, 3);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set source at center (0, 0, 0)
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance = coord.norm();
      if (distance < 0.2)
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    fmm.seed(interface).solve();

    // Verify against analytical solution: u(x,y,z) = sqrt(x^2 + y^2 + z^2)
    Real total_error = 0.0;
    Index count = 0;

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real computed = u[it->getIndex()];
      Real analytical = coord.norm();

      if (std::isfinite(computed) && std::isfinite(analytical) && analytical > 0.05)
      {
        Real error = std::abs(computed - analytical);
        total_error += error;
        count++;
      }
    }

    if (count > 0)
    {
      Real avg_error = total_error / count;
      EXPECT_LT(avg_error, 2.0 / 15.0)
        << "Average error should be within tolerance";
    }
  }

  // Test 6: Anisotropic speed function test
  TEST_F(FMMManufacturedTest, AnisotropicSpeed_2D)
  {
    // Create 2D mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 21, 21 });
    mesh.scale(2.0 / 20.0);  // Scale to [0, 2] x [0, 2]
    mesh.displace(VectorFunction{ -1.0, -1.0 });  // Center at origin: [-1, 1] x [-1, 1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Anisotropic speed: faster in x-direction, slower in y-direction
    auto speed = [](const Geometry::Point& p) -> Real
    {
      const auto& coord = p.getCoordinates();
      Real x = coord[0], y = coord[1];
      return 1.0 + 10.0 * std::abs(x) + 1.0 * std::abs(y);  // Varies with position
    };

    Models::Eikonal::FMM fmm(u, speed);

    // Set source at origin
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance = coord.norm();
      if (distance < 0.05)
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    fmm.seed(interface).solve();

    // Verify solution properties
    Real x_positive_dist = std::numeric_limits<Real>::infinity();
    Real y_positive_dist = std::numeric_limits<Real>::infinity();

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real val = u[it->getIndex()];

      // Check point at (0.5, 0)
      if (std::abs(coord[0] - 0.5) < 0.1 && std::abs(coord[1]) < 0.1)
        x_positive_dist = std::min(x_positive_dist, val);

      // Check point at (0, 0.5)  
      if (std::abs(coord[0]) < 0.1 && std::abs(coord[1] - 0.5) < 0.1)
        y_positive_dist = std::min(y_positive_dist, val);
    }

    // Due to anisotropic speed, travel time should be different in different directions
    EXPECT_GT(y_positive_dist, x_positive_dist)
      << "Travel in y-direction should take longer due to slower speed";
  }

  // Test 7: Line source with constant speed - 2D 
  TEST_F(FMMManufacturedTest, LineSource_ConstantSpeed_2D)
  {
    // Create 2D mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 25, 25 });
    mesh.scale(2.0 / 24.0);  // Scale to [0, 2] x [0, 2]
    mesh.displace(VectorFunction{ -1.0, -1.0 });  // Center at origin: [-1, 1] x [-1, 1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set line source along x-axis: y = 0, x ∈ [-0.5, 0.5]
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0], y = coord[1];
      if (std::abs(y) < 0.05 && std::abs(x) < 0.55)  // Line source
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Line interface should not be empty";

    fmm.seed(interface).solve();

    // For a line source along the x-axis, the distance function should be:
    // u(x,y) = |y| for points near the line
    Real total_error = 0.0;
    Index count = 0;

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0], y = coord[1];
      Real computed = u[it->getIndex()];

      // Test analytical solution for points near the line source
      if (std::abs(x) < 0.4 && std::abs(y) > 0.1 && std::abs(y) < 0.8)
      {
        Real analytical = std::abs(y);  // Distance to line y = 0
        if (std::isfinite(computed) && std::isfinite(analytical))
        {
          Real error = std::abs(computed - analytical);
          total_error += error;
          count++;
        }
      }
    }

    if (count > 0)
    {
      Real avg_error = total_error / count;
      EXPECT_LT(avg_error, 0.1)
        << "Average error should be small for line source";
    }

    // Verify general properties
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      if (std::isfinite(val))
      {
        EXPECT_GE(val, 0.0) << "Distance should be non-negative";
      }
    }
  }

  // Test 8: Circular ring source with constant speed
  TEST_F(FMMManufacturedTest, CircularRingSource_ConstantSpeed_2D)
  {
    // Create 2D mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 33, 33 });
    mesh.scale(2.0 / 32.0);  // Scale to [0, 2] x [0, 2]
    mesh.displace(VectorFunction{ -1.0, -1.0 });  // Center at origin: [-1, 1] x [-1, 1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set circular ring source: radius = 0.5, centered at origin
    const Real ring_radius = 0.5;
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance_to_origin = coord.norm();
      if (std::abs(distance_to_origin - ring_radius) < 0.04)  // Ring source
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Circular ring interface should not be empty";

    fmm.seed(interface).solve();

    // For a circular ring source, the analytical solution is:
    // u(r) = |r - ring_radius| where r is distance from origin
    Real total_error = 0.0;
    Index count = 0;

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real r = coord.norm();
      Real computed = u[it->getIndex()];

      // Test points not too close to boundaries and ring
      if (r > 0.1 && r < 0.9 && std::abs(r - ring_radius) > 0.08)
      {
        Real analytical = std::abs(r - ring_radius);
        if (std::isfinite(computed) && std::isfinite(analytical))
        {
          Real error = std::abs(computed - analytical);
          total_error += error;
          count++;
        }
      }
    }

    if (count > 0)
    {
      Real avg_error = total_error / count;
      EXPECT_LT(avg_error, 0.1)
        << "Average error should be small for circular ring source";
    }

    // Check interior and exterior regions
    Real interior_dist = std::numeric_limits<Real>::infinity();
    Real exterior_dist = std::numeric_limits<Real>::infinity();

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real r = coord.norm();
      Real val = u[it->getIndex()];

      if (std::isfinite(val))
      {
        EXPECT_GE(val, 0.0) << "Distance should be non-negative";

        // Interior points (r < ring_radius)
        if (r < ring_radius - 0.1)
          interior_dist = std::min(interior_dist, val);

        // Exterior points (r > ring_radius)  
        if (r > ring_radius + 0.1)
          exterior_dist = std::min(exterior_dist, val);
      }
    }

    EXPECT_LT(interior_dist, 0.3) << "Interior should have reasonable distance";
    EXPECT_LT(exterior_dist, 0.3) << "Exterior should have reasonable distance";
  }

  // Test 9: Boundary-based initial conditions with analytical solution
  TEST_F(FMMManufacturedTest, BoundaryInitialCondition_SquareDomain)
  {
    // Create 2D mesh on [0,1] x [0,1]
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 21, 21 });
    mesh.scale(1.0 / 20.0);  // Scale to [0, 1] x [0, 1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set interface to left boundary (x = 0)
    std::vector<Index> interface;
    const Real boundary_tolerance = 1e-6;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0];
      if (x < boundary_tolerance)  // Left boundary x = 0
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Boundary interface should not be empty";

    fmm.seed(interface).solve();

    // Analytical solution: u(x,y) = x (distance from left boundary)
    Real total_error = 0.0;
    Index count = 0;

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0];
      Real computed = u[it->getIndex()];

      if (std::isfinite(computed) && x > 0.05)  // Exclude boundary region
      {
        Real analytical = x;  // Distance from left boundary
        Real error = std::abs(computed - analytical);
        total_error += error;
        count++;
      }
    }

    ASSERT_GT(count, 0);
    Real avg_error = total_error / count;
    EXPECT_LT(avg_error, 0.1)
      << "Average error should be small for boundary condition";

    // Verify boundary values
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0];
      Real val = u[it->getIndex()];

      if (x < boundary_tolerance)
      {
        EXPECT_LT(val, TIGHT_TOLERANCE)
          << "Boundary should have near-zero distance";
      }
    }
  }

  // Test 10: Interior interface with known analytical solution
  TEST_F(FMMManufacturedTest, InteriorInterface_AnalyticalValidation)
  {
    // Create 2D mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 25, 25 });
    mesh.scale(2.0 / 24.0);  // Scale to [0, 2] x [0, 2]
    mesh.displace(VectorFunction{ -1.0, -1.0 });  // Center at origin: [-1, 1] x [-1, 1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set square interface in interior: [-0.2, 0.2] x [-0.2, 0.2]
    std::vector<Index> interface;
    const Real square_size = 0.2;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0], y = coord[1];

      // Check if on boundary of square
      bool on_square_boundary = false;
      if ((std::abs(x - square_size) < 0.04 && std::abs(y) <= square_size + 0.02) ||
          (std::abs(x + square_size) < 0.04 && std::abs(y) <= square_size + 0.02) ||
          (std::abs(y - square_size) < 0.04 && std::abs(x) <= square_size + 0.02) ||
          (std::abs(y + square_size) < 0.04 && std::abs(x) <= square_size + 0.02))
      {
        on_square_boundary = true;
      }

      if (on_square_boundary)
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Interior square interface should not be empty";

    fmm.seed(interface).solve();

    // Analytical solution: distance to the nearest point on the square boundary
    Real total_error = 0.0;
    Index count = 0;

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0], y = coord[1];
      Real computed = u[it->getIndex()];

      // Test points outside the square and not too close to domain boundary
      if ((std::abs(x) > square_size + 0.1 || std::abs(y) > square_size + 0.1) &&
          std::abs(x) < 0.8 && std::abs(y) < 0.8)
      {
        // Analytical distance to square boundary
        Real dist_to_x_sides = std::max(0.0, std::abs(x) - square_size);
        Real dist_to_y_sides = std::max(0.0, std::abs(y) - square_size);
        Real analytical = std::sqrt(dist_to_x_sides * dist_to_x_sides + 
                                  dist_to_y_sides * dist_to_y_sides);

        if (std::isfinite(computed) && analytical > 0.05)
        {
          Real error = std::abs(computed - analytical);
          total_error += error;
          count++;
        }
      }
    }

    if (count > 0)
    {
      Real avg_error = total_error / count;
      EXPECT_LT(avg_error, 0.15)
        << "Average error should be reasonable for interior interface";
    }

    // Verify general properties
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      if (std::isfinite(val))
      {
        EXPECT_GE(val, 0.0) << "Distance should be non-negative";
      }
    }
  }

  // Test 11: Mixed speed function with validation
  TEST_F(FMMManufacturedTest, MixedSpeedFunction_LayeredMedium)
  {
    // Create 2D mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 25, 25 });
    mesh.scale(2.0 / 24.0);  // Scale to [0, 2] x [0, 2]
    mesh.displace(VectorFunction{ -1.0, -1.0 });  // Center at origin: [-1, 1] x [-1, 1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Layered speed function: fast on left (x < 0), slow on right (x > 0)
    auto speed = [](const Geometry::Point& p) -> Real
    {
      const auto& coord = p.getCoordinates();
      Real x = coord[0];
      return (x < 0.0) ? 2.0 : 0.5;  // Fast left, slow right
    };

    Models::Eikonal::FMM fmm(u, speed);

    // Set source at left boundary
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0];
      if (std::abs(x + 1.0) < 0.05)  // Left boundary x = -1
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Left boundary interface should not be empty";

    fmm.seed(interface).solve();

    // Verify that travel times are different in different layers
    Real left_travel_time = std::numeric_limits<Real>::infinity();
    Real right_travel_time = std::numeric_limits<Real>::infinity();

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0], y = coord[1];
      Real val = u[it->getIndex()];

      if (std::isfinite(val) && std::abs(y) < 0.2)  // Points near y = 0
      {
        // Left side (fast medium)
        if (std::abs(x + 0.5) < 0.1)
          left_travel_time = std::min(left_travel_time, val);

        // Right side (slow medium)
        if (std::abs(x - 0.5) < 0.1)
          right_travel_time = std::min(right_travel_time, val);
      }
    }

    // In the fast medium (left), travel time should be: distance / speed = 0.5 / 2.0 = 0.25
    // In the slow medium (right), travel time should be: (interface_to_x0 / fast_speed) + (x0_to_point / slow_speed)
    // = 1.0 / 2.0 + 0.5 / 0.5 = 0.5 + 1.0 = 1.5
    if (std::isfinite(left_travel_time) && std::isfinite(right_travel_time))
    {
      EXPECT_LT(left_travel_time, 0.4)
        << "Travel time in fast medium should be small";
      EXPECT_GT(right_travel_time, 1.0)
        << "Travel time in slow medium should be larger";
      EXPECT_GT(right_travel_time, left_travel_time)
        << "Slow medium should have longer travel times";
    }

    // Verify general properties
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      if (std::isfinite(val))
      {
        EXPECT_GE(val, 0.0) << "Distance should be non-negative";
      }
    }
  }
}
