/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Models/Eikonal/FMM.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Basic unit tests for Fast Marching Method (FMM)
   */
  class FMMTest : public ::testing::Test
  {
    protected:
      static constexpr Real TOLERANCE = 1e-3;
  };

  // Test 1: Basic functionality - single point source on 2D triangular mesh
  TEST_F(FMMTest, SinglePointSource_2D_Triangle)
  {
    // Create simple 2D triangular mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    mesh.scale(1.0 / 15.0);  // Scale to [0,1] x [0,1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed function
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set source at center
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance = (coord - Math::SpatialVector<Real>{{0.5, 0.5}}).norm();
      if (distance < 0.05)
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    fmm.seed(interface).solve();

    // Test that solution is reasonable
    // Distance should increase away from center
    Real center_val = std::numeric_limits<Real>::infinity();
    Real corner_val = std::numeric_limits<Real>::infinity();

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real dist_to_center = (coord - Math::SpatialVector<Real>{{0.5, 0.5}}).norm();
      Real dist_to_corner = (coord - Math::SpatialVector<Real>{{0.0, 0.0}}).norm();

      if (dist_to_center < 0.1)
        center_val = std::min(center_val, u[it->getIndex()]);
      if (dist_to_corner < 0.1)
        corner_val = std::min(corner_val, u[it->getIndex()]);
    }

    EXPECT_LT(center_val, corner_val)
      << "Distance should increase away from source";

    EXPECT_GE(center_val, 0.0)
      << "Distance should be non-negative";

    EXPECT_LT(center_val, TOLERANCE)
      << "Source should have near-zero distance";
  }

  // Test 2: Surface mesh - 2D surface embedded in 3D
  TEST_F(FMMTest, SurfaceMesh_Box_Triangle)
  {
    // Create surface mesh using Box function
    Mesh mesh;
    mesh = mesh.Box(Polytope::Type::Triangle, { 8, 8, 8 });
    mesh.scale(1.0 / 7.0);  // Scale to [0,1]^3

    // Extract skin (surface)
    mesh.getConnectivity().compute(2, 0); // face to vertex

    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    EXPECT_EQ(mesh.getDimension(), 2)
      << "Surface mesh should be 2D";

    EXPECT_EQ(mesh.getSpaceDimension(), 3)
      << "Surface mesh should be embedded in 3D";

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed function
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set source at one corner of the box
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance = (coord - Math::SpatialVector<Real>{{0.0, 0.0, 0.0}}).norm();
      if (distance < 0.1)
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    fmm.seed(interface).solve();

    // Verify solution properties
    Real min_dist = std::numeric_limits<Real>::infinity();
    Real max_dist = 0.0;

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val))
        << "Solution should not contain NaN";
      EXPECT_FALSE(std::isinf(val))
        << "Solution should be finite";
      EXPECT_GE(val, 0.0)
        << "Distance should be non-negative";

      min_dist = std::min(min_dist, val);
      max_dist = std::max(max_dist, val);
    }

    EXPECT_LT(min_dist, TOLERANCE)
      << "Minimum distance should be near zero at source";
    EXPECT_GT(max_dist, 0.5)
      << "Maximum distance should be reasonable for unit box";
  }

  // Test 3: Volumetric mesh - 3D tetrahedra
  TEST_F(FMMTest, VolumeMesh_3D_Tetrahedron)
  {
    // Create 3D tetrahedral mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, { 16, 16, 16 });
    mesh.scale(1.0 / 15.0);  // Scale to [0,1]^3
    mesh.getConnectivity().compute(3, 0);
    mesh.getConnectivity().compute(0, 0);

    EXPECT_EQ(mesh.getDimension(), 3)
      << "Volume mesh should be 3D";
    EXPECT_EQ(mesh.getSpaceDimension(), 3)
      << "Volume mesh should be in 3D space";

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed function
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set source at center of cube
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto& coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance = (coord - Math::SpatialVector<Real>{{0.5, 0.5, 0.5}}).norm();
      if (distance < 0.1)
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    fmm.seed(interface).solve();

    // Verify solution properties
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val))
        << "Solution should not contain NaN";
      EXPECT_FALSE(std::isinf(val))
        << "Solution should be finite";
      EXPECT_GE(val, 0.0)
        << "Distance should be non-negative";
    }

    // Check monotonicity: distance should increase away from center
    Real center_distance = std::numeric_limits<Real>::infinity();
    Real corner_distance = 0.0;

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real geometric_dist_to_center = (coord - Math::SpatialVector<Real>{{0.5, 0.5, 0.5}}).norm();
      Real geometric_dist_to_corner = (coord - Math::SpatialVector<Real>{{0.0, 0.0, 0.0}}).norm();

      if (geometric_dist_to_center < 0.1)
        center_distance = std::min(center_distance, u[it->getIndex()]);
      if (geometric_dist_to_corner < 0.1)
        corner_distance = std::max(corner_distance, u[it->getIndex()]);
    }

    EXPECT_LT(center_distance, TOLERANCE) << "Center should have near-zero distance";
    EXPECT_GT(corner_distance, 0.5) << "Corner should have larger distance";
  }

  // Test 4: Variable speed function
  TEST_F(FMMTest, VariableSpeedFunction_2D)
  {
    // Create 2D mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    mesh.scale(1.0 / 15.0);
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Variable speed function - slower in center, faster at edges
    auto speed = [](const Geometry::Point& p) -> Real 
    {
      Real r = (p.getCoordinates() - Math::SpatialVector<Real>{{0.5, 0.5}}).norm();
      return 0.1 + 2.0 * r;  // Speed increases with distance from center
    };

    Models::Eikonal::FMM fmm(u, speed);

    // Set source at center
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance = (coord - Math::SpatialVector<Real>{{0.5, 0.5}}).norm();
      if (distance < 0.05)
        interface.push_back(it->getIndex());
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    fmm.seed(interface).solve();

    // Verify solution is reasonable
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val))
        << "Solution should not contain NaN";
      EXPECT_FALSE(std::isinf(val))
        << "Solution should be finite";
      EXPECT_GE(val, 0.0)
        << "Distance should be non-negative";
    }
  }

  // Test 5: Multiple sources
  TEST_F(FMMTest, MultipleSources_2D)
  {
    // Create 2D mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    mesh.scale(1.0 / 15.0);
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set multiple sources at corners
    std::vector<Index> interface;
    std::vector<Math::SpatialVector<Real>> source_locations =
    {
      Math::SpatialVector<Real>{{0.2, 0.2}},
      Math::SpatialVector<Real>{{0.8, 0.8}}
    };

    for (const auto& source : source_locations)
    {
      for (auto it = mesh.getVertex(); !it.end(); ++it)
      {
        const auto coord = mesh.getVertexCoordinates(it->getIndex());
        Real distance = (coord - source).norm();
        if (distance < 0.1)
          interface.push_back(it->getIndex());
      }
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    fmm.seed(interface).solve();

    // Verify solution
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val))
        << "Solution should not contain NaN";
      EXPECT_FALSE(std::isinf(val))
        << "Solution should be finite";
      EXPECT_GE(val, 0.0)
        << "Distance should be non-negative";
    }

    // Check that points near sources have small distances
    for (const auto& source : source_locations)
    {
      Real min_nearby_dist = std::numeric_limits<Real>::infinity();
      for (auto it = mesh.getVertex(); !it.end(); ++it)
      {
        const auto coord = mesh.getVertexCoordinates(it->getIndex());
        Real geometric_dist = (coord - source).norm();
        if (geometric_dist < 0.15)
          min_nearby_dist = std::min(min_nearby_dist, u[it->getIndex()]);
      }
      EXPECT_LT(min_nearby_dist, TOLERANCE)
        << "Points near sources should have small distances";
    }
  }

  // Test 7: Empty interface handling
  TEST_F(FMMTest, EmptyInterface)
  {
    // Create simple mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 8, 8 });
    mesh.scale(1.0 / 7.0);
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Empty interface
    std::vector<Index> interface;
    fmm.seed(interface).solve();

    // All values should remain infinite
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_TRUE(std::isinf(val))
        << "Without sources, all distances should be infinite";
    }
  }

  // Test 8: Boundary interface - initial front is the boundary of the mesh
  TEST_F(FMMTest, BoundaryInterface_2D)
  {
    // Create 2D triangular mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    mesh.scale(1.0 / 15.0);  // Scale to [0,1] x [0,1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed function
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Set interface to boundary vertices (edges of the domain)
    std::vector<Index> interface;
    const Real boundary_tolerance = 1e-10;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0], y = coord[1];
      // Check if vertex is on boundary (x=0, x=1, y=0, or y=1)
      if (x < boundary_tolerance || x > 1.0 - boundary_tolerance ||
          y < boundary_tolerance || y > 1.0 - boundary_tolerance)
      {
        interface.push_back(it->getIndex());
      }
    }

    ASSERT_FALSE(interface.empty())
      << "Boundary interface should not be empty";
    EXPECT_GT(interface.size(), 50)
      << "Should have many boundary vertices";

    fmm.seed(interface).solve();

    // Verify solution properties
    // Distance should be zero at boundary and increase toward interior
    Real center_dist = std::numeric_limits<Real>::infinity();
    Real boundary_dist = std::numeric_limits<Real>::infinity();

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      Real x = coord[0], y = coord[1];
      Real val = u[it->getIndex()];

      EXPECT_FALSE(std::isnan(val))
        << "Solution should not contain NaN";
      EXPECT_GE(val, 0.0)
        << "Distance should be non-negative";

      // Check boundary vertices
      if (x < boundary_tolerance || x > 1.0 - boundary_tolerance ||
          y < boundary_tolerance || y > 1.0 - boundary_tolerance)
      {
        boundary_dist = std::min(boundary_dist, val);
      }

      // Check center region
      if (std::abs(x - 0.5) < 0.1 && std::abs(y - 0.5) < 0.1)
      {
        center_dist = std::min(center_dist, val);
      }
    }

    EXPECT_LT(boundary_dist, TOLERANCE)
      << "Boundary should have near-zero distance";
    EXPECT_GT(center_dist, 0.2)
      << "Center should be farther from boundary";
  }

  // Test 9: Interior rectangular interface
  TEST_F(FMMTest, InteriorRectangularInterface_2D)
  {
    // Create 2D triangular mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 24, 24 });
    mesh.scale(1.0 / 23.0);  // Scale to [0,1] x [0,1]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    const Real h = 1.0 / 23.0;

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed function
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      const Real x = coord[0], y = coord[1];

      // Check if vertex is on the boundary of the rectangular region
      bool on_rect_boundary = false;

      if (std::abs(x - h * 4) < h / 2 && y >= h * 4 && y <= h * 20)
        on_rect_boundary = true;

      if (std::abs(x - h * 20) < h / 2 && y >= h * 4 && y <= h * 20)
        on_rect_boundary = true;

      if (std::abs(y - h * 4) < h / 2 && x >= h * 4 && x <= h * 20)
        on_rect_boundary = true;

      if (std::abs(y - h * 20) < h / 2 && x >= h * 4 && x <= h * 20)
        on_rect_boundary = true;

      if (on_rect_boundary)
      {
        interface.push_back(it->getIndex());
        mesh.setAttribute({ 0, it->getIndex() }, 666); // Reset attribute
      }
    }

    ASSERT_FALSE(interface.empty())
      << "Interior rectangular interface should not be empty";

    EXPECT_GT(interface.size(), 10)
      << "Should have multiple interface vertices";

    fmm.seed(interface).solve();

    // Seeds must be ~0
    for (Index id : interface)
      EXPECT_LT(std::abs(u[id]), 1e-14);

    // Rectangle [x1,x2]×[y1,y2]
    const Real x1 = 4*h, x2 = 20*h, y1 = 4*h, y2 = 20*h;

    // Analytic distance to the rectangle boundary
    auto rect_dist = [&](Real x, Real y) -> Real {
      const bool inside = (x >= x1 && x <= x2 && y >= y1 && y <= y2);
      if (inside)
        return std::min(std::min(x - x1, x2 - x), std::min(y - y1, y2 - y));
      const Real dx = (x < x1 ? x1 - x : (x > x2 ? x - x2 : 0.0));
      const Real dy = (y < y1 ? y1 - y : (y > y2 ? y - y2 : 0.0));
      return (dx > 0.0 && dy > 0.0) ? std::sqrt(dx*dx + dy*dy) : (dx + dy);
    };

    const Real tol = 3.0*h;

    Real max_err = 0.0;
    Real center_u = std::numeric_limits<Real>::quiet_NaN();
    Real corner00_u = std::numeric_limits<Real>::quiet_NaN();

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      const Real x = coord[0], y = coord[1];
      const Real val = u[it->getIndex()];
      const Real d   = rect_dist(x, y);

      EXPECT_FALSE(std::isnan(val));
      EXPECT_GE(val, 0.0);

      max_err = std::max(max_err, std::abs(val - d));

      if (std::abs(x - 0.5) < h/2 && std::abs(y - 0.5) < h/2)
        center_u = val;

      if (x < h/2 && y < h/2) // picks (0,0)
        corner00_u = val;
    }

    // Global accuracy
    EXPECT_LT(max_err, tol) << "Max nodal error vs analytic distance too large";

    // Center check: distance to nearest side = min(0.5 - x1, x2 - 0.5, 0.5 - y1, y2 - 0.5)
    const Real center_d = std::min(std::min(0.5 - x1, x2 - 0.5), std::min(0.5 - y1, y2 - 0.5));
    ASSERT_TRUE(std::isfinite(center_u));
    EXPECT_NEAR(center_u, center_d, tol);

    // Corner check at (0,0): distance to rectangle = sqrt((x1)^2 + (y1)^2)
    const Real corner00_d = std::sqrt(x1*x1 + y1*y1);
    ASSERT_TRUE(std::isfinite(corner00_u));
    EXPECT_NEAR(corner00_u, corner00_d, tol);
  }

  // Test 10: Large mesh stress test
  TEST_F(FMMTest, LargeMeshStressTest_2D)
  {
    // Create larger 2D mesh to test performance
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 32, 32 });
    mesh.scale(2.0 / 31.0);  // Scale to [0,2] x [0,2]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Variable speed function
    auto speed = [](const Geometry::Point& p) -> Real 
    {
      const auto& coord = p.getCoordinates();
      Real x = coord[0], y = coord[1];
      return 0.5 + 0.5 * std::sin(x) * std::cos(y);
    };

    Models::Eikonal::FMM fmm(u, speed);

    // Set multiple sources scattered across the domain
    std::vector<Index> interface;
    std::vector<Math::SpatialVector<Real>> source_locations = {
      Math::SpatialVector<Real>{{0.5, 0.5}},
      Math::SpatialVector<Real>{{1.5, 0.5}},
      Math::SpatialVector<Real>{{0.5, 1.5}},
      Math::SpatialVector<Real>{{1.5, 1.5}}
    };

    for (const auto& source : source_locations)
    {
      for (auto it = mesh.getVertex(); !it.end(); ++it)
      {
        const auto coord = mesh.getVertexCoordinates(it->getIndex());
        Real distance = (coord - source).norm();
        if (distance < 0.1)
          interface.push_back(it->getIndex());
      }
    }

    ASSERT_FALSE(interface.empty())
      << "Interface should not be empty";

    fmm.seed(interface).solve();

    // Verify solution is reasonable on large mesh
    Index valid_vertices = 0;
    Real max_dist = 0.0;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val))
        << "Solution should not contain NaN";

      if (std::isfinite(val))
      {
        EXPECT_GE(val, 0.0)
          << "Distance should be non-negative";
        max_dist = std::max(max_dist, val);
        valid_vertices++;
      }
    }

    EXPECT_GT(valid_vertices, 800)
      << "Most vertices should have valid solution";
    EXPECT_LT(max_dist, 10.0)
      << "Maximum distance should be reasonable";
  }

  // Test 11: Multiple disconnected interfaces
  TEST_F(FMMTest, MultipleDisconnectedInterfaces_2D)
  {
    // Create 2D mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 32, 32 });
    mesh.scale(2.0 / 31.0);  // Scale to [0,2] x [0,2]
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Constant speed
    auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

    Models::Eikonal::FMM fmm(u, speed);

    // Create multiple small disconnected circular interfaces
    std::vector<Index> interface;
    std::vector<std::pair<Math::SpatialVector<Real>, Real>> circles = {
      {Math::SpatialVector<Real>{{0.5, 0.5}}, 0.1},
      {Math::SpatialVector<Real>{{1.5, 0.5}}, 0.1},
      {Math::SpatialVector<Real>{{0.5, 1.5}}, 0.1},
      {Math::SpatialVector<Real>{{1.5, 1.5}}, 0.1}
    };

    for (const auto& circle : circles)
    {
      const auto& center = circle.first;
      Real radius = circle.second;

      for (auto it = mesh.getVertex(); !it.end(); ++it)
      {
        const auto coord = mesh.getVertexCoordinates(it->getIndex());
        Real distance = (coord - center).norm();
        if (std::abs(distance - radius) < 0.05)  // Points on circle boundary
          interface.push_back(it->getIndex());
      }
    }

    ASSERT_FALSE(interface.empty())
      << "Multiple interfaces should not be empty";

    fmm.seed(interface).solve();

    // Verify solution properties
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val))
        << "Solution should not contain NaN";

      if (std::isfinite(val))
      {
        EXPECT_GE(val, 0.0)
          << "Distance should be non-negative";
      }
    }

    // Check that points inside circles have small distances
    for (const auto& circle : circles)
    {
      const auto& center = circle.first;
      Real radius = circle.second;
      Real min_interior_dist = std::numeric_limits<Real>::infinity();

      for (auto it = mesh.getVertex(); !it.end(); ++it)
      {
        const auto coord = mesh.getVertexCoordinates(it->getIndex());
        Real distance = (coord - center).norm();
        if (distance < radius - 0.02)  // Interior points
        {
          min_interior_dist = std::min(min_interior_dist, u[it->getIndex()]);
        }
      }

      if (std::isfinite(min_interior_dist))
      {
        EXPECT_LT(min_interior_dist, 0.1) 
          << "Interior of circular interfaces should have small distances";
      }
    }
  }
}
