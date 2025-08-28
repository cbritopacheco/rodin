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

    ASSERT_FALSE(interface.empty()) << "Interface should not be empty";

    fmm.setFront(std::move(interface));
    fmm.solve();

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

    EXPECT_LT(center_val, corner_val) << "Distance should increase away from source";
    EXPECT_GE(center_val, 0.0) << "Distance should be non-negative";
    EXPECT_LT(center_val, TOLERANCE) << "Source should have near-zero distance";
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

    EXPECT_EQ(mesh.getDimension(), 2) << "Surface mesh should be 2D";
    EXPECT_EQ(mesh.getSpaceDimension(), 3) << "Surface mesh should be embedded in 3D";

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

    ASSERT_FALSE(interface.empty()) << "Interface should not be empty";

    fmm.setFront(std::move(interface));
    fmm.solve();

    // Verify solution properties
    Real min_dist = std::numeric_limits<Real>::infinity();
    Real max_dist = 0.0;

    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val)) << "Solution should not contain NaN";
      EXPECT_FALSE(std::isinf(val)) << "Solution should be finite";
      EXPECT_GE(val, 0.0) << "Distance should be non-negative";

      min_dist = std::min(min_dist, val);
      max_dist = std::max(max_dist, val);
    }

    EXPECT_LT(min_dist, TOLERANCE) << "Minimum distance should be near zero at source";
    EXPECT_GT(max_dist, 0.5) << "Maximum distance should be reasonable for unit box";
  }

  // Test 3: Volumetric mesh - 3D tetrahedra
  TEST_F(FMMTest, VolumeMesh_3D_Tetrahedron)
  {
    // Create 3D tetrahedral mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, { 8, 8, 8 });
    mesh.scale(1.0 / 7.0);  // Scale to [0,1]^3
    mesh.getConnectivity().compute(3, 0);
    mesh.getConnectivity().compute(0, 0);

    EXPECT_EQ(mesh.getDimension(), 3) << "Volume mesh should be 3D";
    EXPECT_EQ(mesh.getSpaceDimension(), 3) << "Volume mesh should be in 3D space";

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

    ASSERT_FALSE(interface.empty()) << "Interface should not be empty";

    fmm.setFront(std::move(interface));
    fmm.solve();

    // Verify solution properties
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val)) << "Solution should not contain NaN";
      EXPECT_FALSE(std::isinf(val)) << "Solution should be finite";
      EXPECT_GE(val, 0.0) << "Distance should be non-negative";
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

    ASSERT_FALSE(interface.empty()) << "Interface should not be empty";

    fmm.setFront(std::move(interface));
    fmm.solve();

    // Verify solution is reasonable
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val)) << "Solution should not contain NaN";
      EXPECT_FALSE(std::isinf(val)) << "Solution should be finite";
      EXPECT_GE(val, 0.0) << "Distance should be non-negative";
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

    ASSERT_FALSE(interface.empty()) << "Interface should not be empty";

    fmm.setFront(std::move(interface));
    fmm.solve();

    // Verify solution
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val)) << "Solution should not contain NaN";
      EXPECT_FALSE(std::isinf(val)) << "Solution should be finite";
      EXPECT_GE(val, 0.0) << "Distance should be non-negative";
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
      EXPECT_LT(min_nearby_dist, TOLERANCE) << "Points near sources should have small distances";
    }
  }

  // Test 6: Error handling - invalid speed function
  TEST_F(FMMTest, InvalidSpeedFunction)
  {
    // Create simple mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 8, 8 });
    mesh.scale(1.0 / 7.0);
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    // Invalid speed function (negative speed)
    auto invalid_speed = [](const Geometry::Point& p) -> Real { return -1.0; };

    Models::Eikonal::FMM fmm(u, invalid_speed);

    // Set source
    std::vector<Index> interface = {0};
    fmm.setFront(std::move(interface));

    // Will cause assertion failure or exception
    EXPECT_EXIT(fmm.solve();, ::testing::KilledBySignal(SIGABRT), ".*");
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
    fmm.setFront(std::move(interface));
    fmm.solve();

    // All values should remain infinite
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      Real val = u[it->getIndex()];
      EXPECT_TRUE(std::isinf(val)) << "Without sources, all distances should be infinite";
    }
  }
}
