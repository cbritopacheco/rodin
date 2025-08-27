/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>

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

    ASSERT_FALSE(interface.empty()) << "Interface should not be empty";

    Models::Eikonal::FMM fmm(u, speed);
    fmm.setInterface(std::move(interface));
    fmm.solve();

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

    ASSERT_FALSE(interface.empty()) << "Interface should not be empty";

    Models::Eikonal::FMM fmm(u, speed);
    fmm.setInterface(std::move(interface));
    fmm.solve();

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

    EXPECT_LT(center_val, TIGHT_TOLERANCE) << "Center should have near-zero distance";
    EXPECT_GT(edge_val, 0.2) << "Edge should have reasonable distance";
    EXPECT_LT(center_val, edge_val) << "Distance should increase with radius";
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


    EXPECT_EQ(mesh.getDimension(), 2) << "Surface mesh should be 2D";
    EXPECT_EQ(mesh.getSpaceDimension(), 3) << "Surface mesh should be embedded in 3D";

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
      {auto const& v2c = mesh.getConnectivity().getIncidence(0,2).at(it->getIndex());
std::cout << "cells at corner: " << v2c.size() << "\n";
        interface.push_back(it->getIndex());
      }
    }
size_t n = 0;
std::vector<char> seen(mesh.getVertexCount(),0);
for (auto it = mesh.getVertex(); !it.end(); ++it) if (!seen[it->getIndex()]) {
  ++n;
  std::queue<Index> q; q.push(it->getIndex()); seen[it->getIndex()] = 1;
  while (!q.empty()) {
    auto v = q.front(); q.pop();
    for (auto nb : mesh.getConnectivity().getIncidence(0,0).at(v))
      if (!seen[nb]) { seen[nb]=1; q.push(nb); }
  }
}
std::cout << "Connected components (0-0): " << n << "\n";

    std::cout << "Number of source points: " << interface.size() << std::endl;

    ASSERT_FALSE(interface.empty()) << "Interface should not be empty";

    fmm.setInterface(std::move(interface));
    fmm.solve();

    u.save("fmm_box_surface_point_source.gf");
    mesh.save("fmm_box_surface_point_source.mesh");

    // Verify solution properties on surface
    Real min_dist = std::numeric_limits<Real>::infinity();
    Real max_dist = 0.0;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const Real val = u[it->getIndex()];
      EXPECT_FALSE(std::isnan(val)) << "Solution should not contain NaN";
      EXPECT_GE(val, 0.0) << "Distance should be non-negative";
      if (std::isfinite(val))
      {
        min_dist = std::min(min_dist, val);
        max_dist = std::max(max_dist, val);
      }
    }
    EXPECT_LT(min_dist, TIGHT_TOLERANCE) << "Minimum distance should be near zero at source";
    EXPECT_GT(max_dist, 1.0) << "Maximum distance should span reasonable range on surface";
  }

  // Test 4: Sphere mapping test - Map cube to sphere and test geodesic distances
  // TEST_F(FMMManufacturedTest, SphereMappingTest_GeodesicDistance)
  // {
  //   // Create cube mesh
  //   Mesh surface_mesh;
  //   surface_mesh = surface_mesh.Box(Polytope::Type::Triangle, { 12, 12, 12 });
  //   surface_mesh.scale(2.0 / 11.0);  // Scale to [0, 2]^3
  //   surface_mesh.displace(VectorFunction{ -1.0, -1.0, -1.0 });  // Center at origin: [-1, 1]^3

  //   // Extract surface
  //   surface_mesh.getConnectivity().compute(2, 0);

  //   // Map surface vertices to sphere
  //   Math::SpatialVector<Real> coord;
  //   for (auto it = surface_mesh.getVertex(); !it.end(); ++it)
  //   {
  //     coord = surface_mesh.getVertexCoordinates(it->getIndex());
  //     Real norm = coord.norm();
  //     if (norm > 1e-12)  // Avoid division by zero
  //       coord = coord / norm;  // Project to unit sphere
  //     surface_mesh.setVertexCoordinates(it->getIndex(), coord);
  //   }
  //   surface_mesh.getConnectivity().compute(2, 0);
  //   surface_mesh.getConnectivity().compute(0, 0);

  //   // Take midsection of sphere (z ≈ 0)
  //   Mesh<Context::Local>::Builder mesh_builder;
  //   mesh_builder.initialize(3);  // 3D embedding
  //   std::vector<Index> vertex_map;
  //   std::vector<Index> valid_vertices;

  //   // Add vertices in midsection (|z| < 0.3)
  //   Index new_idx = 0;
  //   for (auto it = surface_mesh.getVertex(); !it.end(); ++it)
  //   {
  //     const auto coord = surface_mesh.getVertexCoordinates(it->getIndex());
  //     if (std::abs(coord[2]) < 0.3)  // Midsection condition
  //     {
  //       mesh_builder.vertex(coord);
  //       vertex_map.push_back(new_idx);
  //       valid_vertices.push_back(it->getIndex());
  //       new_idx++;
  //     }
  //     else
  //     {
  //       vertex_map.push_back(Index(-1));  // Invalid marker
  //     }
  //   }

  //   // Add triangles where all vertices are in midsection
  //   for (auto it = surface_mesh.getPolytope(2); !it.end(); ++it)
  //   {
  //     const auto& verts = it->getVertices();
  //     if (verts.size() == 3 &&
  //         vertex_map[verts[0]] != Index(-1) &&
  //         vertex_map[verts[1]] != Index(-1) &&
  //         vertex_map[verts[2]] != Index(-1))
  //     {
  //       mesh_builder.polytope(Polytope::Type::Triangle,
  //         {vertex_map[verts[0]], vertex_map[verts[1]], vertex_map[verts[2]]});
  //     }
  //   }

  //   Mesh mesh = mesh_builder.finalize();
  //   mesh.getConnectivity().compute(2, 0);
  //   mesh.getConnectivity().compute(0, 0);

  //   if (mesh.getVertexCount() < 10)
  //   {
  //     GTEST_SKIP() << "Insufficient vertices in sphere midsection for meaningful test";
  //     return;
  //   }

  //   P1 vh(mesh);
  //   GridFunction u(vh);

  //   // Constant speed on sphere surface
  //   auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

  //   Models::Eikonal::FMM fmm(u, speed);

  //   // Set source at one point on the sphere (e.g., (1, 0, 0))
  //   std::vector<Index> interface;
  //   for (auto it = mesh.getVertex(); !it.end(); ++it)
  //   {
  //     const auto coord = mesh.getVertexCoordinates(it->getIndex());
  //     Real distance = (coord - Math::SpatialVector<Real>{{1.0, 0.0, 0.0}}).norm();
  //     if (distance < 0.2)
  //       interface.push_back(it->getIndex());
  //   }

  //   if (interface.empty())
  //   {
  //     GTEST_SKIP() << "No suitable source point found on sphere midsection";
  //     return;
  //   }

  //   fmm.setInterface(std::move(interface));
  //   fmm.solve();

  //   // For a unit sphere, the geodesic distance between two points is the arc length
  //   // For points (1,0,0) and (-1,0,0), the geodesic distance is π
  //   Real max_expected_distance = M_PI;
  //   Real min_dist = std::numeric_limits<Real>::infinity();
  //   Real max_dist = 0.0;
  //   Real opposite_point_distance = std::numeric_limits<Real>::infinity();
  //   for (auto it = mesh.getVertex(); !it.end(); ++it)
  //   {
  //     const auto coord = mesh.getVertexCoordinates(it->getIndex());
  //     Real val = u[it->getIndex()];
  //     if (std::isfinite(val))
  //     {
  //       min_dist = std::min(min_dist, val);
  //       max_dist = std::max(max_dist, val);
  //       // Check point approximately opposite to source (1,0,0) -> (-1,0,0)
  //       Real dist_to_opposite = (coord - Math::SpatialVector<Real>{{-1.0, 0.0, 0.0}}).norm();
  //       if (dist_to_opposite < 0.2)
  //         opposite_point_distance = std::min(opposite_point_distance, val);
  //     }
  //   }

  //   EXPECT_LT(min_dist, TIGHT_TOLERANCE) << "Minimum distance should be near zero at source";
  //   EXPECT_LT(max_dist, max_expected_distance + TOLERANCE) << "Maximum distance should not exceed sphere circumference";
  //   // Check that opposite point has distance close to π (if it exists in midsection)
  //   if (std::isfinite(opposite_point_distance))
  //   {
  //     Real expected_opposite_dist = M_PI;
  //     Real error = std::abs(opposite_point_distance - expected_opposite_dist);
  //     EXPECT_LT(error, TOLERANCE * 2) << "Opposite point should have geodesic distance ≈ π";
  //   }
  // }

  // // Test 5: 3D volumetric manufactured test - Constant speed in cube
  // TEST_F(FMMManufacturedTest, Volume3D_ConstantSpeed_CubeCenter)
  // {
  //   // Create 3D tetrahedral mesh
  //   Mesh mesh;
  //   mesh = mesh.Box(Polytope::Type::Tetrahedron, { 16, 16, 16 });
  //   mesh.scale(2.0 / 15.0);  // Scale to [0, 2]^3
  //   mesh.displace(VectorFunction{ -1.0, -1.0, -1.0 });  // Center at origin: [-1, 1]^3
  //   mesh.getConnectivity().compute(3, 0);
  //   mesh.getConnectivity().compute(0, 0);

  //   P1 vh(mesh);
  //   GridFunction u(vh);

  //   // Constant speed
  //   auto speed = [](const Geometry::Point& p) -> Real { return 1.0; };

  //   Models::Eikonal::FMM fmm(u, speed);

  //   // Set source at center (0, 0, 0)
  //   std::vector<Index> interface;
  //   for (auto it = mesh.getVertex(); !it.end(); ++it)
  //   {
  //     const auto coord = mesh.getVertexCoordinates(it->getIndex());
  //     Real distance = coord.norm();
  //     if (distance < 0.1)
  //       interface.push_back(it->getIndex());
  //   }

  //   ASSERT_FALSE(interface.empty()) << "Interface should not be empty";

  //   fmm.setInterface(std::move(interface));
  //   fmm.solve();

  //   // Verify against analytical solution: u(x,y,z) = sqrt(x^2 + y^2 + z^2)
  //   Real max_error = 0.0;
  //   Real total_error = 0.0;
  //   Index count = 0;

  //   for (auto it = mesh.getVertex(); !it.end(); ++it)
  //   {
  //     const auto coord = mesh.getVertexCoordinates(it->getIndex());
  //     Real computed = u[it->getIndex()];
  //     Real analytical = coord.norm();

  //     if (std::isfinite(computed) && std::isfinite(analytical) && analytical > 0.05)
  //     {
  //       Real error = std::abs(computed - analytical);
  //       max_error = std::max(max_error, error);
  //       total_error += error;
  //       count++;
  //     }
  //   }

  //   if (count > 0)
  //   {
  //     Real avg_error = total_error / count;
  //     EXPECT_LT(max_error, TOLERANCE) << "Maximum error should be within tolerance";
  //     EXPECT_LT(avg_error, TOLERANCE / 2) << "Average error should be within tolerance";
  //   }
  // }

  // // Test 6: Anisotropic speed function test
  // TEST_F(FMMManufacturedTest, AnisotropicSpeed_2D)
  // {
  //   // Create 2D mesh
  //   Mesh mesh;
  //   mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 20, 20 });
  //   mesh.scale(2.0 / 19.0);  // Scale to [0, 2] x [0, 2]
  //   mesh.displace(VectorFunction{ -1.0, -1.0 });  // Center at origin: [-1, 1] x [-1, 1]
  //   mesh.getConnectivity().compute(2, 0);
  //   mesh.getConnectivity().compute(0, 0);

  //   P1 vh(mesh);
  //   GridFunction u(vh);

  //   // Anisotropic speed: faster in x-direction, slower in y-direction
  //   auto speed = [](const Geometry::Point& p) -> Real 
  //   {
  //     const auto& coord = p.getCoordinates();
  //     Real x = coord[0], y = coord[1];
  //     return 1.0 + 0.5 * std::abs(x) + 2.0 * std::abs(y);  // Varies with position
  //   };

  //   Models::Eikonal::FMM fmm(u, speed);

  //   // Set source at origin
  //   std::vector<Index> interface;
  //   for (auto it = mesh.getVertex(); !it.end(); ++it)
  //   {
  //     const auto coord = mesh.getVertexCoordinates(it->getIndex());
  //     Real distance = coord.norm();
  //     if (distance < 0.05)
  //       interface.push_back(it->getIndex());
  //   }

  //   ASSERT_FALSE(interface.empty()) << "Interface should not be empty";

  //   fmm.setInterface(std::move(interface));
  //   fmm.solve();

  //   // Verify solution properties
  //   Real x_positive_dist = std::numeric_limits<Real>::infinity();
  //   Real y_positive_dist = std::numeric_limits<Real>::infinity();

  //   for (auto it = mesh.getVertex(); !it.end(); ++it)
  //   {
  //     const auto coord = mesh.getVertexCoordinates(it->getIndex());
  //     Real val = u[it->getIndex()];

  //     // Check point at (0.5, 0)
  //     if (std::abs(coord[0] - 0.5) < 0.1 && std::abs(coord[1]) < 0.1)
  //       x_positive_dist = std::min(x_positive_dist, val);

  //     // Check point at (0, 0.5)  
  //     if (std::abs(coord[0]) < 0.1 && std::abs(coord[1] - 0.5) < 0.1)
  //       y_positive_dist = std::min(y_positive_dist, val);
  //   }

  //   // Due to anisotropic speed, travel time should be different in different directions
  //   EXPECT_GT(y_positive_dist, x_positive_dist) << "Travel in y-direction should take longer due to slower speed";
  // }
}
