/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/Variational/H1/GLL.h"
#include "Rodin/Variational/H1/Fekete.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  //==========================================================================
  // Common constants for face extraction
  //==========================================================================

  // Tolerance for geometric comparisons (whether a point is on a face)
  constexpr Real FACE_TOL = 1e-10;

  // Epsilon for division by zero protection
  const Real FACE_EPS = std::numeric_limits<Real>::epsilon() * 100;

  //==========================================================================
  // Helper functions for extracting face nodes
  //==========================================================================

  /**
   * @brief Checks if a 3D point lies on a triangular face and computes its barycentric coordinates.
   *
   * Uses the barycentric coordinate method to determine if point (px, py, pz) lies on the
   * triangle defined by vertices (x0,y0,z0), (x1,y1,z1), (x2,y2,z2).
   *
   * Algorithm:
   * 1. Express the point as P = u*V0 + v*(V1-V0) + w*(V2-V0) where u + v + w = 1
   * 2. Solve for v, w using dot products (Cramer's rule)
   * 3. Point is on triangle if 0 <= u,v,w <= 1 and point lies on the plane
   *
   * @param px, py, pz  The query point coordinates
   * @param x0, y0, z0  First vertex of the triangle (corresponds to barycentric u=1)
   * @param x1, y1, z1  Second vertex of the triangle (corresponds to barycentric v=1)
   * @param x2, y2, z2  Third vertex of the triangle (corresponds to barycentric w=1)
   * @return Tuple of (u, v, w) barycentric coordinates if point is on face;
   *         (NaN, NaN, NaN) if point is not on the triangular face
   */
  inline std::tuple<Real, Real, Real> computeTriangleFaceBarycentric(
      Real px, Real py, Real pz,
      Real x0, Real y0, Real z0,
      Real x1, Real y1, Real z1,
      Real x2, Real y2, Real z2)
  {
    // Compute edge vectors from vertex 0
    Real v0x = x1 - x0, v0y = y1 - y0, v0z = z1 - z0;
    Real v1x = x2 - x0, v1y = y2 - y0, v1z = z2 - z0;
    Real v2x = px - x0, v2y = py - y0, v2z = pz - z0;

    // Compute dot products
    Real d00 = v0x * v0x + v0y * v0y + v0z * v0z;
    Real d01 = v0x * v1x + v0y * v1y + v0z * v1z;
    Real d11 = v1x * v1x + v1y * v1y + v1z * v1z;
    Real d20 = v2x * v0x + v2y * v0y + v2z * v0z;
    Real d21 = v2x * v1x + v2y * v1y + v2z * v1z;

    // Compute barycentric coordinates
    Real denom = d00 * d11 - d01 * d01;
    if (std::abs(denom) < FACE_EPS)
    {
      return {std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN()};
    }

    Real v = (d11 * d20 - d01 * d21) / denom;
    Real w = (d00 * d21 - d01 * d20) / denom;
    Real u = 1.0 - v - w;

    // Check if point is in triangle (barycentric coords in [0,1])
    if (u >= -FACE_TOL && u <= 1.0 + FACE_TOL &&
        v >= -FACE_TOL && v <= 1.0 + FACE_TOL &&
        w >= -FACE_TOL && w <= 1.0 + FACE_TOL)
    {
      // Verify point is actually on the plane
      Real expected_x = x0 + v * v0x + w * v1x;
      Real expected_y = y0 + v * v0y + w * v1y;
      Real expected_z = z0 + v * v0z + w * v1z;

      if (std::abs(px - expected_x) < FACE_TOL &&
          std::abs(py - expected_y) < FACE_TOL &&
          std::abs(pz - expected_z) < FACE_TOL)
      {
        return {u, v, w};
      }
    }

    return {std::numeric_limits<Real>::quiet_NaN(),
            std::numeric_limits<Real>::quiet_NaN(),
            std::numeric_limits<Real>::quiet_NaN()};
  }

  // Axis indices for coordinate projection
  constexpr int AXIS_X = 0;
  constexpr int AXIS_Y = 1;
  constexpr int AXIS_Z = 2;

  /**
   * @brief Checks if a 3D point lies on a planar quadrilateral face and computes parametric coords.
   *
   * For a planar quadrilateral with vertices ordered counter-clockwise, determines if point
   * (px, py, pz) lies on the face and computes its parametric coordinates (s, t) in [0,1]x[0,1].
   *
   * Algorithm:
   * 1. Compute the face normal using cross product of edge vectors
   * 2. Check if point lies on the plane within tolerance
   * 3. Project to 2D by dropping the axis with largest normal component
   * 4. Solve for parametric coordinates in the projected space
   *
   * Note: This simplified algorithm assumes axis-aligned or near-axis-aligned quads,
   * which is appropriate for reference element faces.
   *
   * @param px, py, pz  The query point coordinates
   * @param x0, y0, z0  First vertex (maps to parametric (0,0))
   * @param x1, y1, z1  Second vertex (maps to parametric (1,0))
   * @param x2, y2, z2  Third vertex (maps to parametric (1,1))
   * @param x3, y3, z3  Fourth vertex (maps to parametric (0,1))
   * @return Pair of (s, t) parametric coordinates if point is on face;
   *         (NaN, NaN) if point is not on the quadrilateral face
   */
  inline std::pair<Real, Real> computeQuadFaceParametric(
      Real px, Real py, Real pz,
      Real x0, Real y0, Real z0,
      Real x1, Real y1, Real z1,
      Real x2, Real y2, Real z2,
      Real x3, Real y3, Real z3)
  {
    // Compute edge vectors from vertex 0
    Real v1x = x1 - x0, v1y = y1 - y0, v1z = z1 - z0;
    Real v2x = x3 - x0, v2y = y3 - y0, v2z = z3 - z0;

    // Compute normal vector via cross product
    Real nx = v1y * v2z - v1z * v2y;
    Real ny = v1z * v2x - v1x * v2z;
    Real nz = v1x * v2y - v1y * v2x;

    // Check if point lies on the plane (distance to plane < tolerance)
    Real d = nx * x0 + ny * y0 + nz * z0;
    Real dist = std::abs(nx * px + ny * py + nz * pz - d) /
                std::sqrt(nx * nx + ny * ny + nz * nz + FACE_EPS);

    if (dist > FACE_TOL)
    {
      return {std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN()};
    }

    // Project to 2D by dropping the coordinate with largest normal component
    // This minimizes numerical error in the 2D computation
    int drop_axis = AXIS_X;
    if (std::abs(ny) > std::abs(nx) && std::abs(ny) > std::abs(nz))
      drop_axis = AXIS_Y;
    else if (std::abs(nz) > std::abs(nx))
      drop_axis = AXIS_Z;

    auto project = [drop_axis](Real x, Real y, Real z) {
      if (drop_axis == AXIS_X)
        return std::make_pair(y, z);
      else if (drop_axis == AXIS_Y)
        return std::make_pair(x, z);
      else
        return std::make_pair(x, y);
    };

    auto [p2d_x, p2d_y] = project(px, py, pz);
    auto [q0_x, q0_y] = project(x0, y0, z0);
    auto [q1_x, q1_y] = project(x1, y1, z1);
    auto [q3_x, q3_y] = project(x3, y3, z3);

    // Solve for (s, t) using inverse bilinear mapping
    // Simplified for axis-aligned quads in reference space
    Real dx = q1_x - q0_x;
    Real dy = q3_y - q0_y;

    Real s, t;
    if (std::abs(dx) > FACE_EPS)
      s = (p2d_x - q0_x) / dx;
    else
      s = 0.0;

    if (std::abs(dy) > FACE_EPS)
      t = (p2d_y - q0_y) / dy;
    else
      t = 0.0;

    // Validate
    if (s >= -FACE_TOL && s <= 1.0 + FACE_TOL &&
        t >= -FACE_TOL && t <= 1.0 + FACE_TOL)
    {
      return {std::max(0.0, std::min(1.0, s)),
              std::max(0.0, std::min(1.0, t))};
    }

    return {std::numeric_limits<Real>::quiet_NaN(),
            std::numeric_limits<Real>::quiet_NaN()};
  }

  // Extract 2D coordinates of nodes lying on a specific triangular face of tetrahedron
  // Returns vector of (x, y) in reference triangle coordinates
  template <size_t K>
  std::vector<std::pair<Real, Real>> extractTetrahedronTriangleFaceNodes(
      Real x0, Real y0, Real z0,
      Real x1, Real y1, Real z1,
      Real x2, Real y2, Real z2)
  {
    std::vector<std::pair<Real, Real>> coords;
    const auto& nodes = FeketeTetrahedron<K>::getNodes();

    for (const auto& node : nodes)
    {
      auto [u, v, w] = computeTriangleFaceBarycentric(
          node.x(), node.y(), node.z(),
          x0, y0, z0, x1, y1, z1, x2, y2, z2);

      if (!std::isnan(u))
      {
        // Convert barycentric to 2D reference triangle coordinates
        // v corresponds to x, w corresponds to y (since u = 1 - v - w)
        coords.push_back({std::max(0.0, std::min(1.0, v)),
                          std::max(0.0, std::min(1.0, w))});
      }
    }

    return coords;
  }

  // Extract 2D coordinates of nodes lying on a specific triangular face of wedge
  template <size_t K>
  std::vector<std::pair<Real, Real>> extractWedgeTriangleFaceNodes(
      Real x0, Real y0, Real z0,
      Real x1, Real y1, Real z1,
      Real x2, Real y2, Real z2)
  {
    std::vector<std::pair<Real, Real>> coords;
    const auto& nodes = H1Element<K, Real>::getNodes(Polytope::Type::Wedge);

    for (const auto& node : nodes)
    {
      auto [u, v, w] = computeTriangleFaceBarycentric(
          node.x(), node.y(), node.z(),
          x0, y0, z0, x1, y1, z1, x2, y2, z2);

      if (!std::isnan(u))
      {
        coords.push_back({std::max(0.0, std::min(1.0, v)),
                          std::max(0.0, std::min(1.0, w))});
      }
    }

    return coords;
  }

  // Extract 2D coordinates from 2D Triangle element
  template <size_t K>
  std::vector<std::pair<Real, Real>> extractTriangleNodes()
  {
    std::vector<std::pair<Real, Real>> coords;
    const auto& nodes = FeketeTriangle<K>::getNodes();

    for (const auto& node : nodes)
    {
      coords.push_back({node.x(), node.y()});
    }

    return coords;
  }

  // Extract 2D coordinates of nodes lying on a specific quad face of wedge
  // Returns vector of (s, t) in reference quadrilateral coordinates [0,1]x[0,1]
  template <size_t K>
  std::vector<std::pair<Real, Real>> extractWedgeQuadFaceNodes(
      Real x0, Real y0, Real z0,  // (0,0)
      Real x1, Real y1, Real z1,  // (1,0)
      Real x2, Real y2, Real z2,  // (1,1)
      Real x3, Real y3, Real z3)  // (0,1)
  {
    std::vector<std::pair<Real, Real>> coords;
    const auto& nodes = H1Element<K, Real>::getNodes(Polytope::Type::Wedge);

    for (const auto& node : nodes)
    {
      auto [s, t] = computeQuadFaceParametric(
          node.x(), node.y(), node.z(),
          x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

      if (!std::isnan(s))
      {
        coords.push_back({s, t});
      }
    }

    return coords;
  }

  // Extract 2D coordinates from 2D Quadrilateral element
  template <size_t K>
  std::vector<std::pair<Real, Real>> extractQuadrilateralNodes()
  {
    std::vector<std::pair<Real, Real>> coords;
    const auto& nodes = H1Element<K, Real>::getNodes(Polytope::Type::Quadrilateral);

    for (const auto& node : nodes)
    {
      coords.push_back({node.x(), node.y()});
    }

    return coords;
  }

  // Sort 2D coordinates for comparison (by x then y)
  inline void sortCoords2D(std::vector<std::pair<Real, Real>>& coords)
  {
    std::sort(coords.begin(), coords.end(),
              [](const auto& a, const auto& b) {
                if (std::abs(a.first - b.first) > FACE_TOL)
                  return a.first < b.first;
                return a.second < b.second;
              });
  }

  // Compare two sets of 2D coordinates
  inline void compareFaceCoordinates(
      std::vector<std::pair<Real, Real>> face1,
      std::vector<std::pair<Real, Real>> face2,
      const std::string& name1,
      const std::string& name2)
  {
    sortCoords2D(face1);
    sortCoords2D(face2);

    ASSERT_EQ(face1.size(), face2.size())
        << "Face node counts differ: " << name1 << " has " << face1.size()
        << " nodes, " << name2 << " has " << face2.size() << " nodes";

    for (size_t i = 0; i < face1.size(); ++i)
    {
      EXPECT_NEAR(face1[i].first, face2[i].first, FACE_TOL)
          << "Node " << i << " x-coord differs: " << name1 << "[" << i << "].x=" << face1[i].first
          << ", " << name2 << "[" << i << "].x=" << face2[i].first;

      EXPECT_NEAR(face1[i].second, face2[i].second, FACE_TOL)
          << "Node " << i << " y-coord differs: " << name1 << "[" << i << "].y=" << face1[i].second
          << ", " << name2 << "[" << i << "].y=" << face2[i].second;
    }
  }

  //==========================================================================
  // Tetrahedron triangular faces match Triangle (2D)
  // Reference tetrahedron: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
  // Face 0 (base): vertices 0,1,2 => (0,0,0), (1,0,0), (0,1,0) on z=0 plane
  //==========================================================================

  TEST(FaceConformity, Tetrahedron_BaseFace_MatchesTriangle_K2)
  {
    // Base face of tetrahedron (z=0 plane) should match reference triangle
    auto tet_face = extractTetrahedronTriangleFaceNodes<2>(
        0.0, 0.0, 0.0,  // vertex 0
        1.0, 0.0, 0.0,  // vertex 1
        0.0, 1.0, 0.0); // vertex 2

    auto tri_nodes = extractTriangleNodes<2>();

    compareFaceCoordinates(tet_face, tri_nodes,
                           "Tetrahedron base face", "Triangle");
  }

  TEST(FaceConformity, Tetrahedron_BaseFace_MatchesTriangle_K3)
  {
    auto tet_face = extractTetrahedronTriangleFaceNodes<3>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto tri_nodes = extractTriangleNodes<3>();

    compareFaceCoordinates(tet_face, tri_nodes,
                           "Tetrahedron base face K=3", "Triangle K=3");
  }

  TEST(FaceConformity, Tetrahedron_BaseFace_MatchesTriangle_K4)
  {
    auto tet_face = extractTetrahedronTriangleFaceNodes<4>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto tri_nodes = extractTriangleNodes<4>();

    compareFaceCoordinates(tet_face, tri_nodes,
                           "Tetrahedron base face K=4", "Triangle K=4");
  }

  TEST(FaceConformity, Tetrahedron_BaseFace_MatchesTriangle_K5)
  {
    auto tet_face = extractTetrahedronTriangleFaceNodes<5>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto tri_nodes = extractTriangleNodes<5>();

    compareFaceCoordinates(tet_face, tri_nodes,
                           "Tetrahedron base face K=5", "Triangle K=5");
  }

  TEST(FaceConformity, Tetrahedron_BaseFace_MatchesTriangle_K6)
  {
    auto tet_face = extractTetrahedronTriangleFaceNodes<6>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto tri_nodes = extractTriangleNodes<6>();

    compareFaceCoordinates(tet_face, tri_nodes,
                           "Tetrahedron base face K=6", "Triangle K=6");
  }

  TEST(FaceConformity, Tetrahedron_Face013_MatchesTriangle_K3)
  {
    // Face with vertices 0,1,3: (0,0,0), (1,0,0), (0,0,1) on y=0 plane
    auto tet_face = extractTetrahedronTriangleFaceNodes<3>(
        0.0, 0.0, 0.0,  // vertex 0
        1.0, 0.0, 0.0,  // vertex 1
        0.0, 0.0, 1.0); // vertex 3

    auto tri_nodes = extractTriangleNodes<3>();

    compareFaceCoordinates(tet_face, tri_nodes,
                           "Tetrahedron face (y=0)", "Triangle K=3");
  }

  TEST(FaceConformity, Tetrahedron_Face023_MatchesTriangle_K3)
  {
    // Face with vertices 0,2,3: (0,0,0), (0,1,0), (0,0,1) on x=0 plane
    auto tet_face = extractTetrahedronTriangleFaceNodes<3>(
        0.0, 0.0, 0.0,  // vertex 0
        0.0, 1.0, 0.0,  // vertex 2
        0.0, 0.0, 1.0); // vertex 3

    auto tri_nodes = extractTriangleNodes<3>();

    compareFaceCoordinates(tet_face, tri_nodes,
                           "Tetrahedron face (x=0)", "Triangle K=3");
  }

  TEST(FaceConformity, Tetrahedron_Face123_MatchesTriangle_K3)
  {
    // Face with vertices 1,2,3: (1,0,0), (0,1,0), (0,0,1) diagonal face
    auto tet_face = extractTetrahedronTriangleFaceNodes<3>(
        1.0, 0.0, 0.0,  // vertex 1
        0.0, 1.0, 0.0,  // vertex 2
        0.0, 0.0, 1.0); // vertex 3

    auto tri_nodes = extractTriangleNodes<3>();

    compareFaceCoordinates(tet_face, tri_nodes,
                           "Tetrahedron diagonal face", "Triangle K=3");
  }

  //==========================================================================
  // Wedge triangular faces match Triangle (2D)
  // Reference wedge: triangle (0,0), (1,0), (0,1) extruded from z=0 to z=1
  // Bottom face (z=0): (0,0,0), (1,0,0), (0,1,0)
  // Top face (z=1): (0,0,1), (1,0,1), (0,1,1)
  //==========================================================================

  TEST(FaceConformity, Wedge_BottomFace_MatchesTriangle_K2)
  {
    // Bottom triangular face of wedge (z=0)
    auto wedge_face = extractWedgeTriangleFaceNodes<2>(
        0.0, 0.0, 0.0,  // vertex 0
        1.0, 0.0, 0.0,  // vertex 1
        0.0, 1.0, 0.0); // vertex 2

    auto tri_nodes = extractTriangleNodes<2>();

    compareFaceCoordinates(wedge_face, tri_nodes,
                           "Wedge bottom face", "Triangle");
  }

  TEST(FaceConformity, Wedge_BottomFace_MatchesTriangle_K3)
  {
    auto wedge_face = extractWedgeTriangleFaceNodes<3>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto tri_nodes = extractTriangleNodes<3>();

    compareFaceCoordinates(wedge_face, tri_nodes,
                           "Wedge bottom face K=3", "Triangle K=3");
  }

  TEST(FaceConformity, Wedge_BottomFace_MatchesTriangle_K4)
  {
    auto wedge_face = extractWedgeTriangleFaceNodes<4>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto tri_nodes = extractTriangleNodes<4>();

    compareFaceCoordinates(wedge_face, tri_nodes,
                           "Wedge bottom face K=4", "Triangle K=4");
  }

  TEST(FaceConformity, Wedge_TopFace_MatchesTriangle_K3)
  {
    // Top triangular face of wedge (z=1)
    auto wedge_face = extractWedgeTriangleFaceNodes<3>(
        0.0, 0.0, 1.0,  // vertex 3
        1.0, 0.0, 1.0,  // vertex 4
        0.0, 1.0, 1.0); // vertex 5

    auto tri_nodes = extractTriangleNodes<3>();

    compareFaceCoordinates(wedge_face, tri_nodes,
                           "Wedge top face K=3", "Triangle K=3");
  }

  TEST(FaceConformity, Wedge_TopFace_MatchesTriangle_K4)
  {
    auto wedge_face = extractWedgeTriangleFaceNodes<4>(
        0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,
        0.0, 1.0, 1.0);

    auto tri_nodes = extractTriangleNodes<4>();

    compareFaceCoordinates(wedge_face, tri_nodes,
                           "Wedge top face K=4", "Triangle K=4");
  }

  //==========================================================================
  // Tetrahedron triangular faces match Wedge triangular faces
  //==========================================================================

  TEST(FaceConformity, Tetrahedron_Wedge_TriangleFaces_Match_K2)
  {
    // Both have triangular faces - they should have the same node distribution
    auto tet_face = extractTetrahedronTriangleFaceNodes<2>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto wedge_face = extractWedgeTriangleFaceNodes<2>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    compareFaceCoordinates(tet_face, wedge_face,
                           "Tetrahedron triangle face", "Wedge triangle face");
  }

  TEST(FaceConformity, Tetrahedron_Wedge_TriangleFaces_Match_K3)
  {
    auto tet_face = extractTetrahedronTriangleFaceNodes<3>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto wedge_face = extractWedgeTriangleFaceNodes<3>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    compareFaceCoordinates(tet_face, wedge_face,
                           "Tetrahedron triangle face K=3", "Wedge triangle face K=3");
  }

  TEST(FaceConformity, Tetrahedron_Wedge_TriangleFaces_Match_K4)
  {
    auto tet_face = extractTetrahedronTriangleFaceNodes<4>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto wedge_face = extractWedgeTriangleFaceNodes<4>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    compareFaceCoordinates(tet_face, wedge_face,
                           "Tetrahedron triangle face K=4", "Wedge triangle face K=4");
  }

  TEST(FaceConformity, Tetrahedron_Wedge_TriangleFaces_Match_K5)
  {
    auto tet_face = extractTetrahedronTriangleFaceNodes<5>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto wedge_face = extractWedgeTriangleFaceNodes<5>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    compareFaceCoordinates(tet_face, wedge_face,
                           "Tetrahedron triangle face K=5", "Wedge triangle face K=5");
  }

  TEST(FaceConformity, Tetrahedron_Wedge_TriangleFaces_Match_K6)
  {
    auto tet_face = extractTetrahedronTriangleFaceNodes<6>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto wedge_face = extractWedgeTriangleFaceNodes<6>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    compareFaceCoordinates(tet_face, wedge_face,
                           "Tetrahedron triangle face K=6", "Wedge triangle face K=6");
  }

  //==========================================================================
  // Wedge quad faces match Quadrilateral (2D)
  // Reference wedge has 3 quadrilateral faces:
  // Face 1 (y=0): (0,0,0)-(1,0,0)-(1,0,1)-(0,0,1)
  // Face 2 (x=0): (0,0,0)-(0,1,0)-(0,1,1)-(0,0,1)
  // Face 3 (diagonal): (1,0,0)-(0,1,0)-(0,1,1)-(1,0,1)
  //==========================================================================

  TEST(FaceConformity, Wedge_QuadFaceY0_MatchesQuadrilateral_K2)
  {
    // Quad face at y=0
    auto wedge_face = extractWedgeQuadFaceNodes<2>(
        0.0, 0.0, 0.0,  // (0,0)
        1.0, 0.0, 0.0,  // (1,0)
        1.0, 0.0, 1.0,  // (1,1)
        0.0, 0.0, 1.0); // (0,1)

    auto quad_nodes = extractQuadrilateralNodes<2>();

    compareFaceCoordinates(wedge_face, quad_nodes,
                           "Wedge quad face (y=0)", "Quadrilateral");
  }

  TEST(FaceConformity, Wedge_QuadFaceY0_MatchesQuadrilateral_K3)
  {
    auto wedge_face = extractWedgeQuadFaceNodes<3>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 0.0, 1.0);

    auto quad_nodes = extractQuadrilateralNodes<3>();

    compareFaceCoordinates(wedge_face, quad_nodes,
                           "Wedge quad face (y=0) K=3", "Quadrilateral K=3");
  }

  TEST(FaceConformity, Wedge_QuadFaceY0_MatchesQuadrilateral_K4)
  {
    auto wedge_face = extractWedgeQuadFaceNodes<4>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 0.0, 1.0);

    auto quad_nodes = extractQuadrilateralNodes<4>();

    compareFaceCoordinates(wedge_face, quad_nodes,
                           "Wedge quad face (y=0) K=4", "Quadrilateral K=4");
  }

  TEST(FaceConformity, Wedge_QuadFaceX0_MatchesQuadrilateral_K3)
  {
    // Quad face at x=0
    auto wedge_face = extractWedgeQuadFaceNodes<3>(
        0.0, 0.0, 0.0,  // (0,0)
        0.0, 1.0, 0.0,  // (1,0)
        0.0, 1.0, 1.0,  // (1,1)
        0.0, 0.0, 1.0); // (0,1)

    auto quad_nodes = extractQuadrilateralNodes<3>();

    compareFaceCoordinates(wedge_face, quad_nodes,
                           "Wedge quad face (x=0) K=3", "Quadrilateral K=3");
  }

  TEST(FaceConformity, Wedge_QuadFaceX0_MatchesQuadrilateral_K4)
  {
    auto wedge_face = extractWedgeQuadFaceNodes<4>(
        0.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 1.0, 1.0,
        0.0, 0.0, 1.0);

    auto quad_nodes = extractQuadrilateralNodes<4>();

    compareFaceCoordinates(wedge_face, quad_nodes,
                           "Wedge quad face (x=0) K=4", "Quadrilateral K=4");
  }

  //==========================================================================
  // High-order conformity tests
  //==========================================================================

  TEST(FaceConformity, WedgeTriangleFaces_Match_K5)
  {
    // Reference triangle nodes
    auto tri_nodes = extractTriangleNodes<5>();

    // Wedge bottom face
    auto wedge_bottom = extractWedgeTriangleFaceNodes<5>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);
    compareFaceCoordinates(wedge_bottom, tri_nodes,
                           "Wedge bottom face K=5", "Triangle K=5");

    // Wedge top face
    auto wedge_top = extractWedgeTriangleFaceNodes<5>(
        0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,
        0.0, 1.0, 1.0);
    compareFaceCoordinates(wedge_top, tri_nodes,
                           "Wedge top face K=5", "Triangle K=5");
  }

  TEST(FaceConformity, WedgeTriangleFaces_Match_K6)
  {
    auto tri_nodes = extractTriangleNodes<6>();

    auto wedge_bottom = extractWedgeTriangleFaceNodes<6>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);
    compareFaceCoordinates(wedge_bottom, tri_nodes,
                           "Wedge bottom face K=6", "Triangle K=6");
  }

  TEST(FaceConformity, AllQuadFaces_Match_K5)
  {
    auto quad_nodes = extractQuadrilateralNodes<5>();

    // Wedge y=0 face
    auto wedge_y0 = extractWedgeQuadFaceNodes<5>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 0.0, 1.0);
    compareFaceCoordinates(wedge_y0, quad_nodes,
                           "Wedge quad face (y=0) K=5", "Quadrilateral K=5");

    // Wedge x=0 face
    auto wedge_x0 = extractWedgeQuadFaceNodes<5>(
        0.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 1.0, 1.0,
        0.0, 0.0, 1.0);
    compareFaceCoordinates(wedge_x0, quad_nodes,
                           "Wedge quad face (x=0) K=5", "Quadrilateral K=5");
  }

  TEST(FaceConformity, AllQuadFaces_Match_K6)
  {
    auto quad_nodes = extractQuadrilateralNodes<6>();

    auto wedge_y0 = extractWedgeQuadFaceNodes<6>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 0.0, 1.0);
    compareFaceCoordinates(wedge_y0, quad_nodes,
                           "Wedge quad face (y=0) K=6", "Quadrilateral K=6");
  }

  //==========================================================================
  // Cross-geometry comprehensive tests
  //==========================================================================

  TEST(FaceConformity, TriangleFaceNodeCount_K2)
  {
    // Triangle should have (K+1)(K+2)/2 = 6 nodes at K=2
    auto tri = extractTriangleNodes<2>();
    EXPECT_EQ(tri.size(), 6u);

    auto tet_face = extractTetrahedronTriangleFaceNodes<2>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);
    EXPECT_EQ(tet_face.size(), 6u);

    auto wedge_face = extractWedgeTriangleFaceNodes<2>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);
    EXPECT_EQ(wedge_face.size(), 6u);
  }

  TEST(FaceConformity, TriangleFaceNodeCount_K3)
  {
    // Triangle should have (K+1)(K+2)/2 = 10 nodes at K=3
    auto tri = extractTriangleNodes<3>();
    EXPECT_EQ(tri.size(), 10u);

    auto tet_face = extractTetrahedronTriangleFaceNodes<3>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);
    EXPECT_EQ(tet_face.size(), 10u);

    auto wedge_face = extractWedgeTriangleFaceNodes<3>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);
    EXPECT_EQ(wedge_face.size(), 10u);
  }

  TEST(FaceConformity, QuadFaceNodeCount_K2)
  {
    // Quadrilateral should have (K+1)^2 = 9 nodes at K=2
    auto quad = extractQuadrilateralNodes<2>();
    EXPECT_EQ(quad.size(), 9u);

    auto wedge_face = extractWedgeQuadFaceNodes<2>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 0.0, 1.0);
    EXPECT_EQ(wedge_face.size(), 9u);
  }

  TEST(FaceConformity, QuadFaceNodeCount_K3)
  {
    // Quadrilateral should have (K+1)^2 = 16 nodes at K=3
    auto quad = extractQuadrilateralNodes<3>();
    EXPECT_EQ(quad.size(), 16u);

    auto wedge_face = extractWedgeQuadFaceNodes<3>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 0.0, 1.0);
    EXPECT_EQ(wedge_face.size(), 16u);
  }

  //==========================================================================
  // Very high-order conformity tests (K = 15)
  //==========================================================================

  TEST(FaceConformity, Tetrahedron_BaseFace_MatchesTriangle_K15)
  {
    auto tet_face = extractTetrahedronTriangleFaceNodes<15>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto tri_nodes = extractTriangleNodes<15>();

    compareFaceCoordinates(tet_face, tri_nodes,
                           "Tetrahedron base face K=15", "Triangle K=15");
  }

  TEST(FaceConformity, Tetrahedron_AllFaces_MatchTriangle_K15)
  {
    auto tri_nodes = extractTriangleNodes<15>();

    // Base face (z=0)
    auto tet_base = extractTetrahedronTriangleFaceNodes<15>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);
    compareFaceCoordinates(tet_base, tri_nodes,
                           "Tetrahedron base face K=15", "Triangle K=15");

    // Face (y=0)
    auto tet_y0 = extractTetrahedronTriangleFaceNodes<15>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 0.0, 1.0);
    compareFaceCoordinates(tet_y0, tri_nodes,
                           "Tetrahedron face (y=0) K=15", "Triangle K=15");

    // Face (x=0)
    auto tet_x0 = extractTetrahedronTriangleFaceNodes<15>(
        0.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0);
    compareFaceCoordinates(tet_x0, tri_nodes,
                           "Tetrahedron face (x=0) K=15", "Triangle K=15");

    // Diagonal face
    auto tet_diag = extractTetrahedronTriangleFaceNodes<15>(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0);
    compareFaceCoordinates(tet_diag, tri_nodes,
                           "Tetrahedron diagonal face K=15", "Triangle K=15");
  }

  TEST(FaceConformity, Wedge_TriangleFaces_MatchTriangle_K15)
  {
    auto tri_nodes = extractTriangleNodes<15>();

    // Bottom face (z=0)
    auto wedge_bottom = extractWedgeTriangleFaceNodes<15>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);
    compareFaceCoordinates(wedge_bottom, tri_nodes,
                           "Wedge bottom face K=15", "Triangle K=15");

    // Top face (z=1)
    auto wedge_top = extractWedgeTriangleFaceNodes<15>(
        0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,
        0.0, 1.0, 1.0);
    compareFaceCoordinates(wedge_top, tri_nodes,
                           "Wedge top face K=15", "Triangle K=15");
  }

  TEST(FaceConformity, Tetrahedron_Wedge_TriangleFaces_Match_K15)
  {
    auto tet_face = extractTetrahedronTriangleFaceNodes<15>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    auto wedge_face = extractWedgeTriangleFaceNodes<15>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    compareFaceCoordinates(tet_face, wedge_face,
                           "Tetrahedron triangle face K=15", "Wedge triangle face K=15");
  }

  TEST(FaceConformity, Wedge_QuadFaces_MatchQuadrilateral_K15)
  {
    auto quad_nodes = extractQuadrilateralNodes<15>();

    // y=0 face
    auto wedge_y0 = extractWedgeQuadFaceNodes<15>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 0.0, 1.0);
    compareFaceCoordinates(wedge_y0, quad_nodes,
                           "Wedge quad face (y=0) K=15", "Quadrilateral K=15");

    // x=0 face
    auto wedge_x0 = extractWedgeQuadFaceNodes<15>(
        0.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 1.0, 1.0,
        0.0, 0.0, 1.0);
    compareFaceCoordinates(wedge_x0, quad_nodes,
                           "Wedge quad face (x=0) K=15", "Quadrilateral K=15");
  }

  TEST(FaceConformity, TriangleFaceNodeCount_K15)
  {
    // Triangle should have (K+1)(K+2)/2 = 136 nodes at K=15
    auto tri = extractTriangleNodes<15>();
    EXPECT_EQ(tri.size(), 136u);

    auto tet_face = extractTetrahedronTriangleFaceNodes<15>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);
    EXPECT_EQ(tet_face.size(), 136u);

    auto wedge_face = extractWedgeTriangleFaceNodes<15>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0);
    EXPECT_EQ(wedge_face.size(), 136u);
  }

  TEST(FaceConformity, QuadFaceNodeCount_K15)
  {
    // Quadrilateral should have (K+1)^2 = 256 nodes at K=15
    auto quad = extractQuadrilateralNodes<15>();
    EXPECT_EQ(quad.size(), 256u);

    auto wedge_face = extractWedgeQuadFaceNodes<15>(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 0.0, 1.0);
    EXPECT_EQ(wedge_face.size(), 256u);
  }
}
