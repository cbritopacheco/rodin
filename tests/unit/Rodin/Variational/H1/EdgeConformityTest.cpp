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
  // Common constants for edge extraction
  //==========================================================================

  // Tolerance for geometric comparisons (whether a point is on an edge)
  constexpr Real EDGE_TOL = 1e-10;

  // Epsilon for division by zero protection
  const Real EDGE_EPS = std::numeric_limits<Real>::epsilon() * 100;

  //==========================================================================
  // Helper functions for extracting edge nodes
  //==========================================================================

  // Compute parameter t along edge from (start) to (end) for point p
  // Returns t in [0,1] if point is on edge, or NaN if edge is degenerate
  inline Real computeEdgeParameter2D(
      Real px, Real py,
      Real x0, Real y0,
      Real x1, Real y1)
  {
    Real dx = x1 - x0;
    Real dy = y1 - y0;

    Real t;
    if (std::abs(dx) > EDGE_EPS && std::abs(dx) >= std::abs(dy))
    {
      t = (px - x0) / dx;
    }
    else if (std::abs(dy) > EDGE_EPS)
    {
      t = (py - y0) / dy;
    }
    else
    {
      return std::numeric_limits<Real>::quiet_NaN();
    }

    // Verify point is actually on the line
    Real expected_x = x0 + t * dx;
    Real expected_y = y0 + t * dy;

    if (std::abs(px - expected_x) < EDGE_TOL && std::abs(py - expected_y) < EDGE_TOL)
    {
      return t;
    }
    return std::numeric_limits<Real>::quiet_NaN();
  }

  // Compute parameter t along edge in 3D
  inline Real computeEdgeParameter3D(
      Real px, Real py, Real pz,
      Real x0, Real y0, Real z0,
      Real x1, Real y1, Real z1)
  {
    Real dx = x1 - x0;
    Real dy = y1 - y0;
    Real dz = z1 - z0;

    Real t;
    if (std::abs(dx) > EDGE_EPS && std::abs(dx) >= std::abs(dy) && std::abs(dx) >= std::abs(dz))
    {
      t = (px - x0) / dx;
    }
    else if (std::abs(dy) > EDGE_EPS && std::abs(dy) >= std::abs(dz))
    {
      t = (py - y0) / dy;
    }
    else if (std::abs(dz) > EDGE_EPS)
    {
      t = (pz - z0) / dz;
    }
    else
    {
      return std::numeric_limits<Real>::quiet_NaN();
    }

    // Verify point is actually on the line
    Real expected_x = x0 + t * dx;
    Real expected_y = y0 + t * dy;
    Real expected_z = z0 + t * dz;

    if (std::abs(px - expected_x) < EDGE_TOL &&
        std::abs(py - expected_y) < EDGE_TOL &&
        std::abs(pz - expected_z) < EDGE_TOL)
    {
      return t;
    }
    return std::numeric_limits<Real>::quiet_NaN();
  }

  // Extract 1D coordinates of nodes lying on a specific edge
  // Returns sorted coordinates in [0,1]
  template <size_t K>
  std::vector<Real> extractSegmentNodes()
  {
    std::vector<Real> coords;
    const auto& nodes = H1Element<K, Real>::getNodes(Polytope::Type::Segment);
    for (const auto& node : nodes)
    {
      coords.push_back(node.x());
    }
    std::sort(coords.begin(), coords.end());
    return coords;
  }

  // Extract nodes on Triangle edge from vertex (x0,y0) to (x1,y1)
  template <size_t K>
  std::vector<Real> extractTriangleEdgeNodes(Real x0, Real y0, Real x1, Real y1)
  {
    std::vector<Real> coords;
    const auto& nodes = FeketeTriangle<K>::getNodes();

    for (const auto& node : nodes)
    {
      Real t = computeEdgeParameter2D(node.x(), node.y(), x0, y0, x1, y1);
      if (!std::isnan(t) && t >= -EDGE_TOL && t <= 1.0 + EDGE_TOL)
      {
        coords.push_back(std::max(0.0, std::min(1.0, t)));
      }
    }

    std::sort(coords.begin(), coords.end());
    return coords;
  }

  // Extract nodes on Quadrilateral edge
  template <size_t K>
  std::vector<Real> extractQuadrilateralEdgeNodes(Real x0, Real y0, Real x1, Real y1)
  {
    std::vector<Real> coords;
    const auto& nodes = H1Element<K, Real>::getNodes(Polytope::Type::Quadrilateral);

    for (const auto& node : nodes)
    {
      Real t = computeEdgeParameter2D(node.x(), node.y(), x0, y0, x1, y1);
      if (!std::isnan(t) && t >= -EDGE_TOL && t <= 1.0 + EDGE_TOL)
      {
        coords.push_back(std::max(0.0, std::min(1.0, t)));
      }
    }

    std::sort(coords.begin(), coords.end());
    return coords;
  }

  // Extract nodes on Tetrahedron edge
  template <size_t K>
  std::vector<Real> extractTetrahedronEdgeNodes(
      Real x0, Real y0, Real z0,
      Real x1, Real y1, Real z1)
  {
    std::vector<Real> coords;
    const auto& nodes = FeketeTetrahedron<K>::getNodes();

    for (const auto& node : nodes)
    {
      Real t = computeEdgeParameter3D(node.x(), node.y(), node.z(), x0, y0, z0, x1, y1, z1);
      if (!std::isnan(t) && t >= -EDGE_TOL && t <= 1.0 + EDGE_TOL)
      {
        coords.push_back(std::max(0.0, std::min(1.0, t)));
      }
    }

    std::sort(coords.begin(), coords.end());
    return coords;
  }

  // Extract nodes on Wedge edge
  template <size_t K>
  std::vector<Real> extractWedgeEdgeNodes(
      Real x0, Real y0, Real z0,
      Real x1, Real y1, Real z1)
  {
    std::vector<Real> coords;
    const auto& nodes = H1Element<K, Real>::getNodes(Polytope::Type::Wedge);

    for (const auto& node : nodes)
    {
      Real t = computeEdgeParameter3D(node.x(), node.y(), node.z(), x0, y0, z0, x1, y1, z1);
      if (!std::isnan(t) && t >= -EDGE_TOL && t <= 1.0 + EDGE_TOL)
      {
        coords.push_back(std::max(0.0, std::min(1.0, t)));
      }
    }

    std::sort(coords.begin(), coords.end());
    return coords;
  }

  // Compare two vectors of 1D coordinates
  inline void compareEdgeCoordinates(
      const std::vector<Real>& edge1,
      const std::vector<Real>& edge2,
      const std::string& name1,
      const std::string& name2)
  {
    ASSERT_EQ(edge1.size(), edge2.size())
        << "Edge node counts differ: " << name1 << " has " << edge1.size()
        << " nodes, " << name2 << " has " << edge2.size() << " nodes";

    for (size_t i = 0; i < edge1.size(); ++i)
    {
      EXPECT_NEAR(edge1[i], edge2[i], EDGE_TOL)
          << "Node " << i << " differs: " << name1 << "[" << i << "]=" << edge1[i]
          << ", " << name2 << "[" << i << "]=" << edge2[i];
    }
  }

  //==========================================================================
  // Segment matches GLL01 nodes
  //==========================================================================

  TEST(EdgeConformity, Segment_MatchesGLL01_K2)
  {
    auto seg_nodes = extractSegmentNodes<2>();
    const auto& gll = GLL01<2>::getNodes();

    ASSERT_EQ(seg_nodes.size(), 3);
    for (size_t i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(seg_nodes[i], gll[i], 1e-14);
    }
  }

  TEST(EdgeConformity, Segment_MatchesGLL01_K3)
  {
    auto seg_nodes = extractSegmentNodes<3>();
    const auto& gll = GLL01<3>::getNodes();

    ASSERT_EQ(seg_nodes.size(), 4);
    for (size_t i = 0; i < 4; ++i)
    {
      EXPECT_NEAR(seg_nodes[i], gll[i], 1e-14);
    }
  }

  TEST(EdgeConformity, Segment_MatchesGLL01_K5)
  {
    auto seg_nodes = extractSegmentNodes<5>();
    const auto& gll = GLL01<5>::getNodes();

    ASSERT_EQ(seg_nodes.size(), 6);
    for (size_t i = 0; i < 6; ++i)
    {
      EXPECT_NEAR(seg_nodes[i], gll[i], 1e-14);
    }
  }

  //==========================================================================
  // Triangle edges match GLL01 nodes
  //==========================================================================

  TEST(EdgeConformity, Triangle_Edge01_MatchesGLL01_K2)
  {
    // Edge from (0,0) to (1,0)
    auto tri_edge = extractTriangleEdgeNodes<2>(0.0, 0.0, 1.0, 0.0);
    const auto& gll = GLL01<2>::getNodes();

    std::vector<Real> gll_vec(gll.begin(), gll.end());
    compareEdgeCoordinates(tri_edge, gll_vec, "Triangle edge (0,0)-(1,0)", "GLL01");
  }

  TEST(EdgeConformity, Triangle_Edge02_MatchesGLL01_K2)
  {
    // Edge from (0,0) to (0,1)
    auto tri_edge = extractTriangleEdgeNodes<2>(0.0, 0.0, 0.0, 1.0);
    const auto& gll = GLL01<2>::getNodes();

    std::vector<Real> gll_vec(gll.begin(), gll.end());
    compareEdgeCoordinates(tri_edge, gll_vec, "Triangle edge (0,0)-(0,1)", "GLL01");
  }

  TEST(EdgeConformity, Triangle_Edge12_MatchesGLL01_K2)
  {
    // Edge from (1,0) to (0,1)
    auto tri_edge = extractTriangleEdgeNodes<2>(1.0, 0.0, 0.0, 1.0);
    const auto& gll = GLL01<2>::getNodes();

    std::vector<Real> gll_vec(gll.begin(), gll.end());
    compareEdgeCoordinates(tri_edge, gll_vec, "Triangle edge (1,0)-(0,1)", "GLL01");
  }

  TEST(EdgeConformity, Triangle_AllEdges_MatchGLL01_K3)
  {
    const auto& gll = GLL01<3>::getNodes();
    std::vector<Real> gll_vec(gll.begin(), gll.end());

    // Edge (0,0) to (1,0)
    auto edge01 = extractTriangleEdgeNodes<3>(0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(edge01, gll_vec, "Triangle edge (0,0)-(1,0)", "GLL01 K=3");

    // Edge (0,0) to (0,1)
    auto edge02 = extractTriangleEdgeNodes<3>(0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge02, gll_vec, "Triangle edge (0,0)-(0,1)", "GLL01 K=3");

    // Edge (1,0) to (0,1)
    auto edge12 = extractTriangleEdgeNodes<3>(1.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge12, gll_vec, "Triangle edge (1,0)-(0,1)", "GLL01 K=3");
  }

  TEST(EdgeConformity, Triangle_AllEdges_MatchGLL01_K5)
  {
    const auto& gll = GLL01<5>::getNodes();
    std::vector<Real> gll_vec(gll.begin(), gll.end());

    auto edge01 = extractTriangleEdgeNodes<5>(0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(edge01, gll_vec, "Triangle edge (0,0)-(1,0)", "GLL01 K=5");

    auto edge02 = extractTriangleEdgeNodes<5>(0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge02, gll_vec, "Triangle edge (0,0)-(0,1)", "GLL01 K=5");

    auto edge12 = extractTriangleEdgeNodes<5>(1.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge12, gll_vec, "Triangle edge (1,0)-(0,1)", "GLL01 K=5");
  }

  //==========================================================================
  // Quadrilateral edges match GLL01 nodes
  //==========================================================================

  TEST(EdgeConformity, Quadrilateral_AllEdges_MatchGLL01_K2)
  {
    const auto& gll = GLL01<2>::getNodes();
    std::vector<Real> gll_vec(gll.begin(), gll.end());

    // Bottom edge: (0,0) to (1,0)
    auto edge_bottom = extractQuadrilateralEdgeNodes<2>(0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(edge_bottom, gll_vec, "Quad bottom edge", "GLL01 K=2");

    // Right edge: (1,0) to (1,1)
    auto edge_right = extractQuadrilateralEdgeNodes<2>(1.0, 0.0, 1.0, 1.0);
    compareEdgeCoordinates(edge_right, gll_vec, "Quad right edge", "GLL01 K=2");

    // Top edge: (0,1) to (1,1)
    auto edge_top = extractQuadrilateralEdgeNodes<2>(0.0, 1.0, 1.0, 1.0);
    compareEdgeCoordinates(edge_top, gll_vec, "Quad top edge", "GLL01 K=2");

    // Left edge: (0,1) to (0,0)
    auto edge_left = extractQuadrilateralEdgeNodes<2>(0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge_left, gll_vec, "Quad left edge", "GLL01 K=2");
  }

  TEST(EdgeConformity, Quadrilateral_AllEdges_MatchGLL01_K4)
  {
    const auto& gll = GLL01<4>::getNodes();
    std::vector<Real> gll_vec(gll.begin(), gll.end());

    auto edge_bottom = extractQuadrilateralEdgeNodes<4>(0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(edge_bottom, gll_vec, "Quad bottom edge K=4", "GLL01");

    auto edge_right = extractQuadrilateralEdgeNodes<4>(1.0, 0.0, 1.0, 1.0);
    compareEdgeCoordinates(edge_right, gll_vec, "Quad right edge K=4", "GLL01");

    auto edge_top = extractQuadrilateralEdgeNodes<4>(0.0, 1.0, 1.0, 1.0);
    compareEdgeCoordinates(edge_top, gll_vec, "Quad top edge K=4", "GLL01");

    auto edge_left = extractQuadrilateralEdgeNodes<4>(0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge_left, gll_vec, "Quad left edge K=4", "GLL01");
  }

  //==========================================================================
  // Tetrahedron edges match GLL01 nodes
  //==========================================================================

  TEST(EdgeConformity, Tetrahedron_AllEdges_MatchGLL01_K2)
  {
    const auto& gll = GLL01<2>::getNodes();
    std::vector<Real> gll_vec(gll.begin(), gll.end());

    // Reference tetrahedron vertices: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    // Edge 01: (0,0,0) to (1,0,0)
    auto edge01 = extractTetrahedronEdgeNodes<2>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(edge01, gll_vec, "Tet edge 01", "GLL01 K=2");

    // Edge 02: (0,0,0) to (0,1,0)
    auto edge02 = extractTetrahedronEdgeNodes<2>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(edge02, gll_vec, "Tet edge 02", "GLL01 K=2");

    // Edge 03: (0,0,0) to (0,0,1)
    auto edge03 = extractTetrahedronEdgeNodes<2>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge03, gll_vec, "Tet edge 03", "GLL01 K=2");

    // Edge 12: (1,0,0) to (0,1,0)
    auto edge12 = extractTetrahedronEdgeNodes<2>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(edge12, gll_vec, "Tet edge 12", "GLL01 K=2");

    // Edge 13: (1,0,0) to (0,0,1)
    auto edge13 = extractTetrahedronEdgeNodes<2>(1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge13, gll_vec, "Tet edge 13", "GLL01 K=2");

    // Edge 23: (0,1,0) to (0,0,1)
    auto edge23 = extractTetrahedronEdgeNodes<2>(0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge23, gll_vec, "Tet edge 23", "GLL01 K=2");
  }

  TEST(EdgeConformity, Tetrahedron_AllEdges_MatchGLL01_K4)
  {
    const auto& gll = GLL01<4>::getNodes();
    std::vector<Real> gll_vec(gll.begin(), gll.end());

    auto edge01 = extractTetrahedronEdgeNodes<4>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(edge01, gll_vec, "Tet edge 01 K=4", "GLL01");

    auto edge02 = extractTetrahedronEdgeNodes<4>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(edge02, gll_vec, "Tet edge 02 K=4", "GLL01");

    auto edge03 = extractTetrahedronEdgeNodes<4>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge03, gll_vec, "Tet edge 03 K=4", "GLL01");

    auto edge12 = extractTetrahedronEdgeNodes<4>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(edge12, gll_vec, "Tet edge 12 K=4", "GLL01");

    auto edge13 = extractTetrahedronEdgeNodes<4>(1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge13, gll_vec, "Tet edge 13 K=4", "GLL01");

    auto edge23 = extractTetrahedronEdgeNodes<4>(0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge23, gll_vec, "Tet edge 23 K=4", "GLL01");
  }

  //==========================================================================
  // Wedge edges match appropriately
  //==========================================================================

  TEST(EdgeConformity, Wedge_VerticalEdges_MatchGLL01_K2)
  {
    // Wedge: triangle at z=0 and z=1, reference triangle (0,0)-(1,0)-(0,1)
    // Vertical edges (z direction) should match GLL01 nodes
    const auto& gll = GLL01<2>::getNodes();
    std::vector<Real> gll_vec(gll.begin(), gll.end());

    // Vertical edge at (0,0)
    auto edge_v1 = extractWedgeEdgeNodes<2>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(edge_v1, gll_vec, "Wedge vertical edge at (0,0)", "GLL01 K=2");

    // Vertical edge at (1,0)
    auto edge_v2 = extractWedgeEdgeNodes<2>(1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
    compareEdgeCoordinates(edge_v2, gll_vec, "Wedge vertical edge at (1,0)", "GLL01 K=2");

    // Vertical edge at (0,1)
    auto edge_v3 = extractWedgeEdgeNodes<2>(0.0, 1.0, 0.0, 0.0, 1.0, 1.0);
    compareEdgeCoordinates(edge_v3, gll_vec, "Wedge vertical edge at (0,1)", "GLL01 K=2");
  }

  TEST(EdgeConformity, Wedge_TriangleEdges_MatchTriangle_K3)
  {
    // Triangle edges on z=0 face should match the triangle Fekete edge nodes
    const auto& gll = GLL01<3>::getNodes();
    std::vector<Real> gll_vec(gll.begin(), gll.end());

    // Edge (0,0,0) to (1,0,0) on bottom triangle
    auto edge01 = extractWedgeEdgeNodes<3>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(edge01, gll_vec, "Wedge bottom triangle edge 01", "GLL01 K=3");

    // Edge (0,0,0) to (0,1,0) on bottom triangle
    auto edge02 = extractWedgeEdgeNodes<3>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(edge02, gll_vec, "Wedge bottom triangle edge 02", "GLL01 K=3");

    // Edge (1,0,0) to (0,1,0) on bottom triangle
    auto edge12 = extractWedgeEdgeNodes<3>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(edge12, gll_vec, "Wedge bottom triangle edge 12", "GLL01 K=3");

    // Same edges on top triangle (z=1)
    auto top_edge01 = extractWedgeEdgeNodes<3>(0.0, 0.0, 1.0, 1.0, 0.0, 1.0);
    compareEdgeCoordinates(top_edge01, gll_vec, "Wedge top triangle edge 01", "GLL01 K=3");
  }

  //==========================================================================
  // Cross-geometry conformity tests
  //==========================================================================

  TEST(EdgeConformity, Triangle_Quadrilateral_SharedEdge_K3)
  {
    // Both triangle and quadrilateral have edge from (0,0) to (1,0)
    // which should match and use GLL01 nodes
    auto tri_edge = extractTriangleEdgeNodes<3>(0.0, 0.0, 1.0, 0.0);
    auto quad_edge = extractQuadrilateralEdgeNodes<3>(0.0, 0.0, 1.0, 0.0);

    compareEdgeCoordinates(tri_edge, quad_edge, "Triangle edge", "Quadrilateral edge");
  }

  TEST(EdgeConformity, Triangle_Tetrahedron_SharedEdge_K3)
  {
    // Triangle edge (0,0)-(1,0) should match tetrahedron edge (0,0,0)-(1,0,0)
    auto tri_edge = extractTriangleEdgeNodes<3>(0.0, 0.0, 1.0, 0.0);
    auto tet_edge = extractTetrahedronEdgeNodes<3>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);

    compareEdgeCoordinates(tri_edge, tet_edge, "Triangle edge", "Tetrahedron edge");
  }

  TEST(EdgeConformity, Triangle_Wedge_SharedEdge_K4)
  {
    // Triangle edge (0,0)-(1,0) should match wedge edge (0,0,0)-(1,0,0)
    auto tri_edge = extractTriangleEdgeNodes<4>(0.0, 0.0, 1.0, 0.0);
    auto wedge_edge = extractWedgeEdgeNodes<4>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);

    compareEdgeCoordinates(tri_edge, wedge_edge, "Triangle edge", "Wedge triangle edge");
  }

  TEST(EdgeConformity, Segment_Quadrilateral_SharedEdge_K4)
  {
    // Segment nodes should match quadrilateral bottom edge
    auto seg = extractSegmentNodes<4>();
    auto quad_edge = extractQuadrilateralEdgeNodes<4>(0.0, 0.0, 1.0, 0.0);

    compareEdgeCoordinates(seg, quad_edge, "Segment", "Quadrilateral bottom edge");
  }

  TEST(EdgeConformity, Segment_Wedge_VerticalEdge_K5)
  {
    // Segment nodes should match wedge vertical edge
    auto seg = extractSegmentNodes<5>();
    auto wedge_edge = extractWedgeEdgeNodes<5>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

    compareEdgeCoordinates(seg, wedge_edge, "Segment", "Wedge vertical edge");
  }

  //==========================================================================
  // Higher order conformity tests (K = 6)
  //==========================================================================

  TEST(EdgeConformity, AllGeometries_EdgeNodes_Match_K6)
  {
    const auto& gll = GLL01<6>::getNodes();
    std::vector<Real> gll_vec(gll.begin(), gll.end());

    // Segment
    auto seg = extractSegmentNodes<6>();
    compareEdgeCoordinates(seg, gll_vec, "Segment K=6", "GLL01");

    // Triangle
    auto tri_edge = extractTriangleEdgeNodes<6>(0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tri_edge, gll_vec, "Triangle edge K=6", "GLL01");

    // Quadrilateral
    auto quad_edge = extractQuadrilateralEdgeNodes<6>(0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(quad_edge, gll_vec, "Quadrilateral edge K=6", "GLL01");

    // Tetrahedron
    auto tet_edge = extractTetrahedronEdgeNodes<6>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(tet_edge, gll_vec, "Tetrahedron edge K=6", "GLL01");

    // Wedge - triangle edge
    auto wedge_tri_edge = extractWedgeEdgeNodes<6>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(wedge_tri_edge, gll_vec, "Wedge triangle edge K=6", "GLL01");

    // Wedge - vertical edge
    auto wedge_vert_edge = extractWedgeEdgeNodes<6>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(wedge_vert_edge, gll_vec, "Wedge vertical edge K=6", "GLL01");
  }

  //==========================================================================
  // 3D-3D Cross-geometry conformity tests (Tetrahedron vs Wedge)
  //==========================================================================

  TEST(EdgeConformity3D, Tetrahedron_Wedge_SharedEdge_K3)
  {
    // Both tetrahedron and wedge share edge (0,0,0)-(1,0,0) along x-axis
    auto tet_edge = extractTetrahedronEdgeNodes<3>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    auto wedge_edge = extractWedgeEdgeNodes<3>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);

    compareEdgeCoordinates(tet_edge, wedge_edge, "Tetrahedron edge", "Wedge edge");
  }

  TEST(EdgeConformity3D, Tetrahedron_Wedge_SharedEdge_YAxis_K3)
  {
    // Both tetrahedron and wedge share edge (0,0,0)-(0,1,0) along y-axis
    auto tet_edge = extractTetrahedronEdgeNodes<3>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    auto wedge_edge = extractWedgeEdgeNodes<3>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    compareEdgeCoordinates(tet_edge, wedge_edge, "Tetrahedron edge (y-axis)", "Wedge edge (y-axis)");
  }

  TEST(EdgeConformity3D, Tetrahedron_Wedge_DiagonalEdge_K3)
  {
    // Both share the diagonal edge (1,0,0)-(0,1,0) on the xy-plane face
    auto tet_edge = extractTetrahedronEdgeNodes<3>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    auto wedge_edge = extractWedgeEdgeNodes<3>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    compareEdgeCoordinates(tet_edge, wedge_edge, "Tetrahedron diagonal edge", "Wedge diagonal edge");
  }

  TEST(EdgeConformity3D, Tetrahedron_Wedge_AllSharedEdges_K4)
  {
    // Test all edges that can be shared between tetrahedron and wedge base face
    // Reference tetrahedron: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    // Reference wedge base: triangle at z=0 with vertices (0,0,0), (1,0,0), (0,1,0)

    // Edge 01: (0,0,0) to (1,0,0)
    auto tet_01 = extractTetrahedronEdgeNodes<4>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    auto wedge_01 = extractWedgeEdgeNodes<4>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(tet_01, wedge_01, "Tet-Wedge edge 01 K=4", "same");

    // Edge 02: (0,0,0) to (0,1,0)
    auto tet_02 = extractTetrahedronEdgeNodes<4>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    auto wedge_02 = extractWedgeEdgeNodes<4>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tet_02, wedge_02, "Tet-Wedge edge 02 K=4", "same");

    // Edge 12: (1,0,0) to (0,1,0)
    auto tet_12 = extractTetrahedronEdgeNodes<4>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    auto wedge_12 = extractWedgeEdgeNodes<4>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tet_12, wedge_12, "Tet-Wedge edge 12 K=4", "same");
  }

  TEST(EdgeConformity3D, Tetrahedron_Wedge_HighOrder_K5)
  {
    // High-order test for tetrahedron-wedge shared edges
    auto tet_edge = extractTetrahedronEdgeNodes<5>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    auto wedge_edge = extractWedgeEdgeNodes<5>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);

    compareEdgeCoordinates(tet_edge, wedge_edge, "Tetrahedron edge K=5", "Wedge edge K=5");
  }

  TEST(EdgeConformity3D, Tetrahedron_Wedge_HighOrder_K6)
  {
    // All shared edges at K=6
    auto tet_01 = extractTetrahedronEdgeNodes<6>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    auto wedge_01 = extractWedgeEdgeNodes<6>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(tet_01, wedge_01, "Tet-Wedge edge 01 K=6", "same");

    auto tet_02 = extractTetrahedronEdgeNodes<6>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    auto wedge_02 = extractWedgeEdgeNodes<6>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tet_02, wedge_02, "Tet-Wedge edge 02 K=6", "same");

    auto tet_12 = extractTetrahedronEdgeNodes<6>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    auto wedge_12 = extractWedgeEdgeNodes<6>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tet_12, wedge_12, "Tet-Wedge edge 12 K=6", "same");
  }

  //==========================================================================
  // 3D-2D Cross-geometry conformity tests
  //==========================================================================

  TEST(EdgeConformity3D2D, Triangle_Tetrahedron_AllSharedEdges_K3)
  {
    // Triangle can be a face of tetrahedron - all 3 edges should match
    // Reference triangle: (0,0), (1,0), (0,1)
    // Tetrahedron base face: (0,0,0), (1,0,0), (0,1,0)

    // Edge 01
    auto tri_01 = extractTriangleEdgeNodes<3>(0.0, 0.0, 1.0, 0.0);
    auto tet_01 = extractTetrahedronEdgeNodes<3>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(tri_01, tet_01, "Triangle edge 01", "Tetrahedron edge 01");

    // Edge 02
    auto tri_02 = extractTriangleEdgeNodes<3>(0.0, 0.0, 0.0, 1.0);
    auto tet_02 = extractTetrahedronEdgeNodes<3>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tri_02, tet_02, "Triangle edge 02", "Tetrahedron edge 02");

    // Edge 12
    auto tri_12 = extractTriangleEdgeNodes<3>(1.0, 0.0, 0.0, 1.0);
    auto tet_12 = extractTetrahedronEdgeNodes<3>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tri_12, tet_12, "Triangle edge 12", "Tetrahedron edge 12");
  }

  TEST(EdgeConformity3D2D, Triangle_Wedge_AllSharedEdges_K3)
  {
    // Triangle can be top or bottom face of wedge - all 3 edges should match
    // Reference triangle: (0,0), (1,0), (0,1)
    // Wedge bottom face: (0,0,0), (1,0,0), (0,1,0)

    // Edge 01
    auto tri_01 = extractTriangleEdgeNodes<3>(0.0, 0.0, 1.0, 0.0);
    auto wedge_01 = extractWedgeEdgeNodes<3>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(tri_01, wedge_01, "Triangle edge 01", "Wedge bottom edge 01");

    // Edge 02
    auto tri_02 = extractTriangleEdgeNodes<3>(0.0, 0.0, 0.0, 1.0);
    auto wedge_02 = extractWedgeEdgeNodes<3>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tri_02, wedge_02, "Triangle edge 02", "Wedge bottom edge 02");

    // Edge 12
    auto tri_12 = extractTriangleEdgeNodes<3>(1.0, 0.0, 0.0, 1.0);
    auto wedge_12 = extractWedgeEdgeNodes<3>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tri_12, wedge_12, "Triangle edge 12", "Wedge bottom edge 12");

    // Also test top face (z=1)
    auto wedge_top_01 = extractWedgeEdgeNodes<3>(0.0, 0.0, 1.0, 1.0, 0.0, 1.0);
    compareEdgeCoordinates(tri_01, wedge_top_01, "Triangle edge 01", "Wedge top edge 01");
  }

  TEST(EdgeConformity3D2D, Quadrilateral_Wedge_QuadFaces_K3)
  {
    // Wedge has 3 quadrilateral faces - edges on these should match quadrilateral edges
    // Wedge quad face 1: (0,0,0)-(1,0,0)-(1,0,1)-(0,0,1) (y=0 plane)

    // Bottom edge of quad face (y=0 plane): (0,0,0) to (1,0,0)
    auto quad_bottom = extractQuadrilateralEdgeNodes<3>(0.0, 0.0, 1.0, 0.0);
    auto wedge_edge = extractWedgeEdgeNodes<3>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(quad_bottom, wedge_edge, "Quad bottom edge", "Wedge face edge (x-dir)");

    // Vertical edge of quad face: (0,0,0) to (0,0,1)
    auto quad_left = extractQuadrilateralEdgeNodes<3>(0.0, 0.0, 0.0, 1.0);
    auto wedge_vert = extractWedgeEdgeNodes<3>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(quad_left, wedge_vert, "Quad left edge", "Wedge vertical edge");
  }

  TEST(EdgeConformity3D2D, Quadrilateral_Tetrahedron_NoSharedFace)
  {
    // Quadrilateral cannot be a face of tetrahedron, but edges can still match
    // along the xy-plane if we consider the same edge positions

    // Bottom edge of quad (0,0)-(1,0) matches tet edge (0,0,0)-(1,0,0)
    auto quad_edge = extractQuadrilateralEdgeNodes<3>(0.0, 0.0, 1.0, 0.0);
    auto tet_edge = extractTetrahedronEdgeNodes<3>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(quad_edge, tet_edge, "Quad edge", "Tet edge");

    // Left edge of quad (0,0)-(0,1) matches tet edge (0,0,0)-(0,1,0)
    auto quad_left = extractQuadrilateralEdgeNodes<3>(0.0, 0.0, 0.0, 1.0);
    auto tet_left = extractTetrahedronEdgeNodes<3>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(quad_left, tet_left, "Quad left edge", "Tet left edge");
  }

  TEST(EdgeConformity3D2D, Segment_Tetrahedron_SharedEdge_K4)
  {
    // 1D segment should match any tetrahedron edge when parameterized
    auto seg = extractSegmentNodes<4>();
    auto tet_edge = extractTetrahedronEdgeNodes<4>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(seg, tet_edge, "Segment", "Tetrahedron edge");

    // Also test another tet edge
    auto tet_edge2 = extractTetrahedronEdgeNodes<4>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(seg, tet_edge2, "Segment", "Tetrahedron z-edge");
  }

  TEST(EdgeConformity3D2D, Segment_Wedge_AllEdgeTypes_K4)
  {
    // Segment should match all edge types on wedge
    auto seg = extractSegmentNodes<4>();

    // Triangle edge (horizontal in xy-plane)
    auto wedge_tri = extractWedgeEdgeNodes<4>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(seg, wedge_tri, "Segment", "Wedge triangle edge");

    // Vertical edge (z-direction)
    auto wedge_vert = extractWedgeEdgeNodes<4>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(seg, wedge_vert, "Segment", "Wedge vertical edge");

    // Diagonal triangle edge
    auto wedge_diag = extractWedgeEdgeNodes<4>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(seg, wedge_diag, "Segment", "Wedge diagonal edge");
  }

  TEST(EdgeConformity3D2D, AllGeometries_HighOrder_K5)
  {
    // Comprehensive test: all geometry combinations at K=5
    auto seg = extractSegmentNodes<5>();
    auto tri_edge = extractTriangleEdgeNodes<5>(0.0, 0.0, 1.0, 0.0);
    auto quad_edge = extractQuadrilateralEdgeNodes<5>(0.0, 0.0, 1.0, 0.0);
    auto tet_edge = extractTetrahedronEdgeNodes<5>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    auto wedge_edge = extractWedgeEdgeNodes<5>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);

    // 1D-2D
    compareEdgeCoordinates(seg, tri_edge, "Segment", "Triangle edge");
    compareEdgeCoordinates(seg, quad_edge, "Segment", "Quadrilateral edge");

    // 1D-3D
    compareEdgeCoordinates(seg, tet_edge, "Segment", "Tetrahedron edge");
    compareEdgeCoordinates(seg, wedge_edge, "Segment", "Wedge edge");

    // 2D-2D
    compareEdgeCoordinates(tri_edge, quad_edge, "Triangle edge", "Quadrilateral edge");

    // 2D-3D
    compareEdgeCoordinates(tri_edge, tet_edge, "Triangle edge", "Tetrahedron edge");
    compareEdgeCoordinates(tri_edge, wedge_edge, "Triangle edge", "Wedge edge");
    compareEdgeCoordinates(quad_edge, tet_edge, "Quadrilateral edge", "Tetrahedron edge");
    compareEdgeCoordinates(quad_edge, wedge_edge, "Quadrilateral edge", "Wedge edge");

    // 3D-3D
    compareEdgeCoordinates(tet_edge, wedge_edge, "Tetrahedron edge", "Wedge edge");
  }

  TEST(EdgeConformity3D2D, AllGeometries_HighOrder_K6)
  {
    // Comprehensive test: all geometry combinations at K=6
    auto seg = extractSegmentNodes<6>();
    auto tri_edge = extractTriangleEdgeNodes<6>(0.0, 0.0, 1.0, 0.0);
    auto quad_edge = extractQuadrilateralEdgeNodes<6>(0.0, 0.0, 1.0, 0.0);
    auto tet_edge = extractTetrahedronEdgeNodes<6>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    auto wedge_edge = extractWedgeEdgeNodes<6>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);

    // All pairs should be equal
    compareEdgeCoordinates(seg, tri_edge, "Segment K=6", "Triangle edge K=6");
    compareEdgeCoordinates(seg, quad_edge, "Segment K=6", "Quadrilateral edge K=6");
    compareEdgeCoordinates(seg, tet_edge, "Segment K=6", "Tetrahedron edge K=6");
    compareEdgeCoordinates(seg, wedge_edge, "Segment K=6", "Wedge edge K=6");
    compareEdgeCoordinates(tri_edge, quad_edge, "Triangle K=6", "Quadrilateral K=6");
    compareEdgeCoordinates(tri_edge, tet_edge, "Triangle K=6", "Tetrahedron K=6");
    compareEdgeCoordinates(tri_edge, wedge_edge, "Triangle K=6", "Wedge K=6");
    compareEdgeCoordinates(tet_edge, wedge_edge, "Tetrahedron K=6", "Wedge K=6");
  }

  //==========================================================================
  // Very high-order conformity tests (K = 15)
  //==========================================================================

  TEST(EdgeConformity3D2D, AllGeometries_HighOrder_K15)
  {
    // Comprehensive test: all geometry combinations at K=15
    auto seg = extractSegmentNodes<15>();
    auto tri_edge = extractTriangleEdgeNodes<15>(0.0, 0.0, 1.0, 0.0);
    auto quad_edge = extractQuadrilateralEdgeNodes<15>(0.0, 0.0, 1.0, 0.0);
    auto tet_edge = extractTetrahedronEdgeNodes<15>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    auto wedge_edge = extractWedgeEdgeNodes<15>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);

    // 1D-2D
    compareEdgeCoordinates(seg, tri_edge, "Segment K=15", "Triangle edge K=15");
    compareEdgeCoordinates(seg, quad_edge, "Segment K=15", "Quadrilateral edge K=15");

    // 1D-3D
    compareEdgeCoordinates(seg, tet_edge, "Segment K=15", "Tetrahedron edge K=15");
    compareEdgeCoordinates(seg, wedge_edge, "Segment K=15", "Wedge edge K=15");

    // 2D-2D
    compareEdgeCoordinates(tri_edge, quad_edge, "Triangle K=15", "Quadrilateral K=15");

    // 2D-3D
    compareEdgeCoordinates(tri_edge, tet_edge, "Triangle K=15", "Tetrahedron K=15");
    compareEdgeCoordinates(tri_edge, wedge_edge, "Triangle K=15", "Wedge K=15");
    compareEdgeCoordinates(quad_edge, tet_edge, "Quadrilateral K=15", "Tetrahedron K=15");
    compareEdgeCoordinates(quad_edge, wedge_edge, "Quadrilateral K=15", "Wedge K=15");

    // 3D-3D
    compareEdgeCoordinates(tet_edge, wedge_edge, "Tetrahedron K=15", "Wedge K=15");
  }

  TEST(EdgeConformity, AllGeometries_EdgeNodes_Match_K15)
  {
    const auto& gll = GLL01<15>::getNodes();
    std::vector<Real> gll_vec(gll.begin(), gll.end());

    // Segment
    auto seg = extractSegmentNodes<15>();
    compareEdgeCoordinates(seg, gll_vec, "Segment K=15", "GLL01");

    // Triangle
    auto tri_edge = extractTriangleEdgeNodes<15>(0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tri_edge, gll_vec, "Triangle edge K=15", "GLL01");

    // Quadrilateral
    auto quad_edge = extractQuadrilateralEdgeNodes<15>(0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(quad_edge, gll_vec, "Quadrilateral edge K=15", "GLL01");

    // Tetrahedron
    auto tet_edge = extractTetrahedronEdgeNodes<15>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(tet_edge, gll_vec, "Tetrahedron edge K=15", "GLL01");

    // Wedge - triangle edge
    auto wedge_tri_edge = extractWedgeEdgeNodes<15>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(wedge_tri_edge, gll_vec, "Wedge triangle edge K=15", "GLL01");

    // Wedge - vertical edge
    auto wedge_vert_edge = extractWedgeEdgeNodes<15>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    compareEdgeCoordinates(wedge_vert_edge, gll_vec, "Wedge vertical edge K=15", "GLL01");
  }

  TEST(EdgeConformity3D, Tetrahedron_Wedge_AllSharedEdges_K15)
  {
    // Test all edges that can be shared between tetrahedron and wedge base face
    // Reference tetrahedron: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    // Reference wedge base: triangle at z=0 with vertices (0,0,0), (1,0,0), (0,1,0)

    // Edge 01: (0,0,0) to (1,0,0)
    auto tet_01 = extractTetrahedronEdgeNodes<15>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    auto wedge_01 = extractWedgeEdgeNodes<15>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    compareEdgeCoordinates(tet_01, wedge_01, "Tet-Wedge edge 01 K=15", "same");

    // Edge 02: (0,0,0) to (0,1,0)
    auto tet_02 = extractTetrahedronEdgeNodes<15>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    auto wedge_02 = extractWedgeEdgeNodes<15>(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tet_02, wedge_02, "Tet-Wedge edge 02 K=15", "same");

    // Edge 12: (1,0,0) to (0,1,0)
    auto tet_12 = extractTetrahedronEdgeNodes<15>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    auto wedge_12 = extractWedgeEdgeNodes<15>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    compareEdgeCoordinates(tet_12, wedge_12, "Tet-Wedge edge 12 K=15", "same");
  }

  TEST(EdgeConformity, EdgeNodeCount_K15)
  {
    // Segment should have K+1 = 16 nodes at K=15
    auto seg = extractSegmentNodes<15>();
    EXPECT_EQ(seg.size(), 16u);

    // Triangle edge should also have 16 nodes
    auto tri_edge = extractTriangleEdgeNodes<15>(0.0, 0.0, 1.0, 0.0);
    EXPECT_EQ(tri_edge.size(), 16u);

    // Quadrilateral edge should also have 16 nodes
    auto quad_edge = extractQuadrilateralEdgeNodes<15>(0.0, 0.0, 1.0, 0.0);
    EXPECT_EQ(quad_edge.size(), 16u);

    // Tetrahedron edge should also have 16 nodes
    auto tet_edge = extractTetrahedronEdgeNodes<15>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    EXPECT_EQ(tet_edge.size(), 16u);

    // Wedge edge should also have 16 nodes
    auto wedge_edge = extractWedgeEdgeNodes<15>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    EXPECT_EQ(wedge_edge.size(), 16u);
  }
}
