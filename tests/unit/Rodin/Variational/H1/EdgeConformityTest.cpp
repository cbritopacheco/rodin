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

#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/Variational/H1/GLL.h"
#include "Rodin/Variational/H1/Fekete.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  //==========================================================================
  // Helper functions for extracting edge nodes
  //==========================================================================

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
    const Real tol = 1e-10;
    std::vector<Real> coords;
    const auto& nodes = FeketeTriangle<K>::getNodes();

    for (const auto& node : nodes)
    {
      Real px = node.x();
      Real py = node.y();

      // Check if point lies on the line segment from (x0,y0) to (x1,y1)
      Real dx = x1 - x0;
      Real dy = y1 - y0;
      Real len = std::sqrt(dx*dx + dy*dy);

      // Parameter t along the edge
      Real t;
      if (std::abs(dx) > std::abs(dy))
      {
        t = (px - x0) / dx;
      }
      else
      {
        t = (py - y0) / dy;
      }

      // Check if t is in [0,1] and point is on the line
      if (t >= -tol && t <= 1.0 + tol)
      {
        Real expected_x = x0 + t * dx;
        Real expected_y = y0 + t * dy;

        if (std::abs(px - expected_x) < tol && std::abs(py - expected_y) < tol)
        {
          coords.push_back(std::max(0.0, std::min(1.0, t)));
        }
      }
    }

    std::sort(coords.begin(), coords.end());
    return coords;
  }

  // Extract nodes on Quadrilateral edge
  template <size_t K>
  std::vector<Real> extractQuadrilateralEdgeNodes(Real x0, Real y0, Real x1, Real y1)
  {
    const Real tol = 1e-10;
    std::vector<Real> coords;
    const auto& nodes = H1Element<K, Real>::getNodes(Polytope::Type::Quadrilateral);

    for (const auto& node : nodes)
    {
      Real px = node.x();
      Real py = node.y();

      Real dx = x1 - x0;
      Real dy = y1 - y0;

      Real t;
      if (std::abs(dx) > std::abs(dy))
      {
        t = (px - x0) / dx;
      }
      else
      {
        t = (py - y0) / dy;
      }

      if (t >= -tol && t <= 1.0 + tol)
      {
        Real expected_x = x0 + t * dx;
        Real expected_y = y0 + t * dy;

        if (std::abs(px - expected_x) < tol && std::abs(py - expected_y) < tol)
        {
          coords.push_back(std::max(0.0, std::min(1.0, t)));
        }
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
    const Real tol = 1e-10;
    std::vector<Real> coords;
    const auto& nodes = FeketeTetrahedron<K>::getNodes();

    for (const auto& node : nodes)
    {
      Real px = node.x();
      Real py = node.y();
      Real pz = node.z();

      Real dx = x1 - x0;
      Real dy = y1 - y0;
      Real dz = z1 - z0;

      // Find the parameter t (use the dimension with largest delta)
      Real t;
      if (std::abs(dx) >= std::abs(dy) && std::abs(dx) >= std::abs(dz))
      {
        t = (px - x0) / dx;
      }
      else if (std::abs(dy) >= std::abs(dz))
      {
        t = (py - y0) / dy;
      }
      else
      {
        t = (pz - z0) / dz;
      }

      if (t >= -tol && t <= 1.0 + tol)
      {
        Real expected_x = x0 + t * dx;
        Real expected_y = y0 + t * dy;
        Real expected_z = z0 + t * dz;

        if (std::abs(px - expected_x) < tol &&
            std::abs(py - expected_y) < tol &&
            std::abs(pz - expected_z) < tol)
        {
          coords.push_back(std::max(0.0, std::min(1.0, t)));
        }
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
    const Real tol = 1e-10;
    std::vector<Real> coords;
    const auto& nodes = H1Element<K, Real>::getNodes(Polytope::Type::Wedge);

    for (const auto& node : nodes)
    {
      Real px = node.x();
      Real py = node.y();
      Real pz = node.z();

      Real dx = x1 - x0;
      Real dy = y1 - y0;
      Real dz = z1 - z0;

      Real t;
      if (std::abs(dx) >= std::abs(dy) && std::abs(dx) >= std::abs(dz))
      {
        t = (std::abs(dx) > tol) ? (px - x0) / dx : 0.0;
      }
      else if (std::abs(dy) >= std::abs(dz))
      {
        t = (std::abs(dy) > tol) ? (py - y0) / dy : 0.0;
      }
      else
      {
        t = (std::abs(dz) > tol) ? (pz - z0) / dz : 0.0;
      }

      if (t >= -tol && t <= 1.0 + tol)
      {
        Real expected_x = x0 + t * dx;
        Real expected_y = y0 + t * dy;
        Real expected_z = z0 + t * dz;

        if (std::abs(px - expected_x) < tol &&
            std::abs(py - expected_y) < tol &&
            std::abs(pz - expected_z) < tol)
        {
          coords.push_back(std::max(0.0, std::min(1.0, t)));
        }
      }
    }

    std::sort(coords.begin(), coords.end());
    return coords;
  }

  // Compare two vectors of 1D coordinates
  void compareEdgeCoordinates(
      const std::vector<Real>& edge1,
      const std::vector<Real>& edge2,
      const std::string& name1,
      const std::string& name2)
  {
    const Real tol = 1e-10;

    ASSERT_EQ(edge1.size(), edge2.size())
        << "Edge node counts differ: " << name1 << " has " << edge1.size()
        << " nodes, " << name2 << " has " << edge2.size() << " nodes";

    for (size_t i = 0; i < edge1.size(); ++i)
    {
      EXPECT_NEAR(edge1[i], edge2[i], tol)
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

    // Top edge: (1,1) to (0,1)
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
}
