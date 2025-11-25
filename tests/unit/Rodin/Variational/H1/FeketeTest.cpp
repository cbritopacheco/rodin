/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational/H1/Fekete.h"

using namespace Rodin;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  //==========================================================================
  // FeketeTriangle Node Count Tests
  //==========================================================================

  TEST(FeketeTriangle, NodeCount_K0)
  {
    // (K+1)(K+2)/2 = (1)(2)/2 = 1
    EXPECT_EQ(FeketeTriangle<0>::Count, 1);
    EXPECT_EQ(FeketeTriangle<0>::getNodes().size(), 1);
  }

  TEST(FeketeTriangle, NodeCount_K1)
  {
    // (K+1)(K+2)/2 = (2)(3)/2 = 3
    EXPECT_EQ(FeketeTriangle<1>::Count, 3);
    EXPECT_EQ(FeketeTriangle<1>::getNodes().size(), 3);
  }

  TEST(FeketeTriangle, NodeCount_K2)
  {
    // (K+1)(K+2)/2 = (3)(4)/2 = 6
    EXPECT_EQ(FeketeTriangle<2>::Count, 6);
    EXPECT_EQ(FeketeTriangle<2>::getNodes().size(), 6);
  }

  TEST(FeketeTriangle, NodeCount_K3)
  {
    // (K+1)(K+2)/2 = (4)(5)/2 = 10
    EXPECT_EQ(FeketeTriangle<3>::Count, 10);
    EXPECT_EQ(FeketeTriangle<3>::getNodes().size(), 10);
  }

  TEST(FeketeTriangle, NodeCount_K5)
  {
    // (K+1)(K+2)/2 = (6)(7)/2 = 21
    EXPECT_EQ(FeketeTriangle<5>::Count, 21);
    EXPECT_EQ(FeketeTriangle<5>::getNodes().size(), 21);
  }

  //==========================================================================
  // FeketeTriangle Node Location Tests
  //==========================================================================

  TEST(FeketeTriangle, NodesInReferenceTriangle_K2)
  {
    // All nodes should be in the reference triangle: x >= 0, y >= 0, x + y <= 1
    const auto& nodes = FeketeTriangle<2>::getNodes();
    for (const auto& node : nodes)
    {
      EXPECT_GE(node.x(), -1e-10);
      EXPECT_GE(node.y(), -1e-10);
      EXPECT_LE(node.x() + node.y(), 1.0 + 1e-10);
    }
  }

  TEST(FeketeTriangle, NodesInReferenceTriangle_K5)
  {
    const auto& nodes = FeketeTriangle<5>::getNodes();
    for (const auto& node : nodes)
    {
      EXPECT_GE(node.x(), -1e-10);
      EXPECT_GE(node.y(), -1e-10);
      EXPECT_LE(node.x() + node.y(), 1.0 + 1e-10);
    }
  }

  TEST(FeketeTriangle, NodesInReferenceTriangle_K10)
  {
    const auto& nodes = FeketeTriangle<10>::getNodes();
    for (const auto& node : nodes)
    {
      EXPECT_GE(node.x(), -1e-10);
      EXPECT_GE(node.y(), -1e-10);
      EXPECT_LE(node.x() + node.y(), 1.0 + 1e-10);
    }
  }

  //==========================================================================
  // FeketeTriangle Vertex Tests
  //==========================================================================

  TEST(FeketeTriangle, ContainsVertices_K1)
  {
    // K=1 should contain the three vertices
    const auto& nodes = FeketeTriangle<1>::getNodes();

    // Check that we have nodes near each vertex
    bool has_origin = false;
    bool has_x1 = false;
    bool has_y1 = false;

    for (const auto& node : nodes)
    {
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y()) < 1e-10)
        has_origin = true;
      if (std::abs(node.x() - 1.0) < 1e-10 && std::abs(node.y()) < 1e-10)
        has_x1 = true;
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y() - 1.0) < 1e-10)
        has_y1 = true;
    }

    EXPECT_TRUE(has_origin);
    EXPECT_TRUE(has_x1);
    EXPECT_TRUE(has_y1);
  }

  TEST(FeketeTriangle, ContainsVertices_K5)
  {
    const auto& nodes = FeketeTriangle<5>::getNodes();

    bool has_origin = false;
    bool has_x1 = false;
    bool has_y1 = false;

    for (const auto& node : nodes)
    {
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y()) < 1e-10)
        has_origin = true;
      if (std::abs(node.x() - 1.0) < 1e-10 && std::abs(node.y()) < 1e-10)
        has_x1 = true;
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y() - 1.0) < 1e-10)
        has_y1 = true;
    }

    EXPECT_TRUE(has_origin);
    EXPECT_TRUE(has_x1);
    EXPECT_TRUE(has_y1);
  }

  //==========================================================================
  // FeketeTetrahedron Node Count Tests
  //==========================================================================

  TEST(FeketeTetrahedron, NodeCount_K0)
  {
    // (K+1)(K+2)(K+3)/6 = (1)(2)(3)/6 = 1
    EXPECT_EQ(FeketeTetrahedron<0>::Count, 1);
    EXPECT_EQ(FeketeTetrahedron<0>::getNodes().size(), 1);
  }

  TEST(FeketeTetrahedron, NodeCount_K1)
  {
    // (K+1)(K+2)(K+3)/6 = (2)(3)(4)/6 = 4
    EXPECT_EQ(FeketeTetrahedron<1>::Count, 4);
    EXPECT_EQ(FeketeTetrahedron<1>::getNodes().size(), 4);
  }

  TEST(FeketeTetrahedron, NodeCount_K2)
  {
    // (K+1)(K+2)(K+3)/6 = (3)(4)(5)/6 = 10
    EXPECT_EQ(FeketeTetrahedron<2>::Count, 10);
    EXPECT_EQ(FeketeTetrahedron<2>::getNodes().size(), 10);
  }

  TEST(FeketeTetrahedron, NodeCount_K3)
  {
    // (K+1)(K+2)(K+3)/6 = (4)(5)(6)/6 = 20
    EXPECT_EQ(FeketeTetrahedron<3>::Count, 20);
    EXPECT_EQ(FeketeTetrahedron<3>::getNodes().size(), 20);
  }

  TEST(FeketeTetrahedron, NodeCount_K5)
  {
    // (K+1)(K+2)(K+3)/6 = (6)(7)(8)/6 = 56
    EXPECT_EQ(FeketeTetrahedron<5>::Count, 56);
    EXPECT_EQ(FeketeTetrahedron<5>::getNodes().size(), 56);
  }

  //==========================================================================
  // FeketeTetrahedron Node Location Tests
  //==========================================================================

  TEST(FeketeTetrahedron, NodesInReferenceTetrahedron_K2)
  {
    // All nodes should be in the reference tetrahedron:
    // x >= 0, y >= 0, z >= 0, x + y + z <= 1
    const auto& nodes = FeketeTetrahedron<2>::getNodes();
    for (const auto& node : nodes)
    {
      EXPECT_GE(node.x(), -1e-10);
      EXPECT_GE(node.y(), -1e-10);
      EXPECT_GE(node.z(), -1e-10);
      EXPECT_LE(node.x() + node.y() + node.z(), 1.0 + 1e-10);
    }
  }

  TEST(FeketeTetrahedron, NodesInReferenceTetrahedron_K5)
  {
    const auto& nodes = FeketeTetrahedron<5>::getNodes();
    for (const auto& node : nodes)
    {
      EXPECT_GE(node.x(), -1e-10);
      EXPECT_GE(node.y(), -1e-10);
      EXPECT_GE(node.z(), -1e-10);
      EXPECT_LE(node.x() + node.y() + node.z(), 1.0 + 1e-10);
    }
  }

  //==========================================================================
  // FeketeTetrahedron Vertex Tests
  //==========================================================================

  TEST(FeketeTetrahedron, ContainsVertices_K1)
  {
    // K=1 should contain the four vertices
    const auto& nodes = FeketeTetrahedron<1>::getNodes();

    bool has_v0 = false;  // (0,0,0)
    bool has_v1 = false;  // (1,0,0)
    bool has_v2 = false;  // (0,1,0)
    bool has_v3 = false;  // (0,0,1)

    for (const auto& node : nodes)
    {
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y()) < 1e-10 && std::abs(node.z()) < 1e-10)
        has_v0 = true;
      if (std::abs(node.x() - 1.0) < 1e-10 && std::abs(node.y()) < 1e-10 && std::abs(node.z()) < 1e-10)
        has_v1 = true;
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y() - 1.0) < 1e-10 && std::abs(node.z()) < 1e-10)
        has_v2 = true;
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y()) < 1e-10 && std::abs(node.z() - 1.0) < 1e-10)
        has_v3 = true;
    }

    EXPECT_TRUE(has_v0);
    EXPECT_TRUE(has_v1);
    EXPECT_TRUE(has_v2);
    EXPECT_TRUE(has_v3);
  }

  TEST(FeketeTetrahedron, ContainsVertices_K5)
  {
    const auto& nodes = FeketeTetrahedron<5>::getNodes();

    bool has_v0 = false;
    bool has_v1 = false;
    bool has_v2 = false;
    bool has_v3 = false;

    for (const auto& node : nodes)
    {
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y()) < 1e-10 && std::abs(node.z()) < 1e-10)
        has_v0 = true;
      if (std::abs(node.x() - 1.0) < 1e-10 && std::abs(node.y()) < 1e-10 && std::abs(node.z()) < 1e-10)
        has_v1 = true;
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y() - 1.0) < 1e-10 && std::abs(node.z()) < 1e-10)
        has_v2 = true;
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y()) < 1e-10 && std::abs(node.z() - 1.0) < 1e-10)
        has_v3 = true;
    }

    EXPECT_TRUE(has_v0);
    EXPECT_TRUE(has_v1);
    EXPECT_TRUE(has_v2);
    EXPECT_TRUE(has_v3);
  }

  //==========================================================================
  // Fekete Node Uniqueness Tests
  //==========================================================================

  TEST(FeketeTriangle, NodesAreUnique_K5)
  {
    const auto& nodes = FeketeTriangle<5>::getNodes();
    const Real tol = 1e-10;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
      for (size_t j = i + 1; j < nodes.size(); ++j)
      {
        Real dist = std::sqrt(
            (nodes[i].x() - nodes[j].x()) * (nodes[i].x() - nodes[j].x()) +
            (nodes[i].y() - nodes[j].y()) * (nodes[i].y() - nodes[j].y()));
        EXPECT_GT(dist, tol) << "Nodes " << i << " and " << j << " are too close";
      }
    }
  }

  TEST(FeketeTetrahedron, NodesAreUnique_K3)
  {
    const auto& nodes = FeketeTetrahedron<3>::getNodes();
    const Real tol = 1e-10;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
      for (size_t j = i + 1; j < nodes.size(); ++j)
      {
        Real dist = std::sqrt(
            (nodes[i].x() - nodes[j].x()) * (nodes[i].x() - nodes[j].x()) +
            (nodes[i].y() - nodes[j].y()) * (nodes[i].y() - nodes[j].y()) +
            (nodes[i].z() - nodes[j].z()) * (nodes[i].z() - nodes[j].z()));
        EXPECT_GT(dist, tol) << "Nodes " << i << " and " << j << " are too close";
      }
    }
  }

  //==========================================================================
  // Higher Order Tests (K = 6)
  //==========================================================================

  TEST(FeketeTriangle, NodeCount_K6)
  {
    // (K+1)(K+2)/2 = (7)(8)/2 = 28
    EXPECT_EQ(FeketeTriangle<6>::Count, 28);
    EXPECT_EQ(FeketeTriangle<6>::getNodes().size(), 28);
  }

  TEST(FeketeTriangle, NodesInReferenceTriangle_K6)
  {
    const auto& nodes = FeketeTriangle<6>::getNodes();
    for (const auto& node : nodes)
    {
      EXPECT_GE(node.x(), -1e-10);
      EXPECT_GE(node.y(), -1e-10);
      EXPECT_LE(node.x() + node.y(), 1.0 + 1e-10);
    }
  }

  TEST(FeketeTriangle, ContainsVertices_K6)
  {
    const auto& nodes = FeketeTriangle<6>::getNodes();

    bool has_origin = false;
    bool has_x1 = false;
    bool has_y1 = false;

    for (const auto& node : nodes)
    {
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y()) < 1e-10)
        has_origin = true;
      if (std::abs(node.x() - 1.0) < 1e-10 && std::abs(node.y()) < 1e-10)
        has_x1 = true;
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y() - 1.0) < 1e-10)
        has_y1 = true;
    }

    EXPECT_TRUE(has_origin);
    EXPECT_TRUE(has_x1);
    EXPECT_TRUE(has_y1);
  }

  TEST(FeketeTriangle, NodesAreUnique_K6)
  {
    const auto& nodes = FeketeTriangle<6>::getNodes();
    const Real tol = 1e-10;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
      for (size_t j = i + 1; j < nodes.size(); ++j)
      {
        Real dist = std::sqrt(
            (nodes[i].x() - nodes[j].x()) * (nodes[i].x() - nodes[j].x()) +
            (nodes[i].y() - nodes[j].y()) * (nodes[i].y() - nodes[j].y()));
        EXPECT_GT(dist, tol) << "Nodes " << i << " and " << j << " are too close";
      }
    }
  }

  TEST(FeketeTetrahedron, NodeCount_K6)
  {
    // (K+1)(K+2)(K+3)/6 = (7)(8)(9)/6 = 84
    EXPECT_EQ(FeketeTetrahedron<6>::Count, 84);
    EXPECT_EQ(FeketeTetrahedron<6>::getNodes().size(), 84);
  }

  TEST(FeketeTetrahedron, NodesInReferenceTetrahedron_K6)
  {
    const auto& nodes = FeketeTetrahedron<6>::getNodes();
    for (const auto& node : nodes)
    {
      EXPECT_GE(node.x(), -1e-10);
      EXPECT_GE(node.y(), -1e-10);
      EXPECT_GE(node.z(), -1e-10);
      EXPECT_LE(node.x() + node.y() + node.z(), 1.0 + 1e-10);
    }
  }

  TEST(FeketeTetrahedron, ContainsVertices_K6)
  {
    const auto& nodes = FeketeTetrahedron<6>::getNodes();

    bool has_v0 = false;
    bool has_v1 = false;
    bool has_v2 = false;
    bool has_v3 = false;

    for (const auto& node : nodes)
    {
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y()) < 1e-10 && std::abs(node.z()) < 1e-10)
        has_v0 = true;
      if (std::abs(node.x() - 1.0) < 1e-10 && std::abs(node.y()) < 1e-10 && std::abs(node.z()) < 1e-10)
        has_v1 = true;
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y() - 1.0) < 1e-10 && std::abs(node.z()) < 1e-10)
        has_v2 = true;
      if (std::abs(node.x()) < 1e-10 && std::abs(node.y()) < 1e-10 && std::abs(node.z() - 1.0) < 1e-10)
        has_v3 = true;
    }

    EXPECT_TRUE(has_v0);
    EXPECT_TRUE(has_v1);
    EXPECT_TRUE(has_v2);
    EXPECT_TRUE(has_v3);
  }

  TEST(FeketeTetrahedron, NodesAreUnique_K5)
  {
    const auto& nodes = FeketeTetrahedron<5>::getNodes();
    const Real tol = 1e-10;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
      for (size_t j = i + 1; j < nodes.size(); ++j)
      {
        Real dist = std::sqrt(
            (nodes[i].x() - nodes[j].x()) * (nodes[i].x() - nodes[j].x()) +
            (nodes[i].y() - nodes[j].y()) * (nodes[i].y() - nodes[j].y()) +
            (nodes[i].z() - nodes[j].z()) * (nodes[i].z() - nodes[j].z()));
        EXPECT_GT(dist, tol) << "Nodes " << i << " and " << j << " are too close";
      }
    }
  }
}
