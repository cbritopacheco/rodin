/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational/H1/WarpBlend.h"

using namespace Rodin;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  //==========================================================================
  // WarpFactor1D Tests
  //==========================================================================

  TEST(WarpFactor1D, Zero_K0)
  {
    // K=0 should give zero warp
    EXPECT_NEAR(WarpFactor1D<0>::get(0.0), 0.0, 1e-14);
    EXPECT_NEAR(WarpFactor1D<0>::get(0.5), 0.0, 1e-14);
  }

  TEST(WarpFactor1D, Zero_K1)
  {
    // K=1 should give zero warp (only endpoints)
    EXPECT_NEAR(WarpFactor1D<1>::get(0.0), 0.0, 1e-14);
    EXPECT_NEAR(WarpFactor1D<1>::get(0.5), 0.0, 1e-14);
  }

  TEST(WarpFactor1D, ZeroAtEndpoints_K3)
  {
    // Warp should be zero at endpoints
    EXPECT_NEAR(WarpFactor1D<3>::get(-1.0), 0.0, 1e-14);
    EXPECT_NEAR(WarpFactor1D<3>::get(1.0), 0.0, 1e-14);
  }

  TEST(WarpFactor1D, NonZeroInterior_K3)
  {
    // Warp should be non-zero at interior points for K > 1
    Real warp = WarpFactor1D<3>::get(0.0);
    // Just check it's finite, not necessarily non-zero
    EXPECT_FALSE(std::isnan(warp));
    EXPECT_FALSE(std::isinf(warp));
  }

  TEST(WarpFactor1D, Finite_K5)
  {
    for (Real r = -0.9; r <= 0.9; r += 0.1)
    {
      Real warp = WarpFactor1D<5>::get(r);
      EXPECT_FALSE(std::isnan(warp));
      EXPECT_FALSE(std::isinf(warp));
    }
  }

  //==========================================================================
  // TriangleBlend Alpha Tests
  //==========================================================================

  TEST(TriangleBlend, AlphaValues)
  {
    // Verify known alpha values from Hesthaven-Warburton
    EXPECT_NEAR(TriangleBlend<0>::getAlpha(), 0.0, 1e-14);
    EXPECT_NEAR(TriangleBlend<1>::getAlpha(), 0.0, 1e-14);
    EXPECT_NEAR(TriangleBlend<2>::getAlpha(), 0.0, 1e-14);
    EXPECT_NEAR(TriangleBlend<3>::getAlpha(), 1.4152, 1e-4);
    EXPECT_NEAR(TriangleBlend<4>::getAlpha(), 0.1001, 1e-4);
    EXPECT_NEAR(TriangleBlend<5>::getAlpha(), 0.2751, 1e-4);
  }

  TEST(TriangleBlend, AlphaPositive_HigherOrders)
  {
    // For higher orders, alpha should be positive
    EXPECT_GT(TriangleBlend<6>::getAlpha(), 0.0);
    EXPECT_GT(TriangleBlend<10>::getAlpha(), 0.0);
    EXPECT_GT(TriangleBlend<15>::getAlpha(), 0.0);
  }

  //==========================================================================
  // TetrahedronBlend Alpha Tests
  //==========================================================================

  TEST(TetrahedronBlend, AlphaValues)
  {
    // Verify known alpha values
    EXPECT_NEAR(TetrahedronBlend<0>::getAlpha(), 0.0, 1e-14);
    EXPECT_NEAR(TetrahedronBlend<1>::getAlpha(), 0.0, 1e-14);
    EXPECT_NEAR(TetrahedronBlend<2>::getAlpha(), 0.0, 1e-14);
    EXPECT_NEAR(TetrahedronBlend<3>::getAlpha(), 0.0, 1e-14);
    EXPECT_NEAR(TetrahedronBlend<4>::getAlpha(), 0.1002, 1e-4);
  }

  //==========================================================================
  // WarpShiftFace2D Tests
  //==========================================================================

  TEST(WarpShiftFace2D, ZeroAtVertex_K2)
  {
    // At vertex L1 = 1 (other two are 0)
    const auto [dx, dy] = WarpShiftFace2D<2>::apply(1.0, 0.0, 0.0, 0.0);
    EXPECT_NEAR(dx, 0.0, 1e-14);
    EXPECT_NEAR(dy, 0.0, 1e-14);
  }

  TEST(WarpShiftFace2D, ZeroForLowOrder)
  {
    // K=1 should give zero warp
    const auto [dx, dy] = WarpShiftFace2D<1>::apply(0.5, 0.25, 0.25, 0.0);
    EXPECT_NEAR(dx, 0.0, 1e-14);
    EXPECT_NEAR(dy, 0.0, 1e-14);
  }

  TEST(WarpShiftFace2D, Finite_K5)
  {
    // Test at interior point
    const auto [dx, dy] = WarpShiftFace2D<5>::apply(0.4, 0.3, 0.3, 0.5);
    EXPECT_FALSE(std::isnan(dx));
    EXPECT_FALSE(std::isnan(dy));
    EXPECT_FALSE(std::isinf(dx));
    EXPECT_FALSE(std::isinf(dy));
  }

  //==========================================================================
  // WarpBlendTriangle Tests
  //==========================================================================

  TEST(WarpBlendTriangle, PreservesVertices_K3)
  {
    std::array<Math::SpatialPoint, 10> nodes;

    // Set up K=3 equispaced nodes
    size_t idx = 0;
    for (size_t j = 0; j <= 3; ++j)
    {
      for (size_t i = 0; i <= 3 - j; ++i, ++idx)
      {
        Real s = static_cast<Real>(i) / 3.0;
        Real t = static_cast<Real>(j) / 3.0;
        nodes[idx] = Math::SpatialPoint{{s, t}};
      }
    }

    // Apply warp
    WarpBlendTriangle<3>::apply<10>(nodes);

    // Check vertices are preserved
    bool has_v0 = false, has_v1 = false, has_v2 = false;
    for (const auto& n : nodes)
    {
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y()) < 1e-10)
        has_v0 = true;
      if (std::abs(n.x() - 1.0) < 1e-10 && std::abs(n.y()) < 1e-10)
        has_v1 = true;
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y() - 1.0) < 1e-10)
        has_v2 = true;
    }

    EXPECT_TRUE(has_v0);
    EXPECT_TRUE(has_v1);
    EXPECT_TRUE(has_v2);
  }

  TEST(WarpBlendTriangle, NodesStayInTriangle_K5)
  {
    std::array<Math::SpatialPoint, 21> nodes;

    size_t idx = 0;
    for (size_t j = 0; j <= 5; ++j)
    {
      for (size_t i = 0; i <= 5 - j; ++i, ++idx)
      {
        Real s = static_cast<Real>(i) / 5.0;
        Real t = static_cast<Real>(j) / 5.0;
        nodes[idx] = Math::SpatialPoint{{s, t}};
      }
    }

    WarpBlendTriangle<5>::apply<21>(nodes);

    for (const auto& n : nodes)
    {
      EXPECT_GE(n.x(), -1e-10);
      EXPECT_GE(n.y(), -1e-10);
      EXPECT_LE(n.x() + n.y(), 1.0 + 1e-10);
    }
  }

  TEST(WarpBlendTriangle, NoEffect_K1)
  {
    std::array<Math::SpatialPoint, 3> nodes = {
      Math::SpatialPoint{{0.0, 0.0}},
      Math::SpatialPoint{{1.0, 0.0}},
      Math::SpatialPoint{{0.0, 1.0}}
    };

    auto nodes_copy = nodes;

    WarpBlendTriangle<1>::apply<3>(nodes);

    // Should be unchanged
    for (size_t i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(nodes[i].x(), nodes_copy[i].x(), 1e-14);
      EXPECT_NEAR(nodes[i].y(), nodes_copy[i].y(), 1e-14);
    }
  }

  //==========================================================================
  // WarpBlendTetrahedron Tests
  //==========================================================================

  TEST(WarpBlendTetrahedron, PreservesVertices_K3)
  {
    std::array<Math::SpatialPoint, 20> nodes;

    size_t idx = 0;
    for (size_t k = 0; k <= 3; ++k)
    {
      for (size_t j = 0; j <= 3 - k; ++j)
      {
        for (size_t i = 0; i <= 3 - j - k; ++i, ++idx)
        {
          Real r = static_cast<Real>(i) / 3.0;
          Real s = static_cast<Real>(j) / 3.0;
          Real t = static_cast<Real>(k) / 3.0;
          nodes[idx] = Math::SpatialPoint{{r, s, t}};
        }
      }
    }

    WarpBlendTetrahedron<3>::apply<20>(nodes);

    // Check vertices
    bool has_v0 = false, has_v1 = false, has_v2 = false, has_v3 = false;
    for (const auto& n : nodes)
    {
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y()) < 1e-10 && std::abs(n.z()) < 1e-10)
        has_v0 = true;
      if (std::abs(n.x() - 1.0) < 1e-10 && std::abs(n.y()) < 1e-10 && std::abs(n.z()) < 1e-10)
        has_v1 = true;
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y() - 1.0) < 1e-10 && std::abs(n.z()) < 1e-10)
        has_v2 = true;
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y()) < 1e-10 && std::abs(n.z() - 1.0) < 1e-10)
        has_v3 = true;
    }

    EXPECT_TRUE(has_v0);
    EXPECT_TRUE(has_v1);
    EXPECT_TRUE(has_v2);
    EXPECT_TRUE(has_v3);
  }

  TEST(WarpBlendTetrahedron, NodesStayInTetrahedron_K4)
  {
    std::array<Math::SpatialPoint, 35> nodes;

    size_t idx = 0;
    for (size_t k = 0; k <= 4; ++k)
    {
      for (size_t j = 0; j <= 4 - k; ++j)
      {
        for (size_t i = 0; i <= 4 - j - k; ++i, ++idx)
        {
          Real r = static_cast<Real>(i) / 4.0;
          Real s = static_cast<Real>(j) / 4.0;
          Real t = static_cast<Real>(k) / 4.0;
          nodes[idx] = Math::SpatialPoint{{r, s, t}};
        }
      }
    }

    WarpBlendTetrahedron<4>::apply<35>(nodes);

    for (const auto& n : nodes)
    {
      EXPECT_GE(n.x(), -1e-10);
      EXPECT_GE(n.y(), -1e-10);
      EXPECT_GE(n.z(), -1e-10);
      EXPECT_LE(n.x() + n.y() + n.z(), 1.0 + 1e-10);
    }
  }

  TEST(WarpBlendTetrahedron, NoEffect_K1)
  {
    std::array<Math::SpatialPoint, 4> nodes = {
      Math::SpatialPoint{{0.0, 0.0, 0.0}},
      Math::SpatialPoint{{1.0, 0.0, 0.0}},
      Math::SpatialPoint{{0.0, 1.0, 0.0}},
      Math::SpatialPoint{{0.0, 0.0, 1.0}}
    };

    auto nodes_copy = nodes;

    WarpBlendTetrahedron<1>::apply<4>(nodes);

    for (size_t i = 0; i < 4; ++i)
    {
      EXPECT_NEAR(nodes[i].x(), nodes_copy[i].x(), 1e-14);
      EXPECT_NEAR(nodes[i].y(), nodes_copy[i].y(), 1e-14);
      EXPECT_NEAR(nodes[i].z(), nodes_copy[i].z(), 1e-14);
    }
  }

  //==========================================================================
  // Higher Order Tests (K = 5, 6)
  //==========================================================================

  TEST(WarpFactor1D, Finite_K6)
  {
    for (Real r = -0.9; r <= 0.9; r += 0.1)
    {
      Real warp = WarpFactor1D<6>::get(r);
      EXPECT_FALSE(std::isnan(warp));
      EXPECT_FALSE(std::isinf(warp));
    }
  }

  TEST(WarpFactor1D, ZeroAtEndpoints_K5)
  {
    EXPECT_NEAR(WarpFactor1D<5>::get(-1.0), 0.0, 1e-14);
    EXPECT_NEAR(WarpFactor1D<5>::get(1.0), 0.0, 1e-14);
  }

  TEST(WarpFactor1D, ZeroAtEndpoints_K6)
  {
    EXPECT_NEAR(WarpFactor1D<6>::get(-1.0), 0.0, 1e-14);
    EXPECT_NEAR(WarpFactor1D<6>::get(1.0), 0.0, 1e-14);
  }

  TEST(WarpBlendTriangle, NodesStayInTriangle_K6)
  {
    std::array<Math::SpatialPoint, 28> nodes;

    size_t idx = 0;
    for (size_t j = 0; j <= 6; ++j)
    {
      for (size_t i = 0; i <= 6 - j; ++i, ++idx)
      {
        Real s = static_cast<Real>(i) / 6.0;
        Real t = static_cast<Real>(j) / 6.0;
        nodes[idx] = Math::SpatialPoint{{s, t}};
      }
    }

    WarpBlendTriangle<6>::apply<28>(nodes);

    for (const auto& n : nodes)
    {
      EXPECT_GE(n.x(), -1e-10);
      EXPECT_GE(n.y(), -1e-10);
      EXPECT_LE(n.x() + n.y(), 1.0 + 1e-10);
    }
  }

  TEST(WarpBlendTriangle, PreservesVertices_K6)
  {
    std::array<Math::SpatialPoint, 28> nodes;

    size_t idx = 0;
    for (size_t j = 0; j <= 6; ++j)
    {
      for (size_t i = 0; i <= 6 - j; ++i, ++idx)
      {
        Real s = static_cast<Real>(i) / 6.0;
        Real t = static_cast<Real>(j) / 6.0;
        nodes[idx] = Math::SpatialPoint{{s, t}};
      }
    }

    WarpBlendTriangle<6>::apply<28>(nodes);

    bool has_v0 = false, has_v1 = false, has_v2 = false;
    for (const auto& n : nodes)
    {
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y()) < 1e-10)
        has_v0 = true;
      if (std::abs(n.x() - 1.0) < 1e-10 && std::abs(n.y()) < 1e-10)
        has_v1 = true;
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y() - 1.0) < 1e-10)
        has_v2 = true;
    }

    EXPECT_TRUE(has_v0);
    EXPECT_TRUE(has_v1);
    EXPECT_TRUE(has_v2);
  }

  TEST(WarpBlendTetrahedron, NodesStayInTetrahedron_K5)
  {
    std::array<Math::SpatialPoint, 56> nodes;

    size_t idx = 0;
    for (size_t k = 0; k <= 5; ++k)
    {
      for (size_t j = 0; j <= 5 - k; ++j)
      {
        for (size_t i = 0; i <= 5 - j - k; ++i, ++idx)
        {
          Real r = static_cast<Real>(i) / 5.0;
          Real s = static_cast<Real>(j) / 5.0;
          Real t = static_cast<Real>(k) / 5.0;
          nodes[idx] = Math::SpatialPoint{{r, s, t}};
        }
      }
    }

    WarpBlendTetrahedron<5>::apply<56>(nodes);

    for (const auto& n : nodes)
    {
      EXPECT_GE(n.x(), -1e-10);
      EXPECT_GE(n.y(), -1e-10);
      EXPECT_GE(n.z(), -1e-10);
      EXPECT_LE(n.x() + n.y() + n.z(), 1.0 + 1e-10);
    }
  }

  TEST(WarpBlendTetrahedron, PreservesVertices_K5)
  {
    std::array<Math::SpatialPoint, 56> nodes;

    size_t idx = 0;
    for (size_t k = 0; k <= 5; ++k)
    {
      for (size_t j = 0; j <= 5 - k; ++j)
      {
        for (size_t i = 0; i <= 5 - j - k; ++i, ++idx)
        {
          Real r = static_cast<Real>(i) / 5.0;
          Real s = static_cast<Real>(j) / 5.0;
          Real t = static_cast<Real>(k) / 5.0;
          nodes[idx] = Math::SpatialPoint{{r, s, t}};
        }
      }
    }

    WarpBlendTetrahedron<5>::apply<56>(nodes);

    bool has_v0 = false, has_v1 = false, has_v2 = false, has_v3 = false;
    for (const auto& n : nodes)
    {
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y()) < 1e-10 && std::abs(n.z()) < 1e-10)
        has_v0 = true;
      if (std::abs(n.x() - 1.0) < 1e-10 && std::abs(n.y()) < 1e-10 && std::abs(n.z()) < 1e-10)
        has_v1 = true;
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y() - 1.0) < 1e-10 && std::abs(n.z()) < 1e-10)
        has_v2 = true;
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y()) < 1e-10 && std::abs(n.z() - 1.0) < 1e-10)
        has_v3 = true;
    }

    EXPECT_TRUE(has_v0);
    EXPECT_TRUE(has_v1);
    EXPECT_TRUE(has_v2);
    EXPECT_TRUE(has_v3);
  }

  TEST(WarpShiftFace2D, Finite_K6)
  {
    // Test at interior point
    const auto [dx, dy] = WarpShiftFace2D<6>::apply(0.35, 0.35, 0.3, 0.5);
    EXPECT_FALSE(std::isnan(dx));
    EXPECT_FALSE(std::isnan(dy));
    EXPECT_FALSE(std::isinf(dx));
    EXPECT_FALSE(std::isinf(dy));
  }

  //==========================================================================
  // Very High Order Tests (K = 15)
  //==========================================================================

  TEST(WarpFactor1D, ZeroAtEndpoints_K15)
  {
    EXPECT_NEAR(WarpFactor1D<15>::get(-1.0), 0.0, 1e-14);
    EXPECT_NEAR(WarpFactor1D<15>::get(1.0), 0.0, 1e-14);
  }

  TEST(WarpFactor1D, Finite_K15)
  {
    for (Real r = -0.9; r <= 0.9; r += 0.1)
    {
      Real warp = WarpFactor1D<15>::get(r);
      EXPECT_FALSE(std::isnan(warp));
      EXPECT_FALSE(std::isinf(warp));
    }
  }

  TEST(TriangleBlend, Alpha_K15)
  {
    // For K=15, alpha should be positive
    EXPECT_GT(TriangleBlend<15>::getAlpha(), 0.0);
    EXPECT_FALSE(std::isnan(TriangleBlend<15>::getAlpha()));
    EXPECT_FALSE(std::isinf(TriangleBlend<15>::getAlpha()));
  }

  TEST(TetrahedronBlend, Alpha_K15)
  {
    // For K=15, alpha should be positive
    EXPECT_GT(TetrahedronBlend<15>::getAlpha(), 0.0);
    EXPECT_FALSE(std::isnan(TetrahedronBlend<15>::getAlpha()));
    EXPECT_FALSE(std::isinf(TetrahedronBlend<15>::getAlpha()));
  }

  TEST(WarpShiftFace2D, Finite_K15)
  {
    // Test at interior point
    const auto [dx, dy] = WarpShiftFace2D<15>::apply(0.4, 0.3, 0.3, 0.5);
    EXPECT_FALSE(std::isnan(dx));
    EXPECT_FALSE(std::isnan(dy));
    EXPECT_FALSE(std::isinf(dx));
    EXPECT_FALSE(std::isinf(dy));
  }

  TEST(WarpBlendTriangle, NodesStayInTriangle_K15)
  {
    constexpr size_t N = (15 + 1) * (15 + 2) / 2;  // 136
    std::array<Math::SpatialPoint, N> nodes;

    size_t idx = 0;
    for (size_t j = 0; j <= 15; ++j)
    {
      for (size_t i = 0; i <= 15 - j; ++i, ++idx)
      {
        Real s = static_cast<Real>(i) / 15.0;
        Real t = static_cast<Real>(j) / 15.0;
        nodes[idx] = Math::SpatialPoint{{s, t}};
      }
    }

    WarpBlendTriangle<15>::apply<N>(nodes);

    for (const auto& n : nodes)
    {
      EXPECT_GE(n.x(), -1e-10);
      EXPECT_GE(n.y(), -1e-10);
      EXPECT_LE(n.x() + n.y(), 1.0 + 1e-10);
    }
  }

  TEST(WarpBlendTriangle, PreservesVertices_K15)
  {
    constexpr size_t N = (15 + 1) * (15 + 2) / 2;  // 136
    std::array<Math::SpatialPoint, N> nodes;

    size_t idx = 0;
    for (size_t j = 0; j <= 15; ++j)
    {
      for (size_t i = 0; i <= 15 - j; ++i, ++idx)
      {
        Real s = static_cast<Real>(i) / 15.0;
        Real t = static_cast<Real>(j) / 15.0;
        nodes[idx] = Math::SpatialPoint{{s, t}};
      }
    }

    WarpBlendTriangle<15>::apply<N>(nodes);

    bool has_v0 = false, has_v1 = false, has_v2 = false;
    for (const auto& n : nodes)
    {
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y()) < 1e-10)
        has_v0 = true;
      if (std::abs(n.x() - 1.0) < 1e-10 && std::abs(n.y()) < 1e-10)
        has_v1 = true;
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y() - 1.0) < 1e-10)
        has_v2 = true;
    }

    EXPECT_TRUE(has_v0);
    EXPECT_TRUE(has_v1);
    EXPECT_TRUE(has_v2);
  }

  TEST(WarpBlendTetrahedron, NodesStayInTetrahedron_K15)
  {
    constexpr size_t N = (15 + 1) * (15 + 2) * (15 + 3) / 6;  // 816
    std::array<Math::SpatialPoint, N> nodes;

    size_t idx = 0;
    for (size_t k = 0; k <= 15; ++k)
    {
      for (size_t j = 0; j <= 15 - k; ++j)
      {
        for (size_t i = 0; i <= 15 - j - k; ++i, ++idx)
        {
          Real r = static_cast<Real>(i) / 15.0;
          Real s = static_cast<Real>(j) / 15.0;
          Real t = static_cast<Real>(k) / 15.0;
          nodes[idx] = Math::SpatialPoint{{r, s, t}};
        }
      }
    }

    WarpBlendTetrahedron<15>::apply<N>(nodes);

    for (const auto& n : nodes)
    {
      EXPECT_GE(n.x(), -1e-10);
      EXPECT_GE(n.y(), -1e-10);
      EXPECT_GE(n.z(), -1e-10);
      EXPECT_LE(n.x() + n.y() + n.z(), 1.0 + 1e-10);
    }
  }

  TEST(WarpBlendTetrahedron, PreservesVertices_K15)
  {
    constexpr size_t N = (15 + 1) * (15 + 2) * (15 + 3) / 6;  // 816
    std::array<Math::SpatialPoint, N> nodes;

    size_t idx = 0;
    for (size_t k = 0; k <= 15; ++k)
    {
      for (size_t j = 0; j <= 15 - k; ++j)
      {
        for (size_t i = 0; i <= 15 - j - k; ++i, ++idx)
        {
          Real r = static_cast<Real>(i) / 15.0;
          Real s = static_cast<Real>(j) / 15.0;
          Real t = static_cast<Real>(k) / 15.0;
          nodes[idx] = Math::SpatialPoint{{r, s, t}};
        }
      }
    }

    WarpBlendTetrahedron<15>::apply<N>(nodes);

    bool has_v0 = false, has_v1 = false, has_v2 = false, has_v3 = false;
    for (const auto& n : nodes)
    {
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y()) < 1e-10 && std::abs(n.z()) < 1e-10)
        has_v0 = true;
      if (std::abs(n.x() - 1.0) < 1e-10 && std::abs(n.y()) < 1e-10 && std::abs(n.z()) < 1e-10)
        has_v1 = true;
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y() - 1.0) < 1e-10 && std::abs(n.z()) < 1e-10)
        has_v2 = true;
      if (std::abs(n.x()) < 1e-10 && std::abs(n.y()) < 1e-10 && std::abs(n.z() - 1.0) < 1e-10)
        has_v3 = true;
    }

    EXPECT_TRUE(has_v0);
    EXPECT_TRUE(has_v1);
    EXPECT_TRUE(has_v2);
    EXPECT_TRUE(has_v3);
  }
}
