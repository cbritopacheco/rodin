/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational/H1/LagrangeBasis.h"

using namespace Rodin;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  //==========================================================================
  // LagrangeBasis1D Tests
  //==========================================================================

  TEST(LagrangeBasis1D, LagrangeProperty_K2)
  {
    // On equispaced nodes [0, 0.5, 1]
    std::array<Real, 3> nodes = {0.0, 0.5, 1.0};

    for (size_t i = 0; i < 3; ++i)
    {
      for (size_t j = 0; j < 3; ++j)
      {
        Real val = LagrangeBasis1D<2>::getBasis(i, nodes[j], nodes);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(val, expected, 1e-14);
      }
    }
  }

  TEST(LagrangeBasis1D, PartitionOfUnity_K2)
  {
    std::array<Real, 3> nodes = {0.0, 0.5, 1.0};

    for (Real x = 0.0; x <= 1.0; x += 0.1)
    {
      Real sum = 0.0;
      for (size_t i = 0; i < 3; ++i)
        sum += LagrangeBasis1D<2>::getBasis(i, x, nodes);
      EXPECT_NEAR(sum, 1.0, 1e-14);
    }
  }

  TEST(LagrangeBasis1D, DerivativeSum_K2)
  {
    // Sum of derivatives should be 0 (derivative of constant 1)
    std::array<Real, 3> nodes = {0.0, 0.5, 1.0};

    for (Real x = 0.1; x <= 0.9; x += 0.1)
    {
      Real sum = 0.0;
      for (size_t i = 0; i < 3; ++i)
        sum += LagrangeBasis1D<2>::getDerivative(i, x, nodes);
      EXPECT_NEAR(sum, 0.0, 1e-13);
    }
  }

  //==========================================================================
  // LagrangeBasisPoint Tests
  //==========================================================================

  TEST(LagrangeBasisPoint, Basis_K0)
  {
    EXPECT_NEAR(LagrangeBasisPoint<0>::getBasis(), 1.0, 1e-14);
  }

  TEST(LagrangeBasisPoint, Derivative_K0)
  {
    EXPECT_NEAR(LagrangeBasisPoint<0>::getDerivative(), 0.0, 1e-14);
  }

  //==========================================================================
  // LagrangeBasisSegment Tests (using GLL01 nodes)
  //==========================================================================

  TEST(LagrangeBasisSegment, LagrangeProperty_K2)
  {
    // GLL01<2> nodes are at 0, 0.5, 1
    const auto& nodes = GLL01<2>::getNodes();

    for (size_t i = 0; i < 3; ++i)
    {
      for (size_t j = 0; j < 3; ++j)
      {
        Real val = LagrangeBasisSegment<2>::getBasis(i, nodes[j]);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(val, expected, 1e-14);
      }
    }
  }

  TEST(LagrangeBasisSegment, LagrangeProperty_K3)
  {
    const auto& nodes = GLL01<3>::getNodes();

    for (size_t i = 0; i <= 3; ++i)
    {
      for (size_t j = 0; j <= 3; ++j)
      {
        Real val = LagrangeBasisSegment<3>::getBasis(i, nodes[j]);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(val, expected, 1e-14);
      }
    }
  }

  TEST(LagrangeBasisSegment, PartitionOfUnity_K3)
  {
    for (Real x = 0.0; x <= 1.0; x += 0.05)
    {
      Real sum = 0.0;
      for (size_t i = 0; i <= 3; ++i)
        sum += LagrangeBasisSegment<3>::getBasis(i, x);
      EXPECT_NEAR(sum, 1.0, 1e-14);
    }
  }

  TEST(LagrangeBasisSegment, DerivativeSum_K3)
  {
    for (Real x = 0.1; x <= 0.9; x += 0.1)
    {
      Real sum = 0.0;
      for (size_t i = 0; i <= 3; ++i)
        sum += LagrangeBasisSegment<3>::getDerivative(i, x);
      EXPECT_NEAR(sum, 0.0, 1e-13);
    }
  }

  //==========================================================================
  // LagrangeBasisTriangle Tests
  //==========================================================================

  TEST(LagrangeBasisTriangle, LagrangeProperty_K1)
  {
    // K=1 has 3 nodes: (0,0), (1,0), (0,1)
    // Corresponding to (i,j) = (0,0), (1,0), (0,1)
    std::vector<std::pair<size_t, size_t>> nodes = {{0,0}, {1,0}, {0,1}};
    std::vector<std::pair<Real, Real>> coords = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};

    for (size_t ni = 0; ni < 3; ++ni)
    {
      auto [i, j] = nodes[ni];
      for (size_t nj = 0; nj < 3; ++nj)
      {
        auto [x, y] = coords[nj];
        Real val = LagrangeBasisTriangle<1>::getBasis(i, j, x, y);
        Real expected = (ni == nj) ? 1.0 : 0.0;
        EXPECT_NEAR(val, expected, 1e-14) << "ni=" << ni << " nj=" << nj;
      }
    }
  }

  TEST(LagrangeBasisTriangle, PartitionOfUnity_K2)
  {
    // Test at interior point
    Real x = 0.3, y = 0.3;
    Real sum = 0.0;

    // Sum over all (i,j) with i+j <= K
    for (size_t j = 0; j <= 2; ++j)
    {
      for (size_t i = 0; i <= 2 - j; ++i)
      {
        sum += LagrangeBasisTriangle<2>::getBasis(i, j, x, y);
      }
    }
    EXPECT_NEAR(sum, 1.0, 1e-14);
  }

  TEST(LagrangeBasisTriangle, PartitionOfUnity_K3)
  {
    // Test at several points
    std::vector<std::pair<Real, Real>> test_points = {
      {0.25, 0.25}, {0.1, 0.1}, {0.5, 0.3}, {0.1, 0.8}
    };

    for (auto [x, y] : test_points)
    {
      if (x + y > 1.0) continue;  // Skip outside triangle

      Real sum = 0.0;
      for (size_t j = 0; j <= 3; ++j)
      {
        for (size_t i = 0; i <= 3 - j; ++i)
        {
          sum += LagrangeBasisTriangle<3>::getBasis(i, j, x, y);
        }
      }
      EXPECT_NEAR(sum, 1.0, 1e-14) << "x=" << x << " y=" << y;
    }
  }

  TEST(LagrangeBasisTriangle, DerivativeSum_K2)
  {
    Real x = 0.3, y = 0.3;
    Real sum_dx = 0.0;
    Real sum_dy = 0.0;

    for (size_t j = 0; j <= 2; ++j)
    {
      for (size_t i = 0; i <= 2 - j; ++i)
      {
        sum_dx += LagrangeBasisTriangle<2>::getDerivative(i, j, 0, x, y);
        sum_dy += LagrangeBasisTriangle<2>::getDerivative(i, j, 1, x, y);
      }
    }
    EXPECT_NEAR(sum_dx, 0.0, 1e-13);
    EXPECT_NEAR(sum_dy, 0.0, 1e-13);
  }

  //==========================================================================
  // LagrangeBasisTetrahedron Tests
  //==========================================================================

  TEST(LagrangeBasisTetrahedron, LagrangeProperty_K1)
  {
    // K=1 has 4 nodes: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    std::vector<std::tuple<size_t, size_t, size_t>> nodes = {
      {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}
    };
    std::vector<std::tuple<Real, Real, Real>> coords = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
    };

    for (size_t ni = 0; ni < 4; ++ni)
    {
      auto [i, j, k] = nodes[ni];
      for (size_t nj = 0; nj < 4; ++nj)
      {
        auto [x, y, z] = coords[nj];
        Real val = LagrangeBasisTetrahedron<1>::getBasis(i, j, k, x, y, z);
        Real expected = (ni == nj) ? 1.0 : 0.0;
        EXPECT_NEAR(val, expected, 1e-14) << "ni=" << ni << " nj=" << nj;
      }
    }
  }

  TEST(LagrangeBasisTetrahedron, PartitionOfUnity_K2)
  {
    Real x = 0.2, y = 0.2, z = 0.2;
    Real sum = 0.0;

    for (size_t k = 0; k <= 2; ++k)
    {
      for (size_t j = 0; j <= 2 - k; ++j)
      {
        for (size_t i = 0; i <= 2 - j - k; ++i)
        {
          sum += LagrangeBasisTetrahedron<2>::getBasis(i, j, k, x, y, z);
        }
      }
    }
    EXPECT_NEAR(sum, 1.0, 1e-14);
  }

  TEST(LagrangeBasisTetrahedron, DerivativeSum_K2)
  {
    Real x = 0.2, y = 0.2, z = 0.2;
    Real sum_dx = 0.0, sum_dy = 0.0, sum_dz = 0.0;

    for (size_t k = 0; k <= 2; ++k)
    {
      for (size_t j = 0; j <= 2 - k; ++j)
      {
        for (size_t i = 0; i <= 2 - j - k; ++i)
        {
          sum_dx += LagrangeBasisTetrahedron<2>::getDerivative(i, j, k, 0, x, y, z);
          sum_dy += LagrangeBasisTetrahedron<2>::getDerivative(i, j, k, 1, x, y, z);
          sum_dz += LagrangeBasisTetrahedron<2>::getDerivative(i, j, k, 2, x, y, z);
        }
      }
    }
    EXPECT_NEAR(sum_dx, 0.0, 1e-12);
    EXPECT_NEAR(sum_dy, 0.0, 1e-12);
    EXPECT_NEAR(sum_dz, 0.0, 1e-12);
  }

  //==========================================================================
  // LagrangeBasisQuadrilateral Tests
  //==========================================================================

  TEST(LagrangeBasisQuadrilateral, LagrangeProperty_K2)
  {
    const auto& nodes = GLL01<2>::getNodes();

    for (size_t i = 0; i <= 2; ++i)
    {
      for (size_t j = 0; j <= 2; ++j)
      {
        for (size_t ii = 0; ii <= 2; ++ii)
        {
          for (size_t jj = 0; jj <= 2; ++jj)
          {
            Real val = LagrangeBasisQuadrilateral<2>::getBasis(i, j, nodes[ii], nodes[jj]);
            Real expected = (i == ii && j == jj) ? 1.0 : 0.0;
            EXPECT_NEAR(val, expected, 1e-14);
          }
        }
      }
    }
  }

  TEST(LagrangeBasisQuadrilateral, PartitionOfUnity_K3)
  {
    for (Real x = 0.1; x <= 0.9; x += 0.2)
    {
      for (Real y = 0.1; y <= 0.9; y += 0.2)
      {
        Real sum = 0.0;
        for (size_t i = 0; i <= 3; ++i)
        {
          for (size_t j = 0; j <= 3; ++j)
          {
            sum += LagrangeBasisQuadrilateral<3>::getBasis(i, j, x, y);
          }
        }
        EXPECT_NEAR(sum, 1.0, 1e-14);
      }
    }
  }

  TEST(LagrangeBasisQuadrilateral, DerivativeSum_K2)
  {
    Real x = 0.3, y = 0.7;
    Real sum_dx = 0.0, sum_dy = 0.0;

    for (size_t i = 0; i <= 2; ++i)
    {
      for (size_t j = 0; j <= 2; ++j)
      {
        sum_dx += LagrangeBasisQuadrilateral<2>::getDerivative(i, j, 0, x, y);
        sum_dy += LagrangeBasisQuadrilateral<2>::getDerivative(i, j, 1, x, y);
      }
    }
    EXPECT_NEAR(sum_dx, 0.0, 1e-13);
    EXPECT_NEAR(sum_dy, 0.0, 1e-13);
  }

  //==========================================================================
  // LagrangeBasisWedge Tests
  //==========================================================================

  TEST(LagrangeBasisWedge, PartitionOfUnity_K2)
  {
    Real x = 0.2, y = 0.2, z = 0.5;  // Point inside wedge
    Real sum = 0.0;

    for (size_t k = 0; k <= 2; ++k)  // z direction
    {
      for (size_t j = 0; j <= 2; ++j)  // triangle j
      {
        for (size_t i = 0; i <= 2 - j; ++i)  // triangle i
        {
          sum += LagrangeBasisWedge<2>::getBasis(i, j, k, x, y, z);
        }
      }
    }
    EXPECT_NEAR(sum, 1.0, 1e-14);
  }

  TEST(LagrangeBasisWedge, DerivativeSum_K2)
  {
    Real x = 0.2, y = 0.2, z = 0.5;
    Real sum_dx = 0.0, sum_dy = 0.0, sum_dz = 0.0;

    for (size_t k = 0; k <= 2; ++k)
    {
      for (size_t j = 0; j <= 2; ++j)
      {
        for (size_t i = 0; i <= 2 - j; ++i)
        {
          sum_dx += LagrangeBasisWedge<2>::getDerivative(i, j, k, 0, x, y, z);
          sum_dy += LagrangeBasisWedge<2>::getDerivative(i, j, k, 1, x, y, z);
          sum_dz += LagrangeBasisWedge<2>::getDerivative(i, j, k, 2, x, y, z);
        }
      }
    }
    EXPECT_NEAR(sum_dx, 0.0, 1e-12);
    EXPECT_NEAR(sum_dy, 0.0, 1e-12);
    EXPECT_NEAR(sum_dz, 0.0, 1e-12);
  }

  //==========================================================================
  // Higher Order Tests (K = 5, 6)
  //==========================================================================

  TEST(LagrangeBasis1D, LagrangeProperty_K5)
  {
    const auto& nodes = GLL01<5>::getNodes();

    for (size_t i = 0; i <= 5; ++i)
    {
      for (size_t j = 0; j <= 5; ++j)
      {
        Real val = LagrangeBasis1D<5>::getBasis(i, nodes[j], nodes);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(val, expected, 1e-13);
      }
    }
  }

  TEST(LagrangeBasis1D, PartitionOfUnity_K6)
  {
    const auto& nodes = GLL01<6>::getNodes();

    for (Real x = 0.0; x <= 1.0; x += 0.1)
    {
      Real sum = 0.0;
      for (size_t i = 0; i <= 6; ++i)
        sum += LagrangeBasis1D<6>::getBasis(i, x, nodes);
      EXPECT_NEAR(sum, 1.0, 1e-13);
    }
  }

  TEST(LagrangeBasisSegment, LagrangeProperty_K5)
  {
    const auto& nodes = GLL01<5>::getNodes();

    for (size_t i = 0; i <= 5; ++i)
    {
      for (size_t j = 0; j <= 5; ++j)
      {
        Real val = LagrangeBasisSegment<5>::getBasis(i, nodes[j]);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(val, expected, 1e-13);
      }
    }
  }

  TEST(LagrangeBasisSegment, PartitionOfUnity_K6)
  {
    for (Real x = 0.0; x <= 1.0; x += 0.05)
    {
      Real sum = 0.0;
      for (size_t i = 0; i <= 6; ++i)
        sum += LagrangeBasisSegment<6>::getBasis(i, x);
      EXPECT_NEAR(sum, 1.0, 1e-13);
    }
  }

  TEST(LagrangeBasisTriangle, PartitionOfUnity_K5)
  {
    std::vector<std::pair<Real, Real>> test_points = {
      {0.2, 0.2}, {0.1, 0.6}, {0.4, 0.1}, {0.33, 0.33}
    };

    for (auto [x, y] : test_points)
    {
      if (x + y > 1.0) continue;

      Real sum = 0.0;
      for (size_t j = 0; j <= 5; ++j)
      {
        for (size_t i = 0; i <= 5 - j; ++i)
        {
          sum += LagrangeBasisTriangle<5>::getBasis(i, j, x, y);
        }
      }
      EXPECT_NEAR(sum, 1.0, 1e-12) << "x=" << x << " y=" << y;
    }
  }

  TEST(LagrangeBasisTriangle, PartitionOfUnity_K6)
  {
    Real x = 0.25, y = 0.25;
    Real sum = 0.0;

    for (size_t j = 0; j <= 6; ++j)
    {
      for (size_t i = 0; i <= 6 - j; ++i)
      {
        sum += LagrangeBasisTriangle<6>::getBasis(i, j, x, y);
      }
    }
    EXPECT_NEAR(sum, 1.0, 1e-11);
  }

  TEST(LagrangeBasisTetrahedron, PartitionOfUnity_K5)
  {
    Real x = 0.15, y = 0.15, z = 0.15;
    Real sum = 0.0;

    for (size_t k = 0; k <= 5; ++k)
    {
      for (size_t j = 0; j <= 5 - k; ++j)
      {
        for (size_t i = 0; i <= 5 - j - k; ++i)
        {
          sum += LagrangeBasisTetrahedron<5>::getBasis(i, j, k, x, y, z);
        }
      }
    }
    EXPECT_NEAR(sum, 1.0, 1e-11);
  }

  TEST(LagrangeBasisQuadrilateral, LagrangeProperty_K5)
  {
    const auto& nodes = GLL01<5>::getNodes();

    for (size_t i = 0; i <= 5; ++i)
    {
      for (size_t j = 0; j <= 5; ++j)
      {
        for (size_t ii = 0; ii <= 5; ++ii)
        {
          for (size_t jj = 0; jj <= 5; ++jj)
          {
            Real val = LagrangeBasisQuadrilateral<5>::getBasis(i, j, nodes[ii], nodes[jj]);
            Real expected = (i == ii && j == jj) ? 1.0 : 0.0;
            EXPECT_NEAR(val, expected, 1e-12);
          }
        }
      }
    }
  }

  TEST(LagrangeBasisQuadrilateral, PartitionOfUnity_K6)
  {
    for (Real x = 0.1; x <= 0.9; x += 0.2)
    {
      for (Real y = 0.1; y <= 0.9; y += 0.2)
      {
        Real sum = 0.0;
        for (size_t i = 0; i <= 6; ++i)
        {
          for (size_t j = 0; j <= 6; ++j)
          {
            sum += LagrangeBasisQuadrilateral<6>::getBasis(i, j, x, y);
          }
        }
        EXPECT_NEAR(sum, 1.0, 1e-12);
      }
    }
  }

  TEST(LagrangeBasisWedge, PartitionOfUnity_K5)
  {
    Real x = 0.15, y = 0.15, z = 0.5;
    Real sum = 0.0;

    for (size_t k = 0; k <= 5; ++k)
    {
      for (size_t j = 0; j <= 5; ++j)
      {
        for (size_t i = 0; i <= 5 - j; ++i)
        {
          sum += LagrangeBasisWedge<5>::getBasis(i, j, k, x, y, z);
        }
      }
    }
    EXPECT_NEAR(sum, 1.0, 1e-11);
  }

  //==========================================================================
  // Very High Order Tests (K = 15)
  //==========================================================================

  TEST(LagrangeBasis1D, LagrangeProperty_K15)
  {
    const auto& nodes = GLL01<15>::getNodes();

    for (size_t i = 0; i <= 15; ++i)
    {
      for (size_t j = 0; j <= 15; ++j)
      {
        Real val = LagrangeBasis1D<15>::getBasis(i, nodes[j], nodes);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(val, expected, 1e-10);
      }
    }
  }

  TEST(LagrangeBasis1D, PartitionOfUnity_K15)
  {
    const auto& nodes = GLL01<15>::getNodes();

    for (Real x = 0.0; x <= 1.0; x += 0.1)
    {
      Real sum = 0.0;
      for (size_t i = 0; i <= 15; ++i)
        sum += LagrangeBasis1D<15>::getBasis(i, x, nodes);
      EXPECT_NEAR(sum, 1.0, 1e-11);
    }
  }

  TEST(LagrangeBasisSegment, LagrangeProperty_K15)
  {
    const auto& nodes = GLL01<15>::getNodes();

    for (size_t i = 0; i <= 15; ++i)
    {
      for (size_t j = 0; j <= 15; ++j)
      {
        Real val = LagrangeBasisSegment<15>::getBasis(i, nodes[j]);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(val, expected, 1e-10);
      }
    }
  }

  TEST(LagrangeBasisSegment, PartitionOfUnity_K15)
  {
    for (Real x = 0.0; x <= 1.0; x += 0.05)
    {
      Real sum = 0.0;
      for (size_t i = 0; i <= 15; ++i)
        sum += LagrangeBasisSegment<15>::getBasis(i, x);
      EXPECT_NEAR(sum, 1.0, 1e-10);
    }
  }

  TEST(LagrangeBasisTriangle, PartitionOfUnity_K15)
  {
    std::vector<std::pair<Real, Real>> test_points = {
      {0.2, 0.2}, {0.1, 0.6}, {0.4, 0.1}, {0.33, 0.33}
    };

    for (auto [x, y] : test_points)
    {
      if (x + y > 1.0) continue;

      Real sum = 0.0;
      for (size_t j = 0; j <= 15; ++j)
      {
        for (size_t i = 0; i <= 15 - j; ++i)
        {
          sum += LagrangeBasisTriangle<15>::getBasis(i, j, x, y);
        }
      }
      EXPECT_NEAR(sum, 1.0, 1e-8) << "x=" << x << " y=" << y;
    }
  }

  TEST(LagrangeBasisTetrahedron, PartitionOfUnity_K15)
  {
    Real x = 0.15, y = 0.15, z = 0.15;
    Real sum = 0.0;

    for (size_t k = 0; k <= 15; ++k)
    {
      for (size_t j = 0; j <= 15 - k; ++j)
      {
        for (size_t i = 0; i <= 15 - j - k; ++i)
        {
          sum += LagrangeBasisTetrahedron<15>::getBasis(i, j, k, x, y, z);
        }
      }
    }
    EXPECT_NEAR(sum, 1.0, 1e-7);
  }

  TEST(LagrangeBasisQuadrilateral, PartitionOfUnity_K15)
  {
    for (Real x = 0.1; x <= 0.9; x += 0.2)
    {
      for (Real y = 0.1; y <= 0.9; y += 0.2)
      {
        Real sum = 0.0;
        for (size_t i = 0; i <= 15; ++i)
        {
          for (size_t j = 0; j <= 15; ++j)
          {
            sum += LagrangeBasisQuadrilateral<15>::getBasis(i, j, x, y);
          }
        }
        EXPECT_NEAR(sum, 1.0, 1e-9);
      }
    }
  }

  TEST(LagrangeBasisWedge, PartitionOfUnity_K15)
  {
    Real x = 0.15, y = 0.15, z = 0.5;
    Real sum = 0.0;

    for (size_t k = 0; k <= 15; ++k)
    {
      for (size_t j = 0; j <= 15; ++j)
      {
        for (size_t i = 0; i <= 15 - j; ++i)
        {
          sum += LagrangeBasisWedge<15>::getBasis(i, j, k, x, y, z);
        }
      }
    }
    EXPECT_NEAR(sum, 1.0, 1e-7);
  }
}
