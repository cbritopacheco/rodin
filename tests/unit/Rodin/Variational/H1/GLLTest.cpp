/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational/H1/GLL.h"

using namespace Rodin;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  //==========================================================================
  // GLL Node Count Tests
  //==========================================================================

  TEST(GLL, NodeCount_K0)
  {
    EXPECT_EQ(GLL<0>::getCount(), 1);
  }

  TEST(GLL, NodeCount_K1)
  {
    EXPECT_EQ(GLL<1>::getCount(), 2);
  }

  TEST(GLL, NodeCount_K2)
  {
    EXPECT_EQ(GLL<2>::getCount(), 3);
  }

  TEST(GLL, NodeCount_K5)
  {
    EXPECT_EQ(GLL<5>::getCount(), 6);
  }

  TEST(GLL, NodeCount_K10)
  {
    EXPECT_EQ(GLL<10>::getCount(), 11);
  }

  //==========================================================================
  // GLL Node Endpoint Tests (K >= 1)
  //==========================================================================

  TEST(GLL, Endpoints_K1)
  {
    const auto& nodes = GLL<1>::getNodes();
    EXPECT_NEAR(nodes[0], -1.0, 1e-14);
    EXPECT_NEAR(nodes[1], 1.0, 1e-14);
  }

  TEST(GLL, Endpoints_K2)
  {
    const auto& nodes = GLL<2>::getNodes();
    EXPECT_NEAR(nodes[0], -1.0, 1e-14);
    EXPECT_NEAR(nodes[2], 1.0, 1e-14);
  }

  TEST(GLL, Endpoints_K5)
  {
    const auto& nodes = GLL<5>::getNodes();
    EXPECT_NEAR(nodes[0], -1.0, 1e-14);
    EXPECT_NEAR(nodes[5], 1.0, 1e-14);
  }

  TEST(GLL, Endpoints_K10)
  {
    const auto& nodes = GLL<10>::getNodes();
    EXPECT_NEAR(nodes[0], -1.0, 1e-14);
    EXPECT_NEAR(nodes[10], 1.0, 1e-14);
  }

  //==========================================================================
  // GLL Node Ordering Tests (ascending)
  //==========================================================================

  TEST(GLL, Ascending_K2)
  {
    const auto& nodes = GLL<2>::getNodes();
    EXPECT_LT(nodes[0], nodes[1]);
    EXPECT_LT(nodes[1], nodes[2]);
  }

  TEST(GLL, Ascending_K5)
  {
    const auto& nodes = GLL<5>::getNodes();
    for (size_t i = 0; i < 5; ++i)
    {
      EXPECT_LT(nodes[i], nodes[i + 1]);
    }
  }

  TEST(GLL, Ascending_K10)
  {
    const auto& nodes = GLL<10>::getNodes();
    for (size_t i = 0; i < 10; ++i)
    {
      EXPECT_LT(nodes[i], nodes[i + 1]);
    }
  }

  //==========================================================================
  // GLL Node Symmetry Tests
  //==========================================================================

  TEST(GLL, Symmetry_K2)
  {
    // GLL nodes are symmetric about 0
    const auto& nodes = GLL<2>::getNodes();
    EXPECT_NEAR(nodes[1], 0.0, 1e-14);  // Middle node at 0
    EXPECT_NEAR(nodes[0], -nodes[2], 1e-14);
  }

  TEST(GLL, Symmetry_K3)
  {
    const auto& nodes = GLL<3>::getNodes();
    EXPECT_NEAR(nodes[0], -nodes[3], 1e-14);
    EXPECT_NEAR(nodes[1], -nodes[2], 1e-14);
  }

  TEST(GLL, Symmetry_K4)
  {
    const auto& nodes = GLL<4>::getNodes();
    EXPECT_NEAR(nodes[2], 0.0, 1e-14);  // Middle node at 0
    EXPECT_NEAR(nodes[0], -nodes[4], 1e-14);
    EXPECT_NEAR(nodes[1], -nodes[3], 1e-14);
  }

  TEST(GLL, Symmetry_K5)
  {
    const auto& nodes = GLL<5>::getNodes();
    EXPECT_NEAR(nodes[0], -nodes[5], 1e-14);
    EXPECT_NEAR(nodes[1], -nodes[4], 1e-14);
    EXPECT_NEAR(nodes[2], -nodes[3], 1e-14);
  }

  //==========================================================================
  // GLL Known Values Tests
  //==========================================================================

  TEST(GLL, KnownValues_K2)
  {
    // K=2: nodes are -1, 0, 1
    const auto& nodes = GLL<2>::getNodes();
    EXPECT_NEAR(nodes[0], -1.0, 1e-14);
    EXPECT_NEAR(nodes[1], 0.0, 1e-14);
    EXPECT_NEAR(nodes[2], 1.0, 1e-14);
  }

  TEST(GLL, KnownValues_K3)
  {
    // K=3: nodes are -1, -1/sqrt(5), 1/sqrt(5), 1
    const auto& nodes = GLL<3>::getNodes();
    const Real sqrt5 = std::sqrt(5.0);
    EXPECT_NEAR(nodes[0], -1.0, 1e-14);
    EXPECT_NEAR(nodes[1], -1.0 / sqrt5, 1e-14);
    EXPECT_NEAR(nodes[2], 1.0 / sqrt5, 1e-14);
    EXPECT_NEAR(nodes[3], 1.0, 1e-14);
  }

  TEST(GLL, KnownValues_K4)
  {
    // K=4: nodes are -1, -sqrt(3/7), 0, sqrt(3/7), 1
    const auto& nodes = GLL<4>::getNodes();
    const Real val = std::sqrt(3.0 / 7.0);
    EXPECT_NEAR(nodes[0], -1.0, 1e-14);
    EXPECT_NEAR(nodes[1], -val, 1e-14);
    EXPECT_NEAR(nodes[2], 0.0, 1e-14);
    EXPECT_NEAR(nodes[3], val, 1e-14);
    EXPECT_NEAR(nodes[4], 1.0, 1e-14);
  }

  //==========================================================================
  // GLL01 (mapped to [0,1]) Tests
  //==========================================================================

  TEST(GLL01, NodeCount)
  {
    EXPECT_EQ(GLL01<0>::getCount(), 1);
    EXPECT_EQ(GLL01<2>::getCount(), 3);
    EXPECT_EQ(GLL01<5>::getCount(), 6);
  }

  TEST(GLL01, Endpoints_K1)
  {
    const auto& nodes = GLL01<1>::getNodes();
    EXPECT_NEAR(nodes[0], 0.0, 1e-14);
    EXPECT_NEAR(nodes[1], 1.0, 1e-14);
  }

  TEST(GLL01, Endpoints_K2)
  {
    const auto& nodes = GLL01<2>::getNodes();
    EXPECT_NEAR(nodes[0], 0.0, 1e-14);
    EXPECT_NEAR(nodes[2], 1.0, 1e-14);
  }

  TEST(GLL01, Endpoints_K5)
  {
    const auto& nodes = GLL01<5>::getNodes();
    EXPECT_NEAR(nodes[0], 0.0, 1e-14);
    EXPECT_NEAR(nodes[5], 1.0, 1e-14);
  }

  TEST(GLL01, KnownValues_K2)
  {
    // K=2 on [0,1]: 0, 0.5, 1
    const auto& nodes = GLL01<2>::getNodes();
    EXPECT_NEAR(nodes[0], 0.0, 1e-14);
    EXPECT_NEAR(nodes[1], 0.5, 1e-14);
    EXPECT_NEAR(nodes[2], 1.0, 1e-14);
  }

  TEST(GLL01, MappingConsistency)
  {
    // GLL01 = (GLL + 1) / 2
    const auto& gll = GLL<5>::getNodes();
    const auto& gll01 = GLL01<5>::getNodes();

    for (size_t i = 0; i <= 5; ++i)
    {
      Real expected = (gll[i] + 1.0) / 2.0;
      EXPECT_NEAR(gll01[i], expected, 1e-14);
    }
  }

  TEST(GLL01, Ascending_K5)
  {
    const auto& nodes = GLL01<5>::getNodes();
    for (size_t i = 0; i < 5; ++i)
    {
      EXPECT_LT(nodes[i], nodes[i + 1]);
    }
  }

  TEST(GLL01, Range_K10)
  {
    // All nodes should be in [0, 1]
    const auto& nodes = GLL01<10>::getNodes();
    for (size_t i = 0; i <= 10; ++i)
    {
      EXPECT_GE(nodes[i], 0.0);
      EXPECT_LE(nodes[i], 1.0);
    }
  }

  //==========================================================================
  // GLL Interior Node Property Tests
  //==========================================================================

  TEST(GLL, InteriorNodesInInterval_K5)
  {
    // All interior nodes should be in (-1, 1)
    const auto& nodes = GLL<5>::getNodes();
    for (size_t i = 1; i < 5; ++i)
    {
      EXPECT_GT(nodes[i], -1.0);
      EXPECT_LT(nodes[i], 1.0);
    }
  }

  TEST(GLL, InteriorNodesInInterval_K10)
  {
    const auto& nodes = GLL<10>::getNodes();
    for (size_t i = 1; i < 10; ++i)
    {
      EXPECT_GT(nodes[i], -1.0);
      EXPECT_LT(nodes[i], 1.0);
    }
  }

  //==========================================================================
  // GLL Static Access Tests
  //==========================================================================

  TEST(GLL, GetNode_K3)
  {
    EXPECT_NEAR(GLL<3>::getNode(0), -1.0, 1e-14);
    EXPECT_NEAR(GLL<3>::getNode(3), 1.0, 1e-14);

    const Real sqrt5 = std::sqrt(5.0);
    EXPECT_NEAR(GLL<3>::getNode(1), -1.0 / sqrt5, 1e-14);
    EXPECT_NEAR(GLL<3>::getNode(2), 1.0 / sqrt5, 1e-14);
  }

  TEST(GLL01, GetNode_K3)
  {
    EXPECT_NEAR(GLL01<3>::getNode(0), 0.0, 1e-14);
    EXPECT_NEAR(GLL01<3>::getNode(3), 1.0, 1e-14);

    const Real sqrt5 = std::sqrt(5.0);
    EXPECT_NEAR(GLL01<3>::getNode(1), (1.0 - 1.0/sqrt5) / 2.0, 1e-14);
    EXPECT_NEAR(GLL01<3>::getNode(2), (1.0 + 1.0/sqrt5) / 2.0, 1e-14);
  }

  //==========================================================================
  // Higher Order Tests (K = 5, 6)
  //==========================================================================

  TEST(GLL, KnownValues_K5)
  {
    // K=5: 6 nodes including -1 and 1
    const auto& nodes = GLL<5>::getNodes();
    EXPECT_NEAR(nodes[0], -1.0, 1e-14);
    EXPECT_NEAR(nodes[5], 1.0, 1e-14);

    // Interior nodes are roots of P'_5(x), should be symmetric
    EXPECT_NEAR(nodes[1], -nodes[4], 1e-14);
    EXPECT_NEAR(nodes[2], -nodes[3], 1e-14);
  }

  TEST(GLL, Symmetry_K6)
  {
    const auto& nodes = GLL<6>::getNodes();

    // Middle node should be at 0
    EXPECT_NEAR(nodes[3], 0.0, 1e-14);

    // Check symmetry
    EXPECT_NEAR(nodes[0], -nodes[6], 1e-14);
    EXPECT_NEAR(nodes[1], -nodes[5], 1e-14);
    EXPECT_NEAR(nodes[2], -nodes[4], 1e-14);
  }

  TEST(GLL, Ascending_K6)
  {
    const auto& nodes = GLL<6>::getNodes();
    for (size_t i = 0; i < 6; ++i)
    {
      EXPECT_LT(nodes[i], nodes[i + 1]);
    }
  }

  TEST(GLL, Endpoints_K6)
  {
    const auto& nodes = GLL<6>::getNodes();
    EXPECT_NEAR(nodes[0], -1.0, 1e-14);
    EXPECT_NEAR(nodes[6], 1.0, 1e-14);
  }

  TEST(GLL, InteriorNodesInInterval_K6)
  {
    const auto& nodes = GLL<6>::getNodes();
    for (size_t i = 1; i < 6; ++i)
    {
      EXPECT_GT(nodes[i], -1.0);
      EXPECT_LT(nodes[i], 1.0);
    }
  }

  TEST(GLL01, Endpoints_K6)
  {
    const auto& nodes = GLL01<6>::getNodes();
    EXPECT_NEAR(nodes[0], 0.0, 1e-14);
    EXPECT_NEAR(nodes[6], 1.0, 1e-14);
  }

  TEST(GLL01, MappingConsistency_K6)
  {
    // GLL01 = (GLL + 1) / 2
    const auto& gll = GLL<6>::getNodes();
    const auto& gll01 = GLL01<6>::getNodes();

    for (size_t i = 0; i <= 6; ++i)
    {
      Real expected = (gll[i] + 1.0) / 2.0;
      EXPECT_NEAR(gll01[i], expected, 1e-14);
    }
  }

  TEST(GLL01, Range_K6)
  {
    // All nodes should be in [0, 1]
    const auto& nodes = GLL01<6>::getNodes();
    for (size_t i = 0; i <= 6; ++i)
    {
      EXPECT_GE(nodes[i], 0.0);
      EXPECT_LE(nodes[i], 1.0);
    }
  }

  TEST(GLL, NodeCount_K6)
  {
    EXPECT_EQ(GLL<6>::getCount(), 7);
  }

  TEST(GLL01, NodeCount_K6)
  {
    EXPECT_EQ(GLL01<6>::getCount(), 7);
  }
}
