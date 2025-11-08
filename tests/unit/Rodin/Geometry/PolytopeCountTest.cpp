/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Geometry/PolytopeCount.h"

using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ============================================================================
  // SimplexCount Tests (PolytopeCount is an alias)
  // ============================================================================

  TEST(Geometry_SimplexCount, DefaultConstruction)
  {
    SimplexCount count;
    // Default constructed, no dimensions initialized
  }

  TEST(Geometry_SimplexCount, ConstructWithMeshDimension)
  {
    SimplexCount count(3);  // 3D mesh: vertices(0), edges(1), faces(2), cells(3)
    
    // Should have 4 dimensions: 0, 1, 2, 3
    EXPECT_EQ(count.at(0), 0);
    EXPECT_EQ(count.at(1), 0);
    EXPECT_EQ(count.at(2), 0);
    EXPECT_EQ(count.at(3), 0);
  }

  TEST(Geometry_SimplexCount, ConstructWithInitializerList)
  {
    SimplexCount count{10, 20, 30};  // 10 vertices, 20 edges, 30 faces
    
    EXPECT_EQ(count.at(0), 10);
    EXPECT_EQ(count.at(1), 20);
    EXPECT_EQ(count.at(2), 30);
  }

  TEST(Geometry_SimplexCount, CopyConstruction)
  {
    SimplexCount count1{5, 10, 15};
    SimplexCount count2(count1);
    
    EXPECT_EQ(count2.at(0), 5);
    EXPECT_EQ(count2.at(1), 10);
    EXPECT_EQ(count2.at(2), 15);
  }

  TEST(Geometry_SimplexCount, MoveConstruction)
  {
    SimplexCount count1{5, 10, 15};
    SimplexCount count2(std::move(count1));
    
    EXPECT_EQ(count2.at(0), 5);
    EXPECT_EQ(count2.at(1), 10);
    EXPECT_EQ(count2.at(2), 15);
  }

  TEST(Geometry_SimplexCount, CopyAssignment)
  {
    SimplexCount count1{5, 10, 15};
    SimplexCount count2;
    
    count2 = count1;
    
    EXPECT_EQ(count2.at(0), 5);
    EXPECT_EQ(count2.at(1), 10);
    EXPECT_EQ(count2.at(2), 15);
  }

  TEST(Geometry_SimplexCount, MoveAssignment)
  {
    SimplexCount count1{5, 10, 15};
    SimplexCount count2;
    
    count2 = std::move(count1);
    
    EXPECT_EQ(count2.at(0), 5);
    EXPECT_EQ(count2.at(1), 10);
    EXPECT_EQ(count2.at(2), 15);
  }

  TEST(Geometry_SimplexCount, Initialize)
  {
    SimplexCount count;
    count.initialize(2);  // 2D mesh: vertices, edges, faces
    
    EXPECT_EQ(count.at(0), 0);
    EXPECT_EQ(count.at(1), 0);
    EXPECT_EQ(count.at(2), 0);
  }

  TEST(Geometry_SimplexCount, AtModification)
  {
    SimplexCount count(2);
    
    count.at(0) = 100;
    count.at(1) = 200;
    count.at(2) = 300;
    
    EXPECT_EQ(count.at(0), 100);
    EXPECT_EQ(count.at(1), 200);
    EXPECT_EQ(count.at(2), 300);
  }

  TEST(Geometry_SimplexCount, AtConst)
  {
    const SimplexCount count{42, 84, 126};
    
    EXPECT_EQ(count.at(0), 42);
    EXPECT_EQ(count.at(1), 84);
    EXPECT_EQ(count.at(2), 126);
  }

  TEST(Geometry_SimplexCount, Reinitialize)
  {
    SimplexCount count{1, 2, 3};
    
    count.initialize(1);  // Re-initialize to 1D mesh (2 dimensions: 0 and 1)
    
    // initialize() resizes the vector but preserves existing values
    // So the first two values remain: 1, 2
    EXPECT_EQ(count.at(0), 1);
    EXPECT_EQ(count.at(1), 2);
  }

  TEST(Geometry_SimplexCount, ChainedInitialize)
  {
    SimplexCount count;
    
    // Test that initialize returns reference for chaining
    count.initialize(3).at(0) = 10;
    
    EXPECT_EQ(count.at(0), 10);
  }
}
