/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Geometry/Region.h"

using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ============================================================================
  // Region Enum Tests
  // ============================================================================

  TEST(Geometry_Region, EnumValues)
  {
    // Test that enum values are distinct
    EXPECT_NE(static_cast<int>(Region::Cells), static_cast<int>(Region::Faces));
    EXPECT_NE(static_cast<int>(Region::Cells), static_cast<int>(Region::Boundary));
    EXPECT_NE(static_cast<int>(Region::Cells), static_cast<int>(Region::Interface));
    EXPECT_NE(static_cast<int>(Region::Faces), static_cast<int>(Region::Boundary));
    EXPECT_NE(static_cast<int>(Region::Faces), static_cast<int>(Region::Interface));
    EXPECT_NE(static_cast<int>(Region::Boundary), static_cast<int>(Region::Interface));
  }

  TEST(Geometry_Region, AssignmentAndComparison)
  {
    Region r1 = Region::Cells;
    Region r2 = Region::Cells;
    Region r3 = Region::Faces;
    
    EXPECT_EQ(r1, r2);
    EXPECT_NE(r1, r3);
  }

  TEST(Geometry_Region, AllValuesAccessible)
  {
    // Test that all enum values can be assigned
    Region cells = Region::Cells;
    Region faces = Region::Faces;
    Region boundary = Region::Boundary;
    Region interface = Region::Interface;
    
    EXPECT_EQ(cells, Region::Cells);
    EXPECT_EQ(faces, Region::Faces);
    EXPECT_EQ(boundary, Region::Boundary);
    EXPECT_EQ(interface, Region::Interface);
  }
}
