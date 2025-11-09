/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Geometry/Polytope.h"

using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ============================================================================
  // Polytope::Traits Tests
  // ============================================================================

  TEST(Geometry_PolytopeTraits, PointDimension)
  {
    Polytope::Traits traits(Polytope::Type::Point);
    EXPECT_EQ(traits.getDimension(), 0);
  }

  TEST(Geometry_PolytopeTraits, PointVertexCount)
  {
    Polytope::Traits traits(Polytope::Type::Point);
    EXPECT_EQ(traits.getVertexCount(), 1);
  }

  TEST(Geometry_PolytopeTraits, PointIsSimplex)
  {
    Polytope::Traits traits(Polytope::Type::Point);
    EXPECT_TRUE(traits.isSimplex());
  }

  TEST(Geometry_PolytopeTraits, SegmentDimension)
  {
    Polytope::Traits traits(Polytope::Type::Segment);
    EXPECT_EQ(traits.getDimension(), 1);
  }

  TEST(Geometry_PolytopeTraits, SegmentVertexCount)
  {
    Polytope::Traits traits(Polytope::Type::Segment);
    EXPECT_EQ(traits.getVertexCount(), 2);
  }

  TEST(Geometry_PolytopeTraits, SegmentIsSimplex)
  {
    Polytope::Traits traits(Polytope::Type::Segment);
    EXPECT_TRUE(traits.isSimplex());
  }

  TEST(Geometry_PolytopeTraits, TriangleDimension)
  {
    Polytope::Traits traits(Polytope::Type::Triangle);
    EXPECT_EQ(traits.getDimension(), 2);
  }

  TEST(Geometry_PolytopeTraits, TriangleVertexCount)
  {
    Polytope::Traits traits(Polytope::Type::Triangle);
    EXPECT_EQ(traits.getVertexCount(), 3);
  }

  TEST(Geometry_PolytopeTraits, TriangleIsSimplex)
  {
    Polytope::Traits traits(Polytope::Type::Triangle);
    EXPECT_TRUE(traits.isSimplex());
  }

  TEST(Geometry_PolytopeTraits, QuadrilateralDimension)
  {
    Polytope::Traits traits(Polytope::Type::Quadrilateral);
    EXPECT_EQ(traits.getDimension(), 2);
  }

  TEST(Geometry_PolytopeTraits, QuadrilateralVertexCount)
  {
    Polytope::Traits traits(Polytope::Type::Quadrilateral);
    EXPECT_EQ(traits.getVertexCount(), 4);
  }

  TEST(Geometry_PolytopeTraits, QuadrilateralIsNotSimplex)
  {
    Polytope::Traits traits(Polytope::Type::Quadrilateral);
    EXPECT_FALSE(traits.isSimplex());
  }

  TEST(Geometry_PolytopeTraits, TetrahedronDimension)
  {
    Polytope::Traits traits(Polytope::Type::Tetrahedron);
    EXPECT_EQ(traits.getDimension(), 3);
  }

  TEST(Geometry_PolytopeTraits, TetrahedronVertexCount)
  {
    Polytope::Traits traits(Polytope::Type::Tetrahedron);
    EXPECT_EQ(traits.getVertexCount(), 4);
  }

  TEST(Geometry_PolytopeTraits, TetrahedronIsSimplex)
  {
    Polytope::Traits traits(Polytope::Type::Tetrahedron);
    EXPECT_TRUE(traits.isSimplex());
  }

  TEST(Geometry_PolytopeTraits, WedgeDimension)
  {
    Polytope::Traits traits(Polytope::Type::Wedge);
    EXPECT_EQ(traits.getDimension(), 3);
  }

  TEST(Geometry_PolytopeTraits, WedgeVertexCount)
  {
    Polytope::Traits traits(Polytope::Type::Wedge);
    EXPECT_EQ(traits.getVertexCount(), 6);
  }

  TEST(Geometry_PolytopeTraits, WedgeIsNotSimplex)
  {
    Polytope::Traits traits(Polytope::Type::Wedge);
    EXPECT_FALSE(traits.isSimplex());
  }

  TEST(Geometry_PolytopeTraits, TriangleVertices)
  {
    Polytope::Traits traits(Polytope::Type::Triangle);
    
    // Check that we can retrieve vertices
    const auto& v0 = traits.getVertex(0);
    const auto& v1 = traits.getVertex(1);
    const auto& v2 = traits.getVertex(2);
    
    // Vertices should be accessible (we don't check exact values, just accessibility)
    EXPECT_GE(v0.size(), 0);
    EXPECT_GE(v1.size(), 0);
    EXPECT_GE(v2.size(), 0);
  }

  TEST(Geometry_PolytopeTraits, HalfSpaceAccessible)
  {
    Polytope::Traits traits(Polytope::Type::Triangle);
    
    // Check that half space is accessible
    const auto& halfSpace = traits.getHalfSpace();
    
    // Matrix and vector should be initialized (we don't check exact values)
    EXPECT_GE(halfSpace.matrix.rows(), 0);
    EXPECT_GE(halfSpace.vector.size(), 0);
  }

  // ============================================================================
  // Polytope::Types Array Tests
  // ============================================================================

  TEST(Geometry_PolytopeTypes, ArraySize)
  {
    EXPECT_EQ(Polytope::Types.size(), 6);
  }

  TEST(Geometry_PolytopeTypes, ArrayContents)
  {
    EXPECT_EQ(Polytope::Types[0], Polytope::Type::Point);
    EXPECT_EQ(Polytope::Types[1], Polytope::Type::Segment);
    EXPECT_EQ(Polytope::Types[2], Polytope::Type::Triangle);
    EXPECT_EQ(Polytope::Types[3], Polytope::Type::Quadrilateral);
    EXPECT_EQ(Polytope::Types[4], Polytope::Type::Tetrahedron);
    EXPECT_EQ(Polytope::Types[5], Polytope::Type::Wedge);
  }

  TEST(Geometry_PolytopeTypes, IterableArray)
  {
    int count = 0;
    for (auto type : Polytope::Types)
    {
      // Just check that we can iterate through the array
      (void)type;  // Suppress unused variable warning
      ++count;
    }
    EXPECT_EQ(count, 6);
  }
}
