/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry/Polytope.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ---- Segment projections ----

  TEST(Geometry_PolytopeProject, Segment_Cell)
  {
    Polytope::Project proj(Polytope::Type::Segment);
    Math::SpatialPoint rc{0.5};
    Math::SpatialPoint out(1);
    proj.cell(out, rc);
    EXPECT_NEAR(out(0), 0.5, 1e-14);
  }

  TEST(Geometry_PolytopeProject, Segment_Vertex0)
  {
    Polytope::Project proj(Polytope::Type::Segment);
    // Vertex 0 of segment is at reference coordinate 0
    Math::SpatialPoint rc(0);
    Math::SpatialPoint out(1);
    proj.vertex(0, out, rc);
    EXPECT_NEAR(out(0), 0.0, 1e-14);
  }

  TEST(Geometry_PolytopeProject, Segment_Vertex1)
  {
    Polytope::Project proj(Polytope::Type::Segment);
    Math::SpatialPoint rc(0);
    Math::SpatialPoint out(1);
    proj.vertex(1, out, rc);
    EXPECT_NEAR(out(0), 1.0, 1e-14);
  }

  // ---- Triangle projections ----

  TEST(Geometry_PolytopeProject, Triangle_Cell)
  {
    Polytope::Project proj(Polytope::Type::Triangle);
    Math::SpatialPoint rc{0.25, 0.25};
    Math::SpatialPoint out(2);
    proj.cell(out, rc);
    EXPECT_NEAR(out(0), 0.25, 1e-14);
    EXPECT_NEAR(out(1), 0.25, 1e-14);
  }

  TEST(Geometry_PolytopeProject, Triangle_Face0)
  {
    Polytope::Project proj(Polytope::Type::Triangle);
    // rc must match the geometry dimension (2D for triangle)
    Math::SpatialPoint rc{0.5, 0.3};
    Math::SpatialPoint out(2);
    proj.face(0, out, rc);
    // Should map to a point on the first face of the reference triangle
    EXPECT_GE(out(0), -1e-14);
    EXPECT_GE(out(1), -1e-14);
    EXPECT_LE(out(0) + out(1), 1.0 + 1e-14);
  }

  TEST(Geometry_PolytopeProject, Triangle_Face1)
  {
    Polytope::Project proj(Polytope::Type::Triangle);
    Math::SpatialPoint rc{0.5, 0.3};
    Math::SpatialPoint out(2);
    proj.face(1, out, rc);
    EXPECT_GE(out(0), -1e-14);
    EXPECT_GE(out(1), -1e-14);
    EXPECT_LE(out(0) + out(1), 1.0 + 1e-14);
  }

  TEST(Geometry_PolytopeProject, Triangle_Face2)
  {
    Polytope::Project proj(Polytope::Type::Triangle);
    Math::SpatialPoint rc{0.5, 0.3};
    Math::SpatialPoint out(2);
    proj.face(2, out, rc);
    EXPECT_GE(out(0), -1e-14);
    EXPECT_GE(out(1), -1e-14);
    EXPECT_LE(out(0) + out(1), 1.0 + 1e-14);
  }

  TEST(Geometry_PolytopeProject, Triangle_Vertex0)
  {
    Polytope::Project proj(Polytope::Type::Triangle);
    Math::SpatialPoint rc(0);
    Math::SpatialPoint out(2);
    proj.vertex(0, out, rc);
    // Vertex 0 of reference triangle is at (0, 0)
    EXPECT_NEAR(out(0), 0.0, 1e-14);
    EXPECT_NEAR(out(1), 0.0, 1e-14);
  }

  TEST(Geometry_PolytopeProject, Triangle_Vertex1)
  {
    Polytope::Project proj(Polytope::Type::Triangle);
    Math::SpatialPoint rc(0);
    Math::SpatialPoint out(2);
    proj.vertex(1, out, rc);
    // Vertex 1 of reference triangle is at (1, 0)
    EXPECT_NEAR(out(0), 1.0, 1e-14);
    EXPECT_NEAR(out(1), 0.0, 1e-14);
  }

  TEST(Geometry_PolytopeProject, Triangle_Vertex2)
  {
    Polytope::Project proj(Polytope::Type::Triangle);
    Math::SpatialPoint rc(0);
    Math::SpatialPoint out(2);
    proj.vertex(2, out, rc);
    // Vertex 2 of reference triangle is at (0, 1)
    EXPECT_NEAR(out(0), 0.0, 1e-14);
    EXPECT_NEAR(out(1), 1.0, 1e-14);
  }

  // ---- Quadrilateral projections ----

  TEST(Geometry_PolytopeProject, Quadrilateral_Cell)
  {
    Polytope::Project proj(Polytope::Type::Quadrilateral);
    Math::SpatialPoint rc{0.5, 0.5};
    Math::SpatialPoint out(2);
    proj.cell(out, rc);
    EXPECT_NEAR(out(0), 0.5, 1e-14);
    EXPECT_NEAR(out(1), 0.5, 1e-14);
  }

  TEST(Geometry_PolytopeProject, Quadrilateral_Vertex0)
  {
    Polytope::Project proj(Polytope::Type::Quadrilateral);
    Math::SpatialPoint rc(0);
    Math::SpatialPoint out(2);
    proj.vertex(0, out, rc);
    // Vertex 0 of reference quad is at (0, 0)
    EXPECT_NEAR(out(0), 0.0, 1e-14);
    EXPECT_NEAR(out(1), 0.0, 1e-14);
  }

  // ---- Tetrahedron projections ----

  TEST(Geometry_PolytopeProject, Tetrahedron_Cell)
  {
    Polytope::Project proj(Polytope::Type::Tetrahedron);
    Math::SpatialPoint rc{0.1, 0.1, 0.1};
    Math::SpatialPoint out(3);
    proj.cell(out, rc);
    EXPECT_NEAR(out(0), 0.1, 1e-14);
    EXPECT_NEAR(out(1), 0.1, 1e-14);
    EXPECT_NEAR(out(2), 0.1, 1e-14);
  }

  TEST(Geometry_PolytopeProject, Tetrahedron_Face0)
  {
    Polytope::Project proj(Polytope::Type::Tetrahedron);
    // rc must be 3D for tetrahedron
    Math::SpatialPoint rc{0.25, 0.25, 0.1};
    Math::SpatialPoint out(3);
    proj.face(0, out, rc);
    // Should be a valid point inside the reference tet
    EXPECT_GE(out(0), -1e-14);
    EXPECT_GE(out(1), -1e-14);
    EXPECT_GE(out(2), -1e-14);
    EXPECT_LE(out(0) + out(1) + out(2), 1.0 + 1e-14);
  }

  TEST(Geometry_PolytopeProject, Tetrahedron_Vertex0)
  {
    Polytope::Project proj(Polytope::Type::Tetrahedron);
    Math::SpatialPoint rc(0);
    Math::SpatialPoint out(3);
    proj.vertex(0, out, rc);
    // Vertex 0 of reference tet is at (0, 0, 0)
    EXPECT_NEAR(out(0), 0.0, 1e-14);
    EXPECT_NEAR(out(1), 0.0, 1e-14);
    EXPECT_NEAR(out(2), 0.0, 1e-14);
  }

  // ---- Hexahedron projections ----

  TEST(Geometry_PolytopeProject, Hexahedron_Cell)
  {
    Polytope::Project proj(Polytope::Type::Hexahedron);
    Math::SpatialPoint rc{0.5, 0.5, 0.5};
    Math::SpatialPoint out(3);
    proj.cell(out, rc);
    EXPECT_NEAR(out(0), 0.5, 1e-14);
    EXPECT_NEAR(out(1), 0.5, 1e-14);
    EXPECT_NEAR(out(2), 0.5, 1e-14);
  }

  TEST(Geometry_PolytopeProject, Hexahedron_Vertex0)
  {
    Polytope::Project proj(Polytope::Type::Hexahedron);
    Math::SpatialPoint rc(0);
    Math::SpatialPoint out(3);
    proj.vertex(0, out, rc);
    // Vertex 0 of reference hex is at (0, 0, 0)
    EXPECT_NEAR(out(0), 0.0, 1e-14);
    EXPECT_NEAR(out(1), 0.0, 1e-14);
    EXPECT_NEAR(out(2), 0.0, 1e-14);
  }

  // ---- Boundary projection ----

  TEST(Geometry_PolytopeProject, Triangle_Boundary)
  {
    Polytope::Project proj(Polytope::Type::Triangle);
    // rc must be 2D for triangle
    Math::SpatialPoint rc{0.5, 0.3};
    Math::SpatialPoint out(2);
    proj.boundary(out, rc);
    // Boundary projection should produce a point in the reference triangle
    EXPECT_GE(out(0), -1e-14);
    EXPECT_GE(out(1), -1e-14);
    EXPECT_LE(out(0) + out(1), 1.0 + 1e-14);
  }
}
