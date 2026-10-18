/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <array>
#include <limits>

#include <Rodin/Geometry/Polytope.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Projection onto a triangle edge preserves every point already on that edge.
  TEST(Geometry_PolytopeProject, TriangleFaceProjectionPreservesFacePoints)
  {
    const Polytope::Project project(Polytope::Type::Triangle);
    for (Real t : {0.0, 0.25, 0.75, 1.0})
    {
      SCOPED_TRACE(t);
      const std::array<Math::SpatialPoint, 3> points{
        Math::SpatialPoint{t, 0}, Math::SpatialPoint{1 - t, t}, Math::SpatialPoint{0, t}};
      for (size_t face = 0; face < points.size(); ++face)
      {
        SCOPED_TRACE(face);
        Math::SpatialPoint result;
        project.face(face, result, points[face]);
        ASSERT_EQ(result.size(), 2);
        EXPECT_DOUBLE_EQ(result[0], points[face][0]);
        EXPECT_DOUBLE_EQ(result[1], points[face][1]);
      }
    }
  }

  /// @brief Projection onto a square edge preserves its parent reference coordinates.
  TEST(Geometry_PolytopeProject, QuadrilateralFaceProjectionPreservesFacePoints)
  {
    const Polytope::Project project(Polytope::Type::Quadrilateral);
    for (Real t : {0.0, 0.25, 0.75, 1.0})
    {
      SCOPED_TRACE(t);
      const std::array<Math::SpatialPoint, 4> points{Math::SpatialPoint{t, 0},
        Math::SpatialPoint{1, t}, Math::SpatialPoint{t, 1}, Math::SpatialPoint{0, t}};
      for (size_t face = 0; face < points.size(); ++face)
      {
        SCOPED_TRACE(face);
        Math::SpatialPoint result;
        project.face(face, result, points[face]);
        ASSERT_EQ(result.size(), 2);
        EXPECT_DOUBLE_EQ(result[0], points[face][0]);
        EXPECT_DOUBLE_EQ(result[1], points[face][1]);
      }
    }
  }

  class GeometryFaceProjectionValidity : public ::testing::TestWithParam<Polytope::Type>
  {};

  /**
   * @brief Exterior inputs project into the requested reference face.
   *
   * Sixteen rounding units allow the short coordinate projections and
   * halfspace products in dimensions one through three.
   */
  TEST_P(GeometryFaceProjectionValidity, ExteriorPointsStayOnRequestedFace)
  {
    constexpr Real tolerance = 16 * std::numeric_limits<Real>::epsilon();
    constexpr std::array<Real, 6> coordinates{-3, -1, 0.25, 0.75, 2, 4};
    const auto geometry = GetParam();
    const Polytope::Traits traits(geometry);
    const size_t dimension = traits.getDimension();
    const auto& halfSpace = traits.getHalfSpace();
    const Polytope::Project project(geometry);
    size_t count = 1;
    for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
      count *= coordinates.size();
    for (size_t face = 0; face < static_cast<size_t>(halfSpace.vector.size()); ++face)
    {
      SCOPED_TRACE(face);
      for (size_t sample = 0; sample < count; ++sample)
      {
        SCOPED_TRACE(sample);
        Math::SpatialPoint input(dimension);
        size_t remaining = sample;
        for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
        {
          input[coordinate] = coordinates[remaining % coordinates.size()];
          remaining /= coordinates.size();
        }
        Math::SpatialPoint result;
        project.face(face, result, input);
        ASSERT_EQ(static_cast<size_t>(result.size()), dimension);
        const Math::Vector<Real> residual =
          halfSpace.matrix * result.getData().head(dimension) - halfSpace.vector;
        ASSERT_LE(residual.maxCoeff(), tolerance);
        ASSERT_NEAR(residual[face], 0, tolerance);
      }
    }
  }

  INSTANTIATE_TEST_SUITE_P(PositiveDimensions, GeometryFaceProjectionValidity,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge));

  // ---- Segment projections ----

  /// @brief Verifies segment cell for geometry polytope project by checking tolerance-based numerical results.
  TEST(Geometry_PolytopeProject, Segment_Cell)
  {
    Polytope::Project proj(Polytope::Type::Segment);
    Math::SpatialPoint rc{0.5};
    Math::SpatialPoint out(1);
    proj.cell(out, rc);
    EXPECT_NEAR(out(0), 0.5, 1e-14);
  }

  /// @brief Verifies segment vertex 0 for geometry polytope project by checking tolerance-based numerical results.
  TEST(Geometry_PolytopeProject, Segment_Vertex0)
  {
    Polytope::Project proj(Polytope::Type::Segment);
    // Vertex 0 of segment is at reference coordinate 0
    Math::SpatialPoint rc(0);
    Math::SpatialPoint out(1);
    proj.vertex(0, out, rc);
    EXPECT_NEAR(out(0), 0.0, 1e-14);
  }

  /// @brief Verifies segment vertex 1 for geometry polytope project by checking tolerance-based numerical results.
  TEST(Geometry_PolytopeProject, Segment_Vertex1)
  {
    Polytope::Project proj(Polytope::Type::Segment);
    Math::SpatialPoint rc(0);
    Math::SpatialPoint out(1);
    proj.vertex(1, out, rc);
    EXPECT_NEAR(out(0), 1.0, 1e-14);
  }

  // ---- Triangle projections ----

  /// @brief Verifies triangle cell for geometry polytope project by checking tolerance-based numerical results.
  TEST(Geometry_PolytopeProject, Triangle_Cell)
  {
    Polytope::Project proj(Polytope::Type::Triangle);
    Math::SpatialPoint rc{0.25, 0.25};
    Math::SpatialPoint out(2);
    proj.cell(out, rc);
    EXPECT_NEAR(out(0), 0.25, 1e-14);
    EXPECT_NEAR(out(1), 0.25, 1e-14);
  }

  /// @brief Verifies triangle face 0 for geometry polytope project.
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

  /// @brief Verifies triangle face 1 for geometry polytope project.
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

  /// @brief Verifies triangle face 2 for geometry polytope project.
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

  /// @brief Verifies triangle vertex 0 for geometry polytope project by checking tolerance-based numerical results.
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

  /// @brief Verifies triangle vertex 1 for geometry polytope project by checking tolerance-based numerical results.
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

  /// @brief Verifies triangle vertex 2 for geometry polytope project by checking tolerance-based numerical results.
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

  /// @brief Verifies quadrilateral cell for geometry polytope project by checking tolerance-based numerical results.
  TEST(Geometry_PolytopeProject, Quadrilateral_Cell)
  {
    Polytope::Project proj(Polytope::Type::Quadrilateral);
    Math::SpatialPoint rc{0.5, 0.5};
    Math::SpatialPoint out(2);
    proj.cell(out, rc);
    EXPECT_NEAR(out(0), 0.5, 1e-14);
    EXPECT_NEAR(out(1), 0.5, 1e-14);
  }

  /// @brief Verifies quadrilateral vertex 0 for geometry polytope project by checking tolerance-based numerical results.
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

  /// @brief Verifies tetrahedron cell for geometry polytope project by checking tolerance-based numerical results.
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

  /// @brief Verifies tetrahedron face 0 for geometry polytope project.
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

  /// @brief Verifies tetrahedron vertex 0 for geometry polytope project by checking tolerance-based numerical results.
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

  // ---- Pyramid projections ----

  /// @brief Verifies pyramid cell for geometry polytope project by checking tolerance-based numerical results.
  TEST(Geometry_PolytopeProject, Pyramid_Cell)
  {
    Polytope::Project proj(Polytope::Type::Pyramid);
    Math::SpatialPoint rc{0.2, 0.3, 0.4};
    Math::SpatialPoint out(3);
    proj.cell(out, rc);
    EXPECT_NEAR(out(0), 0.2, 1e-14);
    EXPECT_NEAR(out(1), 0.3, 1e-14);
    EXPECT_NEAR(out(2), 0.4, 1e-14);

    Math::SpatialPoint outside{2.0, 2.0, 0.5};
    proj.cell(out, outside);
    EXPECT_GE(out(0), -1e-14);
    EXPECT_GE(out(1), -1e-14);
    EXPECT_GE(out(2), -1e-14);
    EXPECT_LE(out(0), 1.0 - out(2) + 1e-14);
    EXPECT_LE(out(1), 1.0 - out(2) + 1e-14);
  }

  /// @brief Verifies pyramid faces for geometry polytope project.
  TEST(Geometry_PolytopeProject, Pyramid_Faces)
  {
    Polytope::Project proj(Polytope::Type::Pyramid);
    Math::SpatialPoint rc{0.2, 0.3, 0.4};
    Math::SpatialPoint out(3);

    for (size_t f = 0; f < 5; ++f)
    {
      proj.face(f, out, rc);
      EXPECT_GE(out(0), -1e-14);
      EXPECT_GE(out(1), -1e-14);
      EXPECT_GE(out(2), -1e-14);
      EXPECT_LE(out(0), 1.0 - out(2) + 1e-14);
      EXPECT_LE(out(1), 1.0 - out(2) + 1e-14);
    }
  }

  // ---- Hexahedron projections ----

  /// @brief Verifies hexahedron cell for geometry polytope project by checking tolerance-based numerical results.
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

  /// @brief Verifies hexahedron vertex 0 for geometry polytope project by checking tolerance-based numerical results.
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

  /// @brief Verifies triangle boundary for geometry polytope project.
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

namespace Rodin::Tests::Unit
{
  /// @brief The point face projection resets reused output storage to dimension zero.
  TEST(Geometry_PolytopeProject, PointFaceProjectionClearsReusedStorage)
  {
    const Polytope::Project project(Polytope::Type::Point);
    const Math::SpatialPoint input;
    Math::SpatialPoint output{2, -3, 4};
    project.face(0, output, input);
    EXPECT_EQ(output.size(), 0);
    project.face(0, output, input);
    EXPECT_EQ(output.size(), 0);
    project.face(0, output, output);
    EXPECT_EQ(output.size(), 0);
  }

  /// @brief Reference vertices are fixed points of cell and boundary projections.
  TEST(Geometry_PolytopeProject, EveryReferenceVertexIsFixed)
  {
    constexpr Real tolerance = 1e-13;
    for (auto geometry : Polytope::Types)
    {
      const Polytope::Traits traits(geometry);
      const Polytope::Project project(geometry);
      for (size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
      {
        SCOPED_TRACE(::testing::Message() << static_cast<int>(geometry) << ":" << vertex);
        const auto& expected = traits.getVertex(vertex);
        ASSERT_EQ(expected.size(), traits.getDimension());
        Math::SpatialPoint cell, boundary;
        project.cell(cell, expected);
        project.boundary(boundary, expected);
        ASSERT_EQ(cell.size(), expected.size());
        ASSERT_EQ(boundary.size(), expected.size());
        for (size_t axis = 0; axis < expected.size(); ++axis)
        {
          EXPECT_NEAR(cell(axis), expected(axis), tolerance);
          EXPECT_NEAR(boundary(axis), expected(axis), tolerance);
        }
      }
    }
  }

  /// @brief Projection fixes feasible points and returns feasible, idempotent results.
  TEST(Geometry_PolytopeProject, ExteriorPointsReturnToReferenceCells)
  {
    constexpr Real tolerance = 1e-12;
    for (auto geometry : Polytope::Types)
    {
      const Polytope::Traits traits(geometry);
      const size_t dimension = traits.getDimension();
      if (dimension == 0)
        continue;
      const Polytope::Project project(geometry);
      const auto& halfSpace = traits.getHalfSpace();
      for (Real first : {-2.0, 0.2, 0.8, 2.0})
      {
        for (Real second : {-1.0, 0.3, 1.5})
        {
          for (Real third : {-1.0, 0.4, 2.0})
          {
            SCOPED_TRACE(::testing::Message() << static_cast<int>(geometry) << ":"
                                              << first << "," << second << "," << third);
            Math::SpatialPoint input(dimension);
            input(0) = first;
            if (dimension >= 2)
              input(1) = second;
            if (dimension == 3)
              input(2) = third;
            Math::SpatialPoint cell, boundary, repeated;
            project.cell(cell, input);
            project.boundary(boundary, input);
            for (const auto& result : {cell, boundary})
            {
              const Math::Vector<Real> residual =
                halfSpace.matrix * result.getData().head(dimension) - halfSpace.vector;
              EXPECT_LE(residual.maxCoeff(), tolerance);
            }
            project.cell(repeated, cell);
            EXPECT_LE((repeated - cell).norm(), tolerance);
            project.boundary(repeated, boundary);
            EXPECT_LE((repeated - boundary).norm(), tolerance);
            const Math::Vector<Real> boundaryResidual =
              halfSpace.matrix * boundary.getData().head(dimension) - halfSpace.vector;
            EXPECT_LE(boundaryResidual.cwiseAbs().minCoeff(), tolerance);
          }
        }
      }
    }
  }
}
