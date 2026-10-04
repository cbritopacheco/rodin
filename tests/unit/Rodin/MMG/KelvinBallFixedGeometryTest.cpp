/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "examples/ShapeOptimization/KelvinBall/Sphere.h"

namespace KelvinBall
{
  class FixedGeometryTest : public ::testing::Test
  {
    protected:
      MMG::Mesh makeTetrahedron(const std::array<std::array<Real, 3>, 4>& points)
      {
        Mesh::Builder builder;
        builder.initialize(3).nodes(4);
        for (const auto& x : points)
        {
          Math::SpatialPoint point(3);
          for (size_t i = 0; i < 3; ++i)
            point(i) = x[i];
          builder.vertex(point);
        }
        Index cell;
        IndexArray vertices(4);
        vertices << 0, 1, 2, 3;
        builder.polytope(Polytope::Type::Tetrahedron, vertices, cell);
        builder.attribute({3, cell}, Fluid);
        MMG::Mesh mesh(builder.finalize());
        mesh.getConnectivity().discover(3, 2);
        mesh.getConnectivity().compute(2, 3);
        return mesh;
      }

      void labelOpposite(MMG::Mesh& mesh, Index vertex, Attribute attribute)
      {
        for (auto face = mesh.getFace(); face; ++face)
        {
          bool contains = false;
          for (const Index v : face->getVertices())
            contains = contains || v == vertex;
          if (!contains)
            mesh.setAttribute({2, face->getIndex()}, attribute);
        }
      }
  };

  TEST_F(FixedGeometryTest, SharedCutsAndWallProjectToTheirIntersection)
  {
    const std::array<std::array<Real, 3>, 4> points{
      {{0, 0, 0}, {2, 0, 0}, {2, 2, 2}, {2, 2, -2}}};
    auto mesh = makeTetrahedron(points);
    labelOpposite(mesh, 0, Outer);
    labelOpposite(mesh, 1, SigmaXYPlus);
    labelOpposite(mesh, 2, SigmaMinus);
    labelOpposite(mesh, 3, SigmaPlus);
    for (Index v = 0; v < 4; ++v)
    {
      auto x = mesh.getVertexCoordinates(v);
      x(0) += 0.003;
      x(1) -= 0.002;
      x(2) += 0.001;
      mesh.setVertexCoordinates(v, x);
    }
    Configuration configuration;
    configuration.outerRadius = 2;
    Sphere sphere(configuration);
    sphere.projectFixedGeometry(mesh);
    for (Index v = 0; v < 4; ++v)
      for (size_t i = 0; i < 3; ++i)
        EXPECT_NEAR(mesh.getVertexCoordinates(v)(i), points[v][i], 1e-14);
    sphere.projectFixedGeometry(mesh);
    for (Index v = 0; v < 4; ++v)
      for (size_t i = 0; i < 3; ++i)
        EXPECT_NEAR(mesh.getVertexCoordinates(v)(i), points[v][i], 1e-14);
  }

  TEST_F(FixedGeometryTest, SingleCutProjectionLeavesInteriorVertexUnchanged)
  {
    auto mesh = makeTetrahedron({{{0, 0.01, 0}, {1, 0, 0}, {0, 1, -1}, {0.2, 0.5, 0}}});
    labelOpposite(mesh, 3, SigmaMinus);
    const auto interior = mesh.getVertexCoordinates(3);
    Configuration configuration;
    Sphere(configuration).projectFixedGeometry(mesh);
    EXPECT_NEAR(mesh.getVertexCoordinates(0)(1), 0.005, 1e-14);
    EXPECT_NEAR(mesh.getVertexCoordinates(0)(2), -0.005, 1e-14);
    EXPECT_EQ((mesh.getVertexCoordinates(3) - interior).norm(), 0);
    for (Index v = 0; v < 3; ++v)
    {
      const auto x = mesh.getVertexCoordinates(v);
      EXPECT_NEAR(x(1) + x(2), 0, 1e-14);
    }
  }

  TEST_F(FixedGeometryTest, SuccessfulAdaptationCommitsCandidate)
  {
    auto mesh = makeTetrahedron({{{0, 0, 0}, {2, 0, 0}, {2, 2, 2}, {2, 2, -2}}});
    labelOpposite(mesh, 0, Outer);
    labelOpposite(mesh, 1, Gamma);
    labelOpposite(mesh, 2, SigmaMinus);
    labelOpposite(mesh, 3, SigmaPlus);
    EXPECT_TRUE(mesh.getRequiredTriangles().empty());
    Configuration configuration;
    configuration.hmin = 0.1;
    configuration.hmax = 1;
    EXPECT_NO_THROW(Sphere(configuration).adapt(mesh, 1));
    EXPECT_FALSE(mesh.getRequiredTriangles().empty());
    EXPECT_GT(mesh.getCellCount(), 0);
    for (auto face = mesh.getFace(); face; ++face)
      if (face->getAttribute() == Outer)
        for (const Index vertex : face->getVertices())
          EXPECT_NEAR(mesh.getVertexCoordinates(vertex)(0), 2, 1e-14);
  }

  TEST_F(FixedGeometryTest, FailedAdaptationRetainsFittedMeshAndFeatureTags)
  {
    auto mesh =
      makeTetrahedron({{{0, 0.01, 0}, {1, 0, 0}, {0, 1, -1}, {0.2, 0.5, -0.5}}});
    labelOpposite(mesh, 0, Gamma);
    labelOpposite(mesh, 3, SigmaMinus);
    const MMG::Mesh original = mesh;
    Configuration configuration;
    // Restoring the labelled cut collapses this thin tetrahedron. All MMG
    // sizes and tolerances remain valid, isolating post-adaptation rejection.
    configuration.hmin = 0.1;
    configuration.hmax = 1;
    try
    {
      Sphere(configuration).adapt(mesh, 1);
      FAIL() << "Expected fixed-boundary projection to reject the adapted mesh.";
    }
    catch (const std::runtime_error& error)
    {
      EXPECT_NE(
        std::string(error.what()).find("Fixed-boundary projection"), std::string::npos);
    }
    EXPECT_EQ(mesh.getCellCount(), original.getCellCount());
    EXPECT_EQ(mesh.getVertexCount(), original.getVertexCount());
    EXPECT_EQ(mesh.getCorners(), original.getCorners());
    EXPECT_EQ(mesh.getRidges(), original.getRidges());
    EXPECT_EQ(mesh.getRequiredTriangles(), original.getRequiredTriangles());
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      EXPECT_EQ(
        (mesh.getVertexCoordinates(vertex) - original.getVertexCoordinates(vertex))
          .norm(),
        0);
    for (auto face = mesh.getFace(); face; ++face)
      EXPECT_EQ(
        face->getAttribute(), original.getPolytope(2, face->getIndex())->getAttribute());
  }

  TEST_F(FixedGeometryTest, CollapsedProjectionIsRejectedWithoutMutation)
  {
    auto mesh =
      makeTetrahedron({{{0, 0.01, 0}, {1, 0, 0}, {0, 1, -1}, {0.2, 0.5, -0.5}}});
    labelOpposite(mesh, 3, SigmaMinus);
    const auto original = mesh.getVertexCoordinates(0);
    Configuration configuration;
    EXPECT_THROW(Sphere(configuration).projectFixedGeometry(mesh), std::runtime_error);
    EXPECT_EQ((mesh.getVertexCoordinates(0) - original).norm(), 0);
  }
}
