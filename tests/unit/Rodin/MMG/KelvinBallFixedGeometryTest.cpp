/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>

#include "examples/ShapeOptimization/KelvinBall/Sphere.h"
#include "examples/ShapeOptimization/KelvinBall/LaplacianSmoothing.h"

namespace KelvinBall
{
  TEST(KelvinBallBackground, EmptyCellAttributesExportAsOneMMGReference)
  {
    MMG::Mesh mesh(Mesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3}));
    mesh.getConnectivity().compute(2, 3);
    for (auto cell = mesh.getCell(); cell; ++cell)
      mesh.setAttribute({3, cell->getIndex()}, {});
    for (auto face = mesh.getBoundary(); face; ++face)
      mesh.setAttribute({2, face->getIndex()}, Outer);
    auto* native = MMG::MMG5::rodinToMesh(mesh);
    ASSERT_NE(native, nullptr);
    EXPECT_EQ(native->ne, mesh.getCellCount());
    for (int i = 1; i <= native->ne; ++i)
      EXPECT_EQ(native->tetra[i].ref, 0);
    for (int i = 1; i <= native->nt; ++i)
      EXPECT_EQ(native->tria[i].ref, Outer);
    MMG::MMG5::destroyMesh(native);
  }

  TEST(
    KelvinBallBackground, AdaptationRemovesMaterialsAndOldInterfaceWithoutMutatingSource)
  {
    Configuration configuration;
    configuration.points = 5;
    configuration.hminFactor = Real(0.5);
    configuration.hmaxFactor = Real(2);
    configuration.finalize();
    Sphere sphere(configuration);
    auto initial = sphere.prepareSWIFTBackground();
    auto& source = initial.mesh;
    source.getConnectivity().compute(2, 3);
    for (auto cell = source.getCell(); cell; ++cell)
    {
      source.setAttribute(
        {3, cell->getIndex()}, centroid(source, *cell).norm() < 1 ? Obstacle : Fluid);
    }
    source.trace({{{Obstacle, Fluid}, Gamma}});
    source.setRequiredTetrahedron(0);
    const MMG::Mesh original(source);
    P1<Real, Mesh> sizeSpace(source);
    MMG::RealGridFunction sizes(sizeSpace);
    sizes = RealFunction(configuration.getGridSpacing());
    const auto adapted = sphere.adaptBackground(source, sizes);
    EXPECT_FALSE(adapted.getRequiredTetrahedra().contains(0));
    for (auto cell = adapted.getCell(); cell; ++cell)
      EXPECT_FALSE(cell->getAttribute());
    for (auto face = adapted.getFace(); face; ++face)
      EXPECT_NE(face->getAttribute(), Gamma);
    EXPECT_EQ(source.getCellCount(), original.getCellCount());
    EXPECT_EQ(source.getRequiredTetrahedra(), original.getRequiredTetrahedra());
    for (auto cell = source.getCell(); cell; ++cell)
    {
      EXPECT_EQ(
        cell->getAttribute(), original.getPolytope(3, cell->getIndex())->getAttribute());
    }
    for (Index i = 0; i < source.getVertexCount(); ++i)
    {
      EXPECT_EQ(
        (source.getVertexCoordinates(i) - original.getVertexCoordinates(i)).norm(), 0);
    }
  }

  TEST(KelvinBallMeshSmoothing, ImprovesInteriorVertexWithoutMovingInterfaceOrWall)
  {
    Mesh mesh = Mesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    mesh.getConnectivity().compute(2, 3);
    for (auto face = mesh.getBoundary(); face; ++face)
    {
      mesh.setAttribute(
        {2, face->getIndex()}, centroid(mesh, *face)(0) == 0 ? Gamma : Outer);
    }
    Index interior = mesh.getVertexCount();
    Math::SpatialPoint center(3);
    center.setConstant(1);
    for (Index i = 0; i < mesh.getVertexCount(); ++i)
    {
      if ((mesh.getVertexCoordinates(i) - center).norm() == 0)
        interior = i;
    }
    ASSERT_LT(interior, mesh.getVertexCount());
    auto perturbed = center;
    perturbed(0) += Real(0.3);
    mesh.setVertexCoordinates(interior, perturbed);
    const Mesh original(mesh);
    // Populate mapped caches before smoothing to test their invalidation too.
    Real volumeBefore = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
      volumeBefore += cell->getMeasure();
    EXPECT_GT(LaplacianSmoothing().execute(mesh), 0);
    EXPECT_LT((mesh.getVertexCoordinates(interior) - center).norm(), Real(0.3));
    for (Index i = 0; i < mesh.getVertexCount(); ++i)
    {
      if (i != interior)
      {
        EXPECT_EQ(
          (mesh.getVertexCoordinates(i) - original.getVertexCoordinates(i)).norm(), 0);
      }
    }
    Real volumeAfter = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
      volumeAfter += cell->getMeasure();
    EXPECT_NEAR(volumeAfter, volumeBefore, 1e-12);
    EXPECT_EQ(mesh.getCellCount(), original.getCellCount());
  }

  TEST(KelvinBallMeshSmoothing, RestrictsMotionToEveryIncidentCutPlane)
  {
    Mesh mesh = Mesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    mesh.getConnectivity().compute(2, 3);
    for (auto face = mesh.getBoundary(); face; ++face)
    {
      const auto x = centroid(mesh, *face);
      mesh.setAttribute(
        {2, face->getIndex()}, x(0) == 0 ? SigmaXYPlus : (x(1) == 0 ? SigmaPlus : Outer));
    }
    const Mesh original(mesh);
    LaplacianSmoothing().execute(mesh);
    for (auto face = mesh.getBoundary(); face; ++face)
    {
      for (const Index vertex : face->getVertices())
      {
        const auto delta =
          mesh.getVertexCoordinates(vertex) - original.getVertexCoordinates(vertex);
        if (face->getAttribute() == SigmaXYPlus)
          EXPECT_NEAR(delta(0) - delta(1), 0, 1e-14);
        else if (face->getAttribute() == SigmaPlus)
          EXPECT_NEAR(delta(1) - delta(2), 0, 1e-14);
        else
          EXPECT_EQ(delta.norm(), 0);
      }
    }
  }

  TEST(KelvinBallMeshOperations, MaterialTracePreservesChamberBoundaryLabels)
  {
    Mesh mesh = Mesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    mesh.getConnectivity().compute(2, 3);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Geometry::Point point(
        *cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      mesh.setAttribute({3, cell->getIndex()}, point.x() < 1 ? Obstacle : Fluid);
    }
    for (auto face = mesh.getBoundary(); face; ++face)
      mesh.setAttribute({2, face->getIndex()}, SigmaPlus);

    mesh.trace({{{Obstacle, Fluid}, Gamma}});
    size_t interface = 0;
    for (auto face = mesh.getFace(); face; ++face)
    {
      const auto& cells = mesh.getConnectivity().getIncidence({2, 3}, face->getIndex());
      if (cells.size() == 1)
        EXPECT_EQ(face->getAttribute(), SigmaPlus);
      else if (mesh.getPolytope(3, cells[0])->getAttribute() !=
        mesh.getPolytope(3, cells[1])->getAttribute())
      {
        EXPECT_EQ(face->getAttribute(), Gamma);
        ++interface;
      }
      else
        EXPECT_FALSE(face->getAttribute());
    }
    EXPECT_GT(interface, 0);
  }

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
    {
      for (size_t i = 0; i < 3; ++i)
        EXPECT_NEAR(mesh.getVertexCoordinates(v)(i), points[v][i], 1e-14);
    }
    sphere.projectFixedGeometry(mesh);
    for (Index v = 0; v < 4; ++v)
    {
      for (size_t i = 0; i < 3; ++i)
        EXPECT_NEAR(mesh.getVertexCoordinates(v)(i), points[v][i], 1e-14);
    }
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
    {
      if (face->getAttribute() == Outer)
      {
        for (const Index vertex : face->getVertices())
        {
          EXPECT_NEAR(mesh.getVertexCoordinates(vertex)(0), 2, 1e-14);
        }
      }
    }
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
    {
      EXPECT_EQ(
        (mesh.getVertexCoordinates(vertex) - original.getVertexCoordinates(vertex))
          .norm(),
        0);
    }
    for (auto face = mesh.getFace(); face; ++face)
    {
      EXPECT_EQ(
        face->getAttribute(), original.getPolytope(2, face->getIndex())->getAttribute());
    }
  }

  TEST_F(FixedGeometryTest, RequiredCutAdaptationPreservesCutTriangulations)
  {
    auto mesh = makeTetrahedron({{{0, 0, 0}, {2, 0, 0}, {2, 2, 2}, {2, 2, -2}}});
    labelOpposite(mesh, 0, Outer);
    labelOpposite(mesh, 1, Gamma);
    labelOpposite(mesh, 2, SigmaMinus);
    labelOpposite(mesh, 3, SigmaPlus);
    const auto cutTriangles = [](const MMG::Mesh& input) {
      FlatSet<std::array<long long, 9>> triangles;
      // Compare coordinate sets, since MMG may re-index the mesh.
      constexpr Real coordinateTolerance = 1e-12;
      for (auto face = input.getFace(); face; ++face)
      {
        if (face->getAttribute() != SigmaPlus && face->getAttribute() != SigmaMinus)
          continue;
        std::array<std::array<long long, 3>, 3> vertices;
        for (size_t i = 0; i < 3; ++i)
        {
          const auto& x = input.getVertexCoordinates(face->getVertices()[i]);
          for (size_t j = 0; j < 3; ++j)
            vertices[i][j] = std::llround(x(j) / coordinateTolerance);
        }
        std::sort(vertices.begin(), vertices.end());
        std::array<long long, 9> key;
        for (size_t i = 0; i < 3; ++i)
        {
          for (size_t j = 0; j < 3; ++j)
            key[3 * i + j] = vertices[i][j];
        }
        triangles.insert(key);
      }
      return triangles;
    };
    const auto original = cutTriangles(mesh);
    ASSERT_EQ(original.size(), 2);
    Configuration configuration;
    configuration.hmin = 0.1;
    configuration.hmax = 1;
    ASSERT_NO_THROW(Sphere(configuration).adapt(mesh, 1, true));
    EXPECT_EQ(cutTriangles(mesh), original);
    for (auto face = mesh.getFace(); face; ++face)
    {
      if (face->getAttribute() == SigmaPlus || face->getAttribute() == SigmaMinus)
      {
        EXPECT_TRUE(mesh.getRequiredTriangles().contains(face->getIndex()));
      }
    }
  }

  TEST_F(FixedGeometryTest, FailedRequiredCutRetryRetainsFittedMeshAndFeatureTags)
  {
    auto mesh =
      makeTetrahedron({{{0, 0.01, 0}, {1, 0, 0}, {0, 1, -1}, {0.2, 0.5, -0.5}}});
    labelOpposite(mesh, 0, Gamma);
    labelOpposite(mesh, 3, SigmaMinus);
    const MMG::Mesh original = mesh;
    Configuration configuration;
    configuration.hmin = 0.1;
    configuration.hmax = 1;
    Sphere sphere(configuration);
    EXPECT_THROW(sphere.adapt(mesh, 1), std::runtime_error);
    EXPECT_THROW(sphere.adapt(mesh, 1, true), std::runtime_error);
    EXPECT_EQ(mesh.getCellCount(), original.getCellCount());
    EXPECT_EQ(mesh.getVertexCount(), original.getVertexCount());
    EXPECT_EQ(mesh.getCorners(), original.getCorners());
    EXPECT_EQ(mesh.getRidges(), original.getRidges());
    EXPECT_EQ(mesh.getRequiredTriangles(), original.getRequiredTriangles());
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      EXPECT_EQ(
        (mesh.getVertexCoordinates(vertex) - original.getVertexCoordinates(vertex))
          .norm(),
        0);
    }
    for (auto face = mesh.getFace(); face; ++face)
    {
      EXPECT_EQ(
        face->getAttribute(), original.getPolytope(2, face->getIndex())->getAttribute());
    }
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
