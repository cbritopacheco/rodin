/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Boundary extraction and submesh ownership contracts. */
#include <gtest/gtest.h>
#include <algorithm>
#include <utility>
#include <vector>
#include "Rodin/Geometry.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  LocalMesh parentMesh(Polytope::Type geometry)
  {
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {3})
      : dimension == 2         ? LocalMesh::UniformGrid(geometry, {2, 2})
                               : LocalMesh::UniformGrid(geometry, {2, 2, 2});
    for (size_t d = 1; d <= dimension; ++d)
      for (size_t lower = 0; lower < d; ++lower)
        mesh.getConnectivity().compute(d, lower);
    mesh.getConnectivity().compute(dimension - 1, dimension);
    for (size_t d = 0; d <= dimension; ++d)
      for (Index i = 0; i < mesh.getPolytopeCount(d); ++i)
        mesh.setAttribute({d, i}, Attribute(100 * d + i));
    return mesh;
  }

  void checkBoundaryMaps(const SubMesh<Context::Local>& skin, const LocalMesh& parent)
  {
    const size_t dimension = parent.getDimension();
    EXPECT_EQ(skin.getDimension(), dimension - 1);
    EXPECT_EQ(skin.getSpaceDimension(), dimension);
    EXPECT_EQ(&skin.getParent(), &parent);
    ASSERT_EQ(skin.getAncestors().size(), 1);
    EXPECT_EQ(&skin.getAncestors().front().get(), &parent);
    size_t boundaryCount = 0;
    for (auto face = parent.getBoundary(); face; ++face)
      ++boundaryCount;
    EXPECT_EQ(skin.getCellCount(), boundaryCount);
    for (size_t d = 0; d < dimension; ++d)
    {
      const auto& map = skin.getPolytopeMap(d);
      ASSERT_EQ(map.left.size(), skin.getPolytopeCount(d));
      for (Index i = 0; i < skin.getPolytopeCount(d); ++i)
      {
        const Index original = map.left[i];
        EXPECT_EQ(map.right.at(original), i);
        EXPECT_EQ(skin.getAttribute(d, i), parent.getAttribute(d, original));
        const auto polytope = skin.getPolytope(d, i);
        const auto parentPolytope = parent.getPolytope(d, original);
        EXPECT_EQ(polytope->getGeometry(), parentPolytope->getGeometry());
        if (d == dimension - 1)
        {
          EXPECT_TRUE(parent.getFace(original)->isBoundary());
        }
      }
    }
    for (Index vertex = 0; vertex < skin.getVertexCount(); ++vertex)
    {
      const auto& actual = skin.getVertexCoordinates(vertex);
      const auto& expected =
        parent.getVertexCoordinates(skin.getPolytopeMap(0).left[vertex]);
      EXPECT_NEAR((actual - expected).norm(), 0, 1e-13);
    }
    const Real measure = dimension == 1 ? Real(2) : dimension == 2 ? Real(4) : Real(6);
    EXPECT_NEAR(skin.getMeasure(dimension - 1), measure, 1e-12);
  }

  SubMesh<Context::Local> independentBoundaryCopy(const LocalMesh& parent)
  {
    auto original = parent.skin();
    SubMesh<Context::Local> copied(original);
    SubMesh<Context::Local> moved(std::move(copied));
    return moved;
  }
}

/// @brief A one-dimensional skin contains only the two endpoint vertices.
TEST(Geometry_SkinLifecycle, SegmentBoundarySurvivesCopies)
{
  const auto mesh = parentMesh(Polytope::Type::Segment);
  const auto skin = independentBoundaryCopy(mesh);
  checkBoundaryMaps(skin, mesh);
  ASSERT_EQ(skin.getVertexCount(), 2);
  std::vector<Real> endpoints;
  for (Index vertex = 0; vertex < skin.getVertexCount(); ++vertex)
    endpoints.push_back(skin.getVertexCoordinates(vertex)[0]);
  std::sort(endpoints.begin(), endpoints.end());
  EXPECT_EQ(endpoints, (std::vector<Real>{0, 2}));
}

/// @brief Higher-dimensional skins preserve boundary geometry, attributes and closure maps.
TEST(Geometry_SkinLifecycle, SurfaceBoundariesSurviveCopies)
{
  for (const auto geometry :
    {Polytope::Type::Triangle, Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(int(geometry));
    const auto mesh = parentMesh(geometry);
    const auto skin = independentBoundaryCopy(mesh);
    checkBoundaryMaps(skin, mesh);
  }
}

/// @brief Nested selections retain both ancestors and reject restriction from unrelated meshes.
TEST(Geometry_SkinLifecycle, NestedBoundarySelectionsKeepAncestry)
{
  const auto mesh = parentMesh(Polytope::Type::Hexahedron);
  const auto skin = independentBoundaryCopy(mesh);
  const Attribute selected = *skin.getAttribute(skin.getDimension(), 0);
  const auto nested = skin.keep(selected);
  EXPECT_EQ(&nested.getParent(), &skin);
  ASSERT_EQ(nested.getAncestors().size(), 2);
  EXPECT_EQ(&nested.getAncestors()[0].get(), &skin);
  EXPECT_EQ(&nested.getAncestors()[1].get(), &mesh);
  ASSERT_EQ(nested.getCellCount(), 1);
  const Index rootFace = skin.getPolytopeMap(2).left[nested.getPolytopeMap(2).left[0]];
  const auto face = mesh.getPolytope(2, rootFace);
  const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
  const auto restricted = nested.restriction(point);
  ASSERT_TRUE(restricted);
  EXPECT_EQ(&restricted->getPolytope().getMesh(), &nested);
  EXPECT_NEAR(
    (restricted->getPhysicalCoordinates() - point.getPhysicalCoordinates()).norm(), 0,
    1e-13);
  const auto unrelated = parentMesh(Polytope::Type::Hexahedron);
  const auto otherFace = unrelated.getFace();
  const Point other(*otherFace, Polytope::Traits(otherFace->getGeometry()).getCentroid());
  EXPECT_FALSE(nested.restriction(other));
}

namespace
{
  LocalMesh pointCloud()
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(1)
                  .nodes(3)
                  .vertex({0})
                  .vertex({1})
                  .vertex({2})
                  .finalize();
    mesh.setAttribute({0, 0}, 7);
    mesh.setAttribute({0, 1}, 9);
    mesh.setAttribute({0, 2}, 7);
    return mesh;
  }
}

/// @brief Keeping a labeled point cloud selects vertices without nonexistent intermediate grades.
TEST(Geometry_SkinLifecycle, PointCloudKeepSkipsIntermediateGrades)
{
  const auto mesh = pointCloud();
  const auto selected = mesh.keep(Attribute(7));
  EXPECT_EQ(selected.getDimension(), 0);
  ASSERT_EQ(selected.getVertexCount(), 2);
  EXPECT_EQ(&selected.getParent(), &mesh);
  for (Index i = 0; i < 2; ++i)
  {
    EXPECT_EQ(selected.getVertexCoordinates(i)[0], Real(2 * i));
    EXPECT_EQ(selected.getPolytopeMap(0).left[i], 2 * i);
    EXPECT_EQ(selected.getAttribute(0, i), 7);
  }
}

/// @brief Trimming a labeled point cloud preserves the complementary vertex and parent map.
TEST(Geometry_SkinLifecycle, PointCloudTrimSkipsIntermediateGrades)
{
  const auto mesh = pointCloud();
  const auto selected = mesh.trim(Attribute(7));
  EXPECT_EQ(selected.getDimension(), 0);
  ASSERT_EQ(selected.getVertexCount(), 1);
  EXPECT_EQ(&selected.getParent(), &mesh);
  EXPECT_EQ(selected.getVertexCoordinates(0)[0], 1);
  EXPECT_EQ(selected.getPolytopeMap(0).left[0], 1);
  EXPECT_EQ(selected.getAttribute(0, 0), 9);
}
