/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Orientation and trace contracts for codimension-one normals. */
#include <gtest/gtest.h>
#include "Rodin/Geometry/ParametricTransformation.h"
#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/Variational/FaceNormal.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  void checkOrientation(
    const FaceNormal& normal, const Point& point, const Polytope& cell)
  {
    const auto value = normal.getValue(point);
    const Point center(cell, Polytope::Traits(cell.getGeometry()).getCentroid());
    EXPECT_NEAR(value.norm(), 1, 1e-12);
    EXPECT_LT(value.dot(center.getCoordinates() - point.getCoordinates()), 0);
    for (Eigen::Index j = 0; j < point.getJacobian().cols(); ++j)
      EXPECT_NEAR(value.dot(point.getJacobian().col(j)), 0, 1e-12);
  }
}

TEST(Rodin_Variational_FaceNormal, BoundaryOrientationAndInteriorTraceAcrossGeometries)
{
  for (const auto geometry :
    {Polytope::Type::Triangle, Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 3));
    mesh.getConnectivity().compute(dimension - 1, dimension);
    for (auto cell = mesh.getPolytope(dimension); cell; ++cell)
      mesh.setAttribute({dimension, cell->getIndex()}, cell->getIndex() + 1);
    const auto& adjacent = mesh.getConnectivity().getIncidence(dimension - 1, dimension);
    FaceNormal normal(mesh);
    EXPECT_EQ(normal.getDimension(), dimension);
    size_t boundary = 0, interior = 0;
    for (auto face = mesh.getPolytope(dimension - 1); face; ++face)
    {
      const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
      const auto& cells = adjacent.at(face->getIndex());
      EXPECT_EQ(normal.getOrder(*face), 0);
      if (cells.size() == 1)
      {
        checkOrientation(normal, point, *mesh.getPolytope(dimension, cells[0]));
        std::unique_ptr<FaceNormal> clone(normal.copy());
        FaceNormal moved(std::move(*clone));
        EXPECT_LT((normal(point) - moved(point)).norm(), 1e-12);
        ++boundary;
      }
      else
      {
        ASSERT_EQ(cells.size(), 2);
        if (interior == 0)
        {
          EXPECT_THROW(normal.getValue(point), Alert::Exception);
        }
        auto first = normal;
        first.traceOf(cells[0] + 1);
        auto second = normal;
        second.traceOf(cells[1] + 1);
        checkOrientation(first, point, *mesh.getPolytope(dimension, cells[0]));
        checkOrientation(second, point, *mesh.getPolytope(dimension, cells[1]));
        EXPECT_LT((first(point) + second(point)).norm(), 1e-12);
        if (interior == 0)
        {
          auto invalid = normal;
          invalid.traceOf(999999);
          EXPECT_THROW(invalid.getValue(point), Alert::Exception);
        }
        ++interior;
      }
    }
    EXPECT_GT(boundary, 0);
    EXPECT_GT(interior, 0);
  }
}

TEST(Rodin_Variational_FaceNormal, UnsupportedOneDimensionalJacobianRaises)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  mesh.getConnectivity().compute(0, 1);
  const auto vertex = mesh.getPolytope(0, 0);
  const Point point(*vertex, Math::SpatialPoint(0));
  FaceNormal normal(mesh);
  EXPECT_THROW(normal.getValue(point), Alert::Exception);
}

/// @brief Embedded and folded surface traces use the selected cell's tangent plane.
TEST(Rodin_Variational_FaceNormal, EmbeddedSurfaceConormalsFollowEachIncidentCell)
{
  using Vector = Math::SpatialVector<Real>;
  for (const bool folded : {false, true})
  {
    SCOPED_TRACE(folded);
    LocalMesh::Builder builder;
    builder.initialize(3)
      .nodes(4)
      .vertex({0, 0, 0})
      .vertex({0, 0, 1})
      .vertex({0.6, 0.8, 0})
      .vertex(folded ? Vector{-0.8, 0.6, 0} : Vector{-0.6, -0.8, 0})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 0, 3});
    auto mesh = builder.finalize();
    mesh.getConnectivity().compute(1, 2);
    for (auto cell = mesh.getCell(); cell; ++cell)
      mesh.setAttribute({2, cell->getIndex()}, cell->getIndex() + 1);
    const auto& adjacency = mesh.getConnectivity().getIncidence(1, 2);
    FaceNormal normal(mesh);
    EXPECT_EQ(normal.getDimension(), 3);
    size_t checked = 0;
    for (auto edge = mesh.getFace(); edge; ++edge)
    {
      const auto& vertices = edge->getVertices();
      const Vector tangent =
        mesh.getVertexCoordinates(vertices[1]) - mesh.getVertexCoordinates(vertices[0]);
      for (const Real sample : {0.25, 0.75})
      {
        const Point point(*edge, Math::SpatialPoint{sample});
        for (const Index index : adjacency.at(edge->getIndex()))
        {
          auto traced = normal;
          traced.traceOf(index + 1);
          const auto cell = mesh.getCell(index);
          const Point center(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
          const Vector inward = center.getCoordinates() - point.getCoordinates();
          const Vector perpendicular =
            inward - tangent * (inward.dot(tangent) / tangent.squaredNorm());
          const Vector expected = -perpendicular.normalized();
          const auto actual = traced.getValue(point);
          EXPECT_LT((actual - expected).norm(), 1e-12);
          checkOrientation(traced, point, *cell);
          if (adjacency.at(edge->getIndex()).size() == 1)
          {
            EXPECT_LT((normal.getValue(point) - expected).norm(), 1e-12);
          }
          ++checked;
        }
      }
    }
    EXPECT_EQ(checked, 12);
  }
}

/// @brief Curved-surface conormals use the tangent plane at the evaluated face point.
TEST(Rodin_Variational_FaceNormal, CurvedSurfaceConormalsFollowLocalTangentPlane)
{
  LocalMesh::Builder builder;
  builder.initialize(3)
    .nodes(3)
    .vertex({0, 0, 0})
    .vertex({1, 0, 0.5})
    .vertex({0, 1, 1})
    .polytope(Polytope::Type::Triangle, {0, 1, 2});
  auto mesh = builder.finalize();
  mesh.getConnectivity().compute(1, 2);
  for (size_t dimension : {size_t(1), size_t(2)})
  {
    for (auto polytope = mesh.getPolytope(dimension); polytope; ++polytope)
    {
      RealH1Element<2> element(polytope->getGeometry());
      PointCloud controls(3, element.getCount());
      for (size_t node = 0; node < element.getCount(); ++node)
      {
        const Point affine(*polytope, element.getNode(node));
        const auto coordinates = affine.getCoordinates();
        controls(0, node) = coordinates.x();
        controls(1, node) = coordinates.y();
        controls(2, node) =
          0.5 * coordinates.x() * coordinates.x() + coordinates.y() * coordinates.y();
      }
      ParametricTransformation transformation(controls, element);
      mesh.setPolytopeTransformation(
        {dimension, polytope->getIndex()}, transformation.copy());
    }
  }
  FaceNormal normal(mesh);
  const auto cell = mesh.getCell(0);
  for (auto edge = mesh.getFace(); edge; ++edge)
  {
    const auto& vertices = edge->getVertices();
    const auto chord =
      mesh.getVertexCoordinates(vertices[1]) - mesh.getVertexCoordinates(vertices[0]);
    for (const Real sample : {0.2, 0.5, 0.8})
    {
      const Point point(*edge, Math::SpatialPoint{sample});
      const auto physical = point.getCoordinates();
      const Math::SpatialVector<Real> surfaceNormal{-physical.x(), -2 * physical.y(), 1};
      const Math::SpatialVector<Real> tangent{
        chord.x(), chord.y(), physical.x() * chord.x() + 2 * physical.y() * chord.y()};
      const auto value = normal.getValue(point);
      EXPECT_NEAR(value.dot(surfaceNormal), 0, 1e-11);
      EXPECT_NEAR(value.dot(tangent), 0, 1e-11);
      checkOrientation(normal, point, *cell);
    }
  }
}
