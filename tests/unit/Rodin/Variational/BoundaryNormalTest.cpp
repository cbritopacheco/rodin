/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file BoundaryNormalTest.cpp
 * @brief Tests for the BoundaryNormal class.
 */
#include <gtest/gtest.h>
#include <cmath>
#include <string>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief Verifies construction for variational boundary normal by checking exact expected values.
TEST(Rodin_Variational_BoundaryNormal, Construction)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  BoundaryNormal n(mesh);

  EXPECT_EQ(n.getDimension(), mesh.getSpaceDimension());
}

/// @brief Verifies copy for variational boundary normal by checking exact expected values, copy semantics.
TEST(Rodin_Variational_BoundaryNormal, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  BoundaryNormal n(mesh);
  auto copy = n;

  EXPECT_EQ(copy.getDimension(), n.getDimension());
}

/// @brief Verifies get order for variational boundary normal by checking exact expected values, true predicates.
TEST(Rodin_Variational_BoundaryNormal, GetOrder)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  BoundaryNormal n(mesh);
  auto order = n.getOrder(*it);
  EXPECT_TRUE(order.has_value());
  EXPECT_EQ(*order, 0u);
}

/// @brief Verifies unit normal on boundary for variational boundary normal by checking tolerance-based numerical results, true predicates.
TEST(Rodin_Variational_BoundaryNormal, UnitNormalOnBoundary)
{
  // Create a simple mesh and check normal at a boundary face
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  mesh.getConnectivity().compute(1, 2);

  BoundaryNormal n(mesh);

  // Find a boundary edge (codimension-1 entity)
  const size_t dimension = mesh.getDimension();
  bool foundBoundary = false;
  for (auto it = mesh.getPolytope(dimension - 1); !it.end(); ++it)
  {
    const auto& polytope = *it;
    if (mesh.isBoundary(polytope.getIndex()))
    {
      // Create a point on this boundary face
      const Math::Vector<Real> rc{{ 0.5 }};
      Point p(polytope, rc);

      auto val = n.getValue(p);

      // Check that it's a unit vector (norm ≈ 1)
      Real norm = 0.0;
      for (size_t i = 0; i < static_cast<size_t>(val.size()); ++i)
        norm += val(i) * val(i);
      norm = std::sqrt(norm);

      EXPECT_NEAR(norm, 1.0, 1e-10);
      foundBoundary = true;
      break;
    }
  }
  EXPECT_TRUE(foundBoundary);
}

/// @brief On the unit cube, outward normals agree with the coordinate planes.
TEST(Rodin_Variational_BoundaryNormal, AllCubeFacesHaveExactOutwardNormals)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Hexahedron, {2, 2, 2});
  mesh.getConnectivity().compute(2, 3);
  BoundaryNormal normal(mesh);
  size_t count = 0;
  for (auto face = mesh.getBoundary(); face; ++face)
  {
    Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
    const auto actual = normal.getValue(point);
    const auto coordinates = point.getCoordinates();
    Math::SpatialVector<Real> expected(3);
    expected.setZero();
    for (size_t axis = 0; axis < 3; ++axis)
    {
      if (std::abs(coordinates[axis]) < 1e-12)
        expected[axis] = -1;
      else if (std::abs(coordinates[axis] - 1) < 1e-12)
        expected[axis] = 1;
    }
    EXPECT_LT((actual - expected).norm(), 1e-12);
    ++count;
  }
  EXPECT_EQ(count, 6);
}

/// @brief Corner normals require a selected incident boundary trace.
TEST(Rodin_Variational_BoundaryNormal, CornerTraceSelectsIncidentEdge)
{
  auto mesh = LocalMesh::Builder()
                .initialize(2)
                .nodes(3)
                .vertex({0, 0})
                .vertex({1, 0})
                .vertex({0, 1})
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .finalize();
  mesh.getConnectivity().compute(0, 1);
  mesh.getConnectivity().compute(1, 2);
  for (auto edge = mesh.getPolytope(1); edge; ++edge)
    mesh.setAttribute({1, edge->getIndex()}, edge->getIndex() + 1);
  auto vertex = mesh.getPolytope(0, 0);
  Point corner(*vertex, Math::SpatialPoint(0));
  const auto& incident = mesh.getConnectivity().getIncidence({0, 1}, 0);
  size_t checked = 0;
  for (auto index : incident)
  {
    if (!mesh.isBoundary(index))
      continue;
    auto edge = mesh.getPolytope(1, index);
    Point midpoint(*edge, Polytope::Traits(edge->getGeometry()).getCentroid());
    BoundaryNormal normal(mesh);
    const auto expected = normal.getValue(midpoint);
    normal.traceOf(index + 1);
    EXPECT_LT((normal.getValue(corner) - expected).norm(), 1e-12);
    ++checked;
  }
  EXPECT_EQ(checked, 2);
  BoundaryNormal invalid(mesh);
  invalid.traceOf(999);
  EXPECT_THROW(invalid.getValue(corner), Alert::Exception);
}

/// @brief Zero-measure boundaries report native diagnostics and recover after geometry repair.
TEST(Rodin_Variational_BoundaryNormal, CollapsedBoundaryRejectsAndRecovers)
{
  for (const size_t dimension : {2u, 3u})
  {
    SCOPED_TRACE(dimension);
    auto builder = LocalMesh::Builder();
    builder.initialize(dimension).nodes(dimension + 1);
    for (size_t vertex = 0; vertex <= dimension; ++vertex)
    {
      Math::SpatialPoint coordinate(dimension);
      coordinate.setZero();
      if (vertex > 1)
        coordinate[vertex - 1] = 1;
      builder.vertex(coordinate);
    }
    if (dimension == 2)
      builder.polytope(Polytope::Type::Triangle, {0, 1, 2});
    else
      builder.polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3});
    auto mesh = builder.finalize();
    mesh.getConnectivity().compute(dimension - 1, dimension);
    Optional<Index> selected;
    for (auto face = mesh.getBoundary(); face; ++face)
    {
      bool containsLastVertex = false;
      for (const Index vertex : face->getVertices())
        containsLastVertex |= vertex == dimension;
      if (!containsLastVertex)
        selected = face->getIndex();
    }
    ASSERT_TRUE(selected);
    const auto face = mesh.getPolytope(dimension - 1, *selected);
    const auto reference = Polytope::Traits(face->getGeometry()).getCentroid();
    const Point collapsed(*face, reference);
    BoundaryNormal normal(mesh);
    Math::SpatialVector<Real> spatial;
    Math::Vector<Real> dynamic;
    const auto reject = [&](auto& output) {
      try
      {
        normal.interpolate(output, collapsed);
        FAIL() << "A collapsed boundary has no unit normal.";
      }
      catch (const Alert::Exception& error)
      {
        const std::string message = error.what();
        EXPECT_NE(message.find("In member function"), std::string::npos);
        EXPECT_NE(message.find("Cannot normalize BoundaryNormal"), std::string::npos);
      }
    };
    reject(spatial);
    reject(dynamic);

    mesh.setVertexCoordinates(1, 1, 0);
    const Point repaired(*face, reference);
    Math::SpatialVector<Real> expected(dimension);
    expected.setZero();
    expected[dimension - 1] = -1;
    normal.interpolate(spatial, repaired);
    normal.interpolate(dynamic, repaired);
    ASSERT_EQ(spatial.size(), dimension);
    ASSERT_EQ(dynamic.size(), dimension);
    for (Index component = 0; component < dimension; ++component)
    {
      EXPECT_NEAR(spatial[component], expected[component], 1e-12);
      EXPECT_NEAR(dynamic[component], expected[component], 1e-12);
    }
  }
}

/// @brief Unsupported endpoint Jacobians are rejected through both output-storage overloads.
TEST(Rodin_Variational_BoundaryNormal, OneDimensionalBoundaryReportsUnsupportedJacobian)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  mesh.getConnectivity().compute(0, 1);
  BoundaryNormal normal(mesh);
  size_t checked = 0;
  for (auto face = mesh.getBoundary(); face; ++face)
  {
    const Point point(*face, Math::SpatialPoint(0));
    Math::SpatialVector<Real> spatial;
    Math::Vector<Real> dynamic;
    const auto reject = [&](auto& output) {
      try
      {
        normal.interpolate(output, point);
        FAIL() << "A one-row boundary Jacobian is unsupported.";
      }
      catch (const Alert::Exception& error)
      {
        const std::string message = error.what();
        EXPECT_NE(message.find("In member function"), std::string::npos);
        EXPECT_NE(message.find("boundary Jacobian with 1 rows"), std::string::npos);
      }
    };
    reject(spatial);
    reject(dynamic);
    ++checked;
  }
  EXPECT_EQ(checked, 2);
}

/// @brief Parent and retained-submesh points give exact normals on their own boundaries.
TEST(Rodin_Variational_BoundaryNormal, ParentAndSubmeshPointTransport)
{
  for (const auto geometry :
    {Polytope::Type::Triangle, Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto parent = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 3));
    parent.getConnectivity().compute(dimension - 1, dimension);
    for (auto cell = parent.getCell(); cell; ++cell)
    {
      const Point center(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      parent.setAttribute({dimension, cell->getIndex()}, center(0) < 1 ? 17 : 23);
    }
    auto retained = parent.keep(17);
    retained.getConnectivity().compute(dimension - 1, dimension);
    BoundaryNormal parentNormal(parent), retainedNormal(retained);
    auto copied = retainedNormal;
    auto moved = std::move(copied);
    std::unique_ptr<BoundaryNormal> cloned(retainedNormal.copy());
    size_t restrictedCount = 0, includedCount = 0;
    for (auto face = retained.getBoundary(); face; ++face)
    {
      const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
      Math::SpatialVector<Real> expected(dimension);
      expected.setZero();
      for (size_t axis = 0; axis < dimension; ++axis)
      {
        const Real upper = axis == 0 ? 1 : 2;
        if (std::abs(point(axis)) < 1e-12)
          expected[axis] = -1;
        else if (std::abs(point(axis) - upper) < 1e-12)
          expected[axis] = 1;
      }
      ASSERT_NEAR(expected.norm(), 1, 1e-12);
      EXPECT_LT((retainedNormal.getValue(point) - expected).norm(), 1e-12);
      const auto included = parent.inclusion(point);
      ASSERT_TRUE(included);
      EXPECT_LT((retainedNormal.getValue(*included) - expected).norm(), 1e-12);
      EXPECT_LT((moved.getValue(*included) - expected).norm(), 1e-12);
      EXPECT_LT((cloned->getValue(*included) - expected).norm(), 1e-12);
      ++restrictedCount;
      if (parent.isBoundary(included->getPolytope().getIndex()))
      {
        EXPECT_LT((parentNormal.getValue(point) - expected).norm(), 1e-12);
        ++includedCount;
      }
    }
    EXPECT_GT(restrictedCount, 0);
    EXPECT_GT(includedCount, 0);
    bool checkedExcluded = false;
    for (auto face = parent.getBoundary(); face; ++face)
    {
      const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
      if (std::abs(point(0) - 2) < 1e-12)
      {
        EXPECT_THROW(retainedNormal.getValue(point), Alert::Exception);
        checkedExcluded = true;
        break;
      }
    }
    EXPECT_TRUE(checkedExcluded);
    auto unrelated =
      LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    unrelated.getConnectivity().compute(dimension - 1, dimension);
    auto face = unrelated.getBoundary();
    ASSERT_TRUE(face);
    const Point unrelatedPoint(
      *face, Polytope::Traits(face->getGeometry()).getCentroid());
    EXPECT_THROW(parentNormal.getValue(unrelatedPoint), Alert::Exception);
    EXPECT_THROW(retainedNormal.getValue(unrelatedPoint), Alert::Exception);
  }
}

/// @brief Missing corner/edge traces raise member diagnostics and permit later trace selection.
TEST(Rodin_Variational_BoundaryNormal, MissingTraceRaisesAndPermitsRecovery)
{
  for (const size_t dimension : {size_t(2), size_t(3)})
  {
    SCOPED_TRACE(dimension);
    auto builder = LocalMesh::Builder();
    builder.initialize(dimension).nodes(dimension + 1);
    for (size_t vertex = 0; vertex <= dimension; ++vertex)
    {
      Math::SpatialPoint coordinate(dimension);
      coordinate.setZero();
      if (vertex != 0)
        coordinate(vertex - 1) = 1;
      builder.vertex(coordinate);
    }
    if (dimension == 2)
      builder.polytope(Polytope::Type::Triangle, {0, 1, 2});
    else
      builder.polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3});
    auto mesh = builder.finalize();
    mesh.getConnectivity().compute(dimension - 2, dimension - 1);
    mesh.getConnectivity().compute(dimension - 1, dimension);
    const auto entity = mesh.getPolytope(dimension - 2, 0);
    const Point point(*entity, Polytope::Traits(entity->getGeometry()).getCentroid());
    const auto& incident =
      mesh.getConnectivity().getIncidence({dimension - 2, dimension - 1}, 0);
    ASSERT_EQ(incident.size(), 2u);
    BoundaryNormal normal(mesh);
    Math::SpatialVector<Real> spatial(dimension);
    Math::Vector<Real> dynamic(dimension);
    const auto requireDiagnostic = [&](auto& output) {
      try
      {
        normal.interpolate(output, point);
        FAIL() << "A normal with two incident traces needs an explicit selection";
      }
      catch (const Alert::Exception& error)
      {
        const std::string message = error.what();
        EXPECT_NE(message.find("member function"), std::string::npos);
        EXPECT_NE(message.find("interpolate"), std::string::npos);
        EXPECT_NE(message.find("No trace domain provided"), std::string::npos);
      }
    };
    requireDiagnostic(spatial);
    requireDiagnostic(dynamic);
    EXPECT_EQ(normal.getTraceDomain().size(), 0u);

    const Index faceIndex = *incident.begin();
    mesh.setAttribute({dimension - 1, faceIndex}, 17);
    normal.traceOf(17);
    const auto face = mesh.getPolytope(dimension - 1, faceIndex);
    size_t missingVertex = 0;
    for (; missingVertex <= dimension; ++missingVertex)
    {
      bool found = false;
      for (const Index vertex : face->getVertices())
        found |= vertex == missingVertex;
      if (!found)
        break;
    }
    ASSERT_LE(missingVertex, dimension);
    Math::SpatialVector<Real> expected(dimension);
    expected.setZero();
    if (missingVertex == 0)
      expected.setConstant(Real(1) / std::sqrt(Real(dimension)));
    else
      expected[missingVertex - 1] = -1;
    normal.interpolate(spatial, point);
    normal.interpolate(dynamic, point);
    for (size_t axis = 0; axis < dimension; ++axis)
    {
      EXPECT_NEAR(spatial[axis], expected[axis], 1e-12);
      EXPECT_NEAR(dynamic[axis], expected[axis], 1e-12);
    }
  }
}
