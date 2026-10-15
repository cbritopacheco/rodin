/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Attribute ownership, selected interfaces and embedded metric contracts. */
#include <gtest/gtest.h>
#include <array>
#include <limits>
#include <memory>
#include <vector>
#include "Rodin/Geometry.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

using namespace Rodin;
using namespace Rodin::Geometry;

/// @brief Repeated attribute initialization preserves values and sparse chart insertion owns clones.
TEST(Geometry_MetadataLifecycle, IdempotentAttributesAndSparseChartClones)
{
  AttributeIndex attributes;
  attributes.initialize(3);
  PolytopeTransformationIndex transformations;
  transformations.initialize(3);
  constexpr Index Slot = 5;
  for (size_t dimension = 0; dimension <= 3; ++dimension)
  {
    attributes.set({dimension, Slot}, Attribute(17 + dimension));
    {
      IdentityTransformation original(dimension);
      IdentityTransformation copied(original);
      IdentityTransformation moved(std::move(copied));
      transformations.set(
        {dimension, Slot}, std::unique_ptr<PolytopeTransformation>(moved.copy()));
    }
  }
  attributes.initialize(3);
  for (size_t dimension = 0; dimension <= 3; ++dimension)
  {
    SCOPED_TRACE(dimension);
    EXPECT_EQ(
      attributes.get({dimension, Slot}, Slot + 1), Optional<Attribute>(17 + dimension));
    EXPECT_FALSE(attributes.get({dimension, 0}, Slot + 1));
    EXPECT_EQ(transformations.find(dimension, 0), nullptr);
    const auto* chart = transformations.find(dimension, Slot);
    ASSERT_NE(chart, nullptr);
    EXPECT_EQ(chart->getOrder(), 1u);
    Math::SpatialPoint reference(dimension), physical;
    for (size_t i = 0; i < dimension; ++i)
      reference[i] = Real(i + 1) / 2;
    chart->transform(physical, reference);
    EXPECT_EQ(physical.size(), dimension);
    for (size_t i = 0; i < dimension; ++i)
      EXPECT_EQ(physical[i], reference[i]);
    Math::SpatialMatrix<Real> jacobian;
    chart->jacobian(jacobian, reference);
    EXPECT_EQ(jacobian.rows(), dimension);
    EXPECT_EQ(jacobian.cols(), dimension);
    for (size_t i = 0; i < dimension; ++i)
    {
      for (size_t j = 0; j < dimension; ++j)
        EXPECT_EQ(jacobian(i, j), Real(i == j));
    }
  }
}

/// @brief Connectivity lookup respects unordered keys while retaining the first orientation.
TEST(Geometry_MetadataLifecycle, ConnectivityKeyLookupAndDuplicateInsertion)
{
  SequentialConnectivity connectivity;
  connectivity.initialize(2).nodes(4);
  const Polytope::Key first{0, 1, 2}, reversed{2, 1, 0}, second{1, 2, 3};
  connectivity.polytope(Polytope::Type::Triangle, first);
  connectivity.polytope(Polytope::Type::Triangle, reversed);
  connectivity.polytope(Polytope::Type::Triangle, second);
  ASSERT_EQ(connectivity.getCount(2), 2u);
  EXPECT_EQ(connectivity.getIndex(2, first), Optional<Index>(0));
  EXPECT_EQ(connectivity.getIndex(2, reversed), Optional<Index>(0));
  EXPECT_EQ(connectivity.getIndex(2, second), Optional<Index>(1));
  EXPECT_FALSE(connectivity.getIndex(2, Polytope::Key{0, 2, 3}));
  const auto& index = connectivity.getIndexMap(2);
  ASSERT_EQ(index.left.size(), 2u);
  EXPECT_EQ(index.right.size(), 2u);
  EXPECT_EQ(
    IndexVector(index.left[0]->begin(), index.left[0]->end()), IndexVector({0, 1, 2}));
  EXPECT_EQ(
    IndexVector(index.left[1]->begin(), index.left[1]->end()), IndexVector({1, 2, 3}));
  connectivity.discover(2, 1).restrict(1, 2);
  EXPECT_EQ(connectivity.getCount(1), 5u);
  size_t shared = 0;
  for (Index i = 0; i < connectivity.getCount(1); ++i)
    shared += connectivity.getIncidence({1, 2}, i).size() == 2;
  EXPECT_EQ(shared, 1u);
}

/// @brief Populated attribute assignment replaces sparse extents and owns independent slots.
TEST(Geometry_MetadataLifecycle, SparseAttributeAssignmentReplacesPopulatedStorage)
{
  constexpr size_t bound = 16;
  const std::array<size_t, 4> sizes{3, 7, 0, 12};
  std::array<std::vector<Optional<Attribute>>, 4> expected;
  AttributeIndex destination;
  destination.initialize(3);
  for (size_t d = 0; d < 4; ++d)
  {
    destination.resize(d, d % 2 ? 2 : 9);
    destination.set({d, 0}, Attribute(999));
    expected[d].resize(bound);
  }
  {
    AttributeIndex source;
    source.initialize(3);
    for (size_t d = 0; d < 4; ++d)
    {
      source.resize(d, sizes[d]);
      for (Index i = 0; i < sizes[d]; ++i)
      {
        if (i % 3 == 0)
          expected[d][i] = Attribute(0);
        else if (i % 3 == 1)
          expected[d][i] = std::numeric_limits<Attribute>::max();
        source.set({d, i}, expected[d][i]);
      }
    }
    EXPECT_EQ(&(destination = source), &destination);
    for (size_t d = 0; d < 4; ++d)
    {
      source.set({d, bound - 1}, Attribute(73));
      source.set({d, 0}, Attribute(71));
    }
  }
  for (size_t d = 0; d < 4; ++d)
  {
    FlatSet<Attribute> attributes;
    for (Index i = 0; i < bound; ++i)
    {
      EXPECT_EQ(destination.get({d, i}, bound), expected[d][i]);
      if (expected[d][i])
        attributes.insert(*expected[d][i]);
    }
    EXPECT_EQ(destination.getAttributes(d), attributes);
  }
  AttributeIndex copied(destination);
  destination.set({1, 1}, Attribute(83));
  EXPECT_EQ(copied.get({1, 1}, bound), expected[1][1]);
  const AttributeIndex& alias = copied;
  copied = alias;
  EXPECT_EQ(copied.get({1, 1}, bound), expected[1][1]);
}

/// @brief Interface relabeling honors selection, unordered material pairs and absent attributes.
TEST(Geometry_MetadataLifecycle, SelectedInterfaceTracePreservesOtherFaces)
{
  LocalMesh::Builder builder;
  builder.initialize(2).nodes(10);
  for (size_t y = 0; y < 2; ++y)
  {
    for (size_t x = 0; x < 5; ++x)
      builder.vertex({Real(x), Real(y)});
  }
  for (Index x = 0; x < 4; ++x)
    builder.polytope(Polytope::Type::Quadrilateral, {x, x + 1, x + 6, x + 5});
  auto mesh = builder.finalize();
  mesh.getConnectivity().compute(1, 2);
  mesh.setAttribute({2, 0}, 17);
  mesh.setAttribute({2, 1}, 23);
  mesh.setAttribute({2, 2}, 17);
  mesh.setAttribute({2, 3}, std::nullopt);
  auto interfaceIndex = [&](const Face& face) -> size_t {
    const Point center(face, Polytope::Traits(face.getGeometry()).getCentroid());
    if (center.y() == Real(0.5) && center.x() > 0 && center.x() < 4)
      return static_cast<size_t>(center.x());
    return 0;
  };
  for (auto face = mesh.getFace(); face; ++face)
  {
    const size_t index = interfaceIndex(*face);
    mesh.setAttribute(face.key(),
      index == 0     ? Optional<Attribute>(41)
        : index == 1 ? Optional<Attribute>(31)
        : index == 2 ? Optional<Attribute>(37)
                     : std::nullopt);
  }
  auto check = [&](const std::array<Optional<Attribute>, 4>& labels) {
    size_t exterior = 0, interior = 0;
    for (auto face = mesh.getFace(); face; ++face)
    {
      const size_t index = interfaceIndex(*face);
      EXPECT_EQ(face->getAttribute(), labels[index]);
      index == 0 ? ++exterior : ++interior;
    }
    EXPECT_EQ(exterior, 10);
    EXPECT_EQ(interior, 3);
  };
  const Map<Pair<Attribute, Attribute>, Attribute> forward{{{17, 23}, 71}};
  const Map<Pair<Attribute, Attribute>, Attribute> reverse{{{23, 17}, 73}};
  EXPECT_EQ(&mesh.trace(forward, FlatSet<Attribute>{31}), &mesh);
  check({41, 71, 37, std::nullopt});
  mesh.trace(reverse, FlatSet<Attribute>{71, 37});
  check({41, 73, 73, std::nullopt});
  mesh.trace(forward);
  check({41, 71, 71, std::nullopt});
  for (auto face = mesh.getFace(); face; ++face)
  {
    if (interfaceIndex(*face) == 3)
      mesh.setAttribute(face.key(), 31);
  }
  mesh.trace(forward, FlatSet<Attribute>{31});
  check({41, 71, 71, 31});
  mesh.setAttribute({2, 3}, 23);
  mesh.trace(reverse, FlatSet<Attribute>{31});
  check({41, 71, 71, 73});
  EXPECT_EQ(mesh.getAttribute(2, 0), 17);
  EXPECT_EQ(mesh.getAttribute(2, 1), 23);
  EXPECT_EQ(mesh.getAttribute(2, 2), 17);
  EXPECT_EQ(mesh.getAttribute(2, 3), 23);
}

/// @brief Owned cached Jacobians retain exact embedded metrics after source destruction.
TEST(Geometry_MetadataLifecycle, EmbeddedCachedJacobiansOwnMetricAndReferenceData)
{
  for (size_t dimension : {size_t(1), size_t(2)})
  {
    SCOPED_TRACE(dimension);
    LocalMesh::Builder builder;
    builder.initialize(3).nodes(dimension + 1).vertex({1, 2, 3});
    if (dimension == 1)
      builder.vertex({3, 5, 9}).polytope(Polytope::Type::Segment, {0, 1});
    else
      builder.vertex({3, 2, 4}).vertex({1, 5, 5}).polytope(
        Polytope::Type::Triangle, {0, 1, 2});
    auto mesh = builder.finalize();
    const auto cell = mesh.getCell();
    Math::SpatialMatrix<Real> jacobian(3, dimension), inverse(dimension, 3);
    if (dimension == 1)
    {
      const std::array<Real, 3> entries{2, 3, 6};
      for (size_t i = 0; i < 3; ++i)
      {
        jacobian(i, 0) = entries[i];
        inverse(0, i) = entries[i] / 49;
      }
    }
    else
    {
      const std::array<Real, 6> entries{2, 0, 0, 3, 1, 2};
      const std::array<Real, 6> inverseEntries{26, -6, 9, -4, 15, 8};
      for (size_t row = 0; row < 3; ++row)
      {
        for (size_t col = 0; col < 2; ++col)
          jacobian(row, col) = entries[2 * row + col];
      }
      for (size_t row = 0; row < 2; ++row)
      {
        for (size_t col = 0; col < 3; ++col)
          inverse(row, col) = inverseEntries[3 * row + col] / 61;
      }
    }
    Math::SpatialPoint reference(static_cast<std::uint8_t>(dimension));
    reference.setConstant(Real(0.25));
    const Point explicitPoint(jacobian, *cell, reference);
    std::unique_ptr<Point> retained;
    Math::SpatialPoint sample;
    {
      QF::PolytopeQuadratureFormula formula(3, cell->getGeometry());
      sample = formula.getPoint(0);
      const Point original(jacobian, *cell, formula, 0);
      Point copied(original);
      retained = std::make_unique<Point>(std::move(copied));
    }
    const auto originalJacobian = jacobian;
    jacobian.setZero();
    const Real metric = dimension == 1 ? Real(7) : Math::sqrt(Real(61));
    auto check = [&](const Point& point, const Math::SpatialPoint& rc) {
      EXPECT_NEAR((point.getJacobian() - originalJacobian).norm(), 0, 1e-13);
      EXPECT_NEAR((point.getJacobianInverse() - inverse).norm(), 0, 1e-13);
      EXPECT_NEAR(point.getJacobianDeterminant(), metric, 1e-13);
      EXPECT_NEAR(point.getDistortion(), metric, 1e-13);
      EXPECT_NEAR(
        (point.getCoordinates(PointBase::Coordinates::Reference) - rc).norm(), 0, 1e-13);
      Math::SpatialPoint physical{1, 2, 3};
      physical += originalJacobian * rc;
      EXPECT_NEAR((point.getPhysicalCoordinates() - physical).norm(), 0, 1e-13);
    };
    check(explicitPoint, reference);
    check(*retained, sample);
  }
}
