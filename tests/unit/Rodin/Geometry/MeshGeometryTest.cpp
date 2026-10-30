/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief Copies retain independently owned curved charts after source destruction.
TEST(MeshGeometry, CopyPreservesCurvature)
{
  const auto mesh = [] {
    auto source = LocalMesh::Builder()
                    .initialize(2)
                    .nodes(3)
                    .vertex({0, 0})
                    .vertex({1, 0})
                    .vertex({0, 1})
                    .polytope(Polytope::Type::Triangle, {0, 1, 2})
                    .finalize();
    RealH1Element<2> element(Polytope::Type::Triangle);
    PointCloud nodes(2, element.getCount());
    for (size_t local = 0; local < element.getCount(); ++local)
    {
      const auto point = element.getNode(local);
      nodes(0, local) = point[0];
      nodes(1, local) = point[1] + 0.4 * point[0] * (1 - point[0] - point[1]);
    }
    source.setPolytopeTransformation({2, 0},
      new ParametricTransformation<RealH1Element<2>>(
        std::move(nodes), std::move(element)));
    const auto cell = source.getCell();
    const auto* original = &cell->getTransformation();
    LocalMesh copied(source);
    EXPECT_NE(&copied.getCell()->getTransformation(), original);
    return copied;
  }();
  const auto cell = mesh.getCell();
  EXPECT_EQ(cell->getTransformation().getDegree(), 2u);
  const Point point(*cell, Math::SpatialPoint{0.5, 0});
  EXPECT_NEAR(point.x(), 0.5, 1e-14);
  EXPECT_NEAR(point.y(), 0.1, 1e-14);
  EXPECT_NEAR(cell->getMeasure(), Real(13) / 30, 1e-13);
}

/// @brief Field-dependent displacements see the same original geometry at every vertex.
TEST(MeshGeometry, DisplacementIgnoresCacheWarmth)
{
  for (const bool warm : {false, true})
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
    mesh.getConnectivity().compute(0, 1);
    for (auto cell = mesh.getCell(); cell; ++cell)
      mesh.setAttribute({1, cell->getIndex()}, 1);
    P1<Real> space(mesh);
    GridFunction field(space);
    field = RealFunction([](const Point& point) { return point.x(); });
    mesh.flush();
    if (warm)
    {
      for (auto cell = mesh.getCell(); cell; ++cell)
        cell->getTransformation();
    }
    const auto original = mesh.getVertices();
    auto gradient = Grad(field);
    gradient.traceOf(1);
    mesh.displace(Real(0.1) * gradient);
    for (Index index = 0; index < mesh.getVertexCount(); ++index)
      EXPECT_NEAR(mesh.getVertexCoordinates(index)[0], original[index][0] + 0.1, 1e-14);
  }
}

/// @brief Failed evaluation or a wrong ambient extent leaves all coordinates unchanged.
TEST(MeshGeometry, FailedDisplacementPreservesCoordinates)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
  const auto original = mesh.getVertices();
  auto displacement = VectorFunction(size_t(1), [](const Point& point) {
    if (point.x() > 0)
      Alert::Exception() << "Deliberate displacement failure." << Alert::Raise;
    return Math::SpatialVector<Real>{0.1};
  });
  EXPECT_THROW(mesh.displace(displacement), Alert::Exception);
  EXPECT_THROW(mesh.displace(VectorFunction{0.1, 0.2}), Alert::Exception);
  for (Index index = 0; index < mesh.getVertexCount(); ++index)
    EXPECT_EQ(mesh.getVertexCoordinates(index)[0], original[index][0]);
}

/// @brief Curved maps invert their actual polynomial image across geometries and scales.
TEST(MeshGeometry, CurvedInverseRoundtrip)
{
  for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Wedge, Polytope::Type::Hexahedron})
  {
    const Polytope::Traits traits(geometry);
    const size_t dimension = traits.getDimension();
    for (const Real scale : {Real(1e-18), Real(1), Real(1e18)})
    {
      for (const bool embedded : {false, true})
      {
        if (embedded && dimension == 3)
          continue;
        RealH1Element<2> element(geometry);
        PointCloud nodes(dimension + embedded, element.getCount());
        for (size_t local = 0; local < element.getCount(); ++local)
        {
          const auto r = element.getNode(local);
          Real sum = 0;
          for (size_t i = 0; i < dimension; ++i)
          {
            nodes(i, local) = scale * (r[i] + 0.15 * r[i] * r[i]);
            sum += r[i] * r[i];
          }
          if (embedded)
            nodes(dimension, local) = scale * 0.2 * sum;
        }
        ParametricTransformation<RealH1Element<2>> map(
          std::move(nodes), std::move(element));
        for (size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
        {
          const Math::SpatialPoint expected =
            (traits.getCentroid() + traits.getVertex(vertex)) / 2;
          Math::SpatialPoint physical, actual;
          map.transform(physical, expected);
          map.inverse(actual, physical);
          EXPECT_LT((actual - expected).norm(), 3e-12)
            << "geometry=" << static_cast<int>(geometry) << " scale=" << scale;
        }
      }
    }
  }
}

/// @brief Even degree-one tensor-product geometry may require nonlinear inversion.
TEST(MeshGeometry, BilinearInverseAndFailure)
{
  RealP1Element element(Polytope::Type::Quadrilateral);
  PointCloud nodes(2, element.getCount());
  for (size_t local = 0; local < element.getCount(); ++local)
  {
    const auto point = element.getNode(local);
    nodes(0, local) = point[0];
    nodes(1, local) = point[1] * (1 + 0.3 * point[0]);
  }
  ParametricTransformation<RealP1Element> map(std::move(nodes), std::move(element));
  const Math::SpatialPoint expected{0.3, 0.7};
  Math::SpatialPoint physical, actual;
  map.transform(physical, expected);
  map.inverse(actual, physical);
  EXPECT_LT((actual - expected).norm(), 1e-13);
  EXPECT_THROW(map.inverse(actual, Math::SpatialPoint{1}), Alert::Exception);
  RealH1Element<2> singularElement(Polytope::Type::Segment);
  PointCloud constant(1, singularElement.getCount());
  constant.setZero();
  ParametricTransformation<RealH1Element<2>> singular(
    std::move(constant), std::move(singularElement));
  EXPECT_THROW(singular.inverse(actual, Math::SpatialPoint{1}), Alert::Exception);
}

/// @brief Cubic tetrahedral volume uses determinant degree rather than map degree.
TEST(MeshGeometry, CubicMeasure)
{
  auto mesh = LocalMesh::Builder()
                .initialize(3)
                .nodes(4)
                .vertex({0, 0, 0})
                .vertex({1, 0, 0})
                .vertex({0, 1, 0})
                .vertex({0, 0, 1})
                .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
                .finalize();
  RealH1Element<3> element(Polytope::Type::Tetrahedron);
  PointCloud nodes(3, element.getCount());
  for (size_t local = 0; local < element.getCount(); ++local)
  {
    const auto r = element.getNode(local);
    for (size_t i = 0; i < 3; ++i)
      nodes(i, local) = r[i] + r[i] * r[i] * r[i];
  }
  mesh.setPolytopeTransformation({3, 0},
    new ParametricTransformation<RealH1Element<3>>(std::move(nodes), std::move(element)));
  const Real exact = Real(569) / 1680;
  const auto cell = mesh.getCell();
  EXPECT_NEAR(cell->getMeasure(), exact, 2e-13);
  EXPECT_NEAR(cell->getMeasure(12), exact, 2e-13);
  EXPECT_NEAR(mesh.getMeasure(3), exact, 2e-13);
  EXPECT_GT(std::abs(cell->getMeasure(3) - exact), 1e-4);
}
