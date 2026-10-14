/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief P0 value types and cellwise degrees of freedom. */

#include <gtest/gtest.h>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  class P0SpaceTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  /**
   * @brief Certifies the cellwise zero gradient and its explicit trace policy.
   *
   * Different coefficients on adjacent cells distinguish the broken gradient
   * from the distributional derivative. The zero vector has ambient dimension.
   * Direct and indexed samples must agree exactly for real and complex fields.
   * Boundary traces are automatic; interior traces require a selected cell.
   */
  TEST_P(P0SpaceTest, ScalarGradientsAndTracesAreZero)
  {
    const auto geometry = GetParam();
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    LocalMesh mesh;
    if (dimension == 1)
      mesh = LocalMesh::UniformGrid(geometry, {3});
    else if (dimension == 2)
      mesh = LocalMesh::UniformGrid(geometry, {3, 3});
    else
      mesh = LocalMesh::UniformGrid(geometry, {3, 3, 3});
    mesh.getConnectivity().compute(dimension - 1, dimension);
    for (auto cell = mesh.getCell(); cell; ++cell)
      mesh.setAttribute(cell.key(), Geometry::Attribute(cell->getIndex() + 1));
    const auto check = [&]<class Scalar>() {
      P0<Scalar> space(mesh);
      GridFunction field(space);
      for (size_t a = 0; a < field.getSize(); ++a)
      {
        if constexpr (std::is_same_v<Scalar, Complex>)
          field[a] = Scalar(a + 1, 2 * a + 1);
        else
          field[a] = Scalar(a + 1);
      }
      const auto checkValue = [&](const auto& gradient, const IntegrationPoint& point) {
        const auto direct = gradient.getValue(point.getPoint());
        const auto indexed = gradient.getValue(point);
        ASSERT_EQ(direct.size(), mesh.getSpaceDimension());
        ASSERT_EQ(indexed.size(), mesh.getSpaceDimension());
        for (size_t axis = 0; axis < mesh.getSpaceDimension(); ++axis)
        {
          EXPECT_EQ(direct(axis), Scalar(0));
          EXPECT_EQ(indexed(axis), Scalar(0));
        }
      };
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto& formula = QF::PolytopeQuadratureFormula::get(2, cell->getGeometry());
        const PolytopeQuadrature quadrature(*cell, formula);
        for (size_t q = 0; q < quadrature.getSize(); ++q)
          checkValue(Grad(field), IntegrationPoint(quadrature, q));
      }
      bool boundary = false, interior = false;
      for (auto face = mesh.getFace(); face; ++face)
      {
        if (face->isBoundary() ? boundary : interior)
          continue;
        const auto& formula = QF::PolytopeQuadratureFormula::get(2, face->getGeometry());
        const PolytopeQuadrature quadrature(*face, formula);
        const IntegrationPoint point(quadrature, 0);
        if (face->isBoundary())
        {
          checkValue(Grad(field), point);
          boundary = true;
        }
        else
        {
          EXPECT_THROW(Grad(field).getValue(point), Alert::Exception);
          EXPECT_THROW(Grad(field).traceOf(0).getValue(point), Alert::Exception);
          const auto& incident = mesh.getConnectivity().getIncidence(
            {dimension - 1, dimension}, face->getIndex());
          for (const auto cell : incident)
          {
            const auto attribute = mesh.getPolytope(dimension, cell)->getAttribute();
            ASSERT_TRUE(attribute);
            checkValue(Grad(field).traceOf(*attribute), point);
          }
          interior = true;
        }
        if (boundary && interior)
          break;
      }
      EXPECT_TRUE(boundary);
      EXPECT_TRUE(interior);
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
  }

  TEST_P(P0SpaceTest, ScalarRealAndComplexAreCellwiseConstant)
  {
    const auto geometry = GetParam();
    const size_t dim = Polytope::Traits(geometry).getDimension();
    LocalMesh mesh;
    if (dim == 1)
      mesh = LocalMesh::UniformGrid(geometry, {3});
    if (dim == 2)
      mesh = LocalMesh::UniformGrid(geometry, {3, 3});
    if (dim == 3)
      mesh = LocalMesh::UniformGrid(geometry, {3, 3, 3});

    P0 realSpace(mesh);
    ComplexP0<LocalMesh> complexSpace(mesh);
    ASSERT_EQ(realSpace.getSize(), mesh.getCellCount());
    ASSERT_EQ(complexSpace.getSize(), mesh.getCellCount());

    GridFunction realField(realSpace);
    GridFunction complexField(complexSpace);
    for (Index i = 0; i < mesh.getCellCount(); ++i)
    {
      realField[i] = Real(i + 1);
      complexField[i] = Complex(Real(i + 1), -Real(i + 2));
      EXPECT_EQ(realSpace.getDOFs(dim, i)[0], i);
      EXPECT_EQ(complexSpace.getDOFs(dim, i)[0], i);
      const auto cell = mesh.getCell(i);
      const Point point(*cell, Polytope::Traits(geometry).getCentroid());
      EXPECT_NEAR(realField(point), Real(i + 1), 1e-12);
      EXPECT_NEAR(
        std::abs(complexField(point) - Complex(Real(i + 1), -Real(i + 2))), 0, 1e-12);
    }
  }

  TEST_P(P0SpaceTest, RealAndComplexVectorsHaveIndependentCellDOFs)
  {
    const auto geometry = GetParam();
    const size_t dim = Polytope::Traits(geometry).getDimension();
    LocalMesh mesh;
    if (dim == 1)
      mesh = LocalMesh::UniformGrid(geometry, {3});
    if (dim == 2)
      mesh = LocalMesh::UniformGrid(geometry, {3, 3});
    if (dim == 3)
      mesh = LocalMesh::UniformGrid(geometry, {3, 3, 3});

    P0<Math::SpatialVector<Real>, LocalMesh> realSpace(mesh, 2);
    P0<Math::SpatialVector<Complex>, LocalMesh> complexSpace(mesh, 2);
    ASSERT_EQ(realSpace.getSize(), 2 * mesh.getCellCount());
    ASSERT_EQ(complexSpace.getSize(), 2 * mesh.getCellCount());
    GridFunction realField(realSpace);
    GridFunction complexField(complexSpace);

    for (Index i = 0; i < mesh.getCellCount(); ++i)
    {
      for (Index component = 0; component < 2; ++component)
      {
        const Index global = 2 * i + component;
        EXPECT_EQ(realSpace.getGlobalIndex({dim, i}, component), global);
        EXPECT_EQ(complexSpace.getGlobalIndex({dim, i}, component), global);
        realField[global] = Real(global + 1);
        complexField[global] = Complex(Real(global + 1), Real(component + 1));
      }
      const auto cell = mesh.getCell(i);
      const Point point(*cell, Polytope::Traits(geometry).getCentroid());
      const auto realValue = realField(point);
      const auto complexValue = complexField(point);
      for (Index component = 0; component < 2; ++component)
      {
        const Index global = 2 * i + component;
        EXPECT_NEAR(realValue(component), Real(global + 1), 1e-12);
        EXPECT_NEAR(std::abs(complexValue(component) -
                      Complex(Real(global + 1), Real(component + 1))),
          0, 1e-12);
      }
    }
  }

  /**
   * @brief Zero gradients retain ambient dimension on embedded entities.
   *
   * A point, segment and surface embedded in three dimensions all produce
   * three zero components, rather than a vector of topological dimension.
   */
  TEST(P0Gradient, EmbeddedAmbientDimension)
  {
    const auto check = []<class Scalar>(LocalMesh& mesh) {
      P0<Scalar> space(mesh);
      GridFunction field(space);
      for (size_t a = 0; a < field.getSize(); ++a)
        field[a] = Scalar(a + 1);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto& formula = QF::PolytopeQuadratureFormula::get(2, cell->getGeometry());
        const PolytopeQuadrature quadrature(*cell, formula);
        for (size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const IntegrationPoint point(quadrature, q);
          const auto direct = Grad(field).getValue(point.getPoint());
          const auto indexed = Grad(field).getValue(point);
          ASSERT_EQ(direct.size(), 3);
          ASSERT_EQ(indexed.size(), 3);
          for (size_t axis = 0; axis < 3; ++axis)
          {
            EXPECT_EQ(direct(axis), Scalar(0));
            EXPECT_EQ(indexed(axis), Scalar(0));
          }
        }
      }
    };
    auto point = LocalMesh::Builder().initialize(3).nodes(1).vertex({1, 2, 3}).finalize();
    auto segment = LocalMesh::Builder()
                     .initialize(3)
                     .nodes(2)
                     .vertex({0, 0, 0})
                     .vertex({1, 2, 2})
                     .polytope(Polytope::Type::Segment, {{0, 1}})
                     .finalize();
    auto surface = LocalMesh::Builder()
                     .initialize(3)
                     .nodes(3)
                     .vertex({0, 0, 0})
                     .vertex({1, 0, 1})
                     .vertex({0, 1, 1})
                     .polytope(Polytope::Type::Triangle, {{0, 1, 2}})
                     .finalize();
    for (auto* mesh : {&point, &segment, &surface})
    {
      check.template operator()<Real>(*mesh);
      check.template operator()<Complex>(*mesh);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, P0SpaceTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge));
}
