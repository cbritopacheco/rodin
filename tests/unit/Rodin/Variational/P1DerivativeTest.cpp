/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Scalar P1 partial derivatives, charts and one-sided traces.
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "../../../convergence/Convergence.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    /**
     * @brief Compares every partial with the independently implemented gradient.
     *
     * Nonconstant real/complex coefficients exercise every ambient coordinate.
     * Direct and canonical/legacy integration-point paths must agree exactly,
     * including after coefficient updates at unchanged samples. Boundary and
     * selected interior traces use the adjacent cell's complete derivative.
     * Unselected and unmatched interior traces must raise a Rodin exception.
     *
     * @tparam Scalar Real or complex field coefficients.
     * @param mesh Volume, embedded or zero-dimensional mesh.
     */
    template <class Scalar>
    void checkPartials(LocalMesh& mesh)
    {
      P1<Scalar> space(mesh);
      GridFunction field(space);
      const auto checkPoint = [&](const IntegrationPoint& point) {
        for (size_t axis = 0; axis < mesh.getSpaceDimension(); ++axis)
        {
          auto partial = Derivative(axis, field);
          const auto expected = Grad(field).getValue(point.getPoint());
          EXPECT_EQ(partial.getValue(point.getPoint()), expected(axis));
          EXPECT_EQ(partial.getValue(point), expected(axis));
          const IntegrationPoint legacy(
            point.getPoint(), point.getQuadratureFormula(), point.getIndex());
          EXPECT_EQ(partial.getValue(legacy), expected(axis));
        }
      };
      for (size_t update = 0; update < 2; ++update)
      {
        for (size_t a = 0; a < field.getSize(); ++a)
        {
          if constexpr (std::is_same_v<Scalar, Complex>)
            field[a] = Scalar((update + 1) * (a + 1), a + 2 * update + 1);
          else
            field[a] = Scalar((update + 1) * (a + 1));
        }
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const auto& formula =
            QF::PolytopeQuadratureFormula::get(4, cell->getGeometry());
          const PolytopeQuadrature quadrature(*cell, formula);
          for (size_t q = 0; q < quadrature.getSize(); ++q)
            checkPoint(IntegrationPoint(quadrature, q));
        }
      }
      const size_t dimension = mesh.getDimension();
      if (dimension == 0)
        return;
      mesh.getConnectivity().compute(dimension - 1, dimension);
      for (auto cell = mesh.getCell(); cell; ++cell)
        mesh.setAttribute(cell.key(), Geometry::Attribute(cell->getIndex() + 1));
      bool boundary = false, interior = false;
      for (auto face = mesh.getFace(); face; ++face)
      {
        const auto& formula = QF::PolytopeQuadratureFormula::get(2, face->getGeometry());
        const PolytopeQuadrature quadrature(*face, formula);
        const IntegrationPoint point(quadrature, 0);
        if (face->isBoundary())
        {
          checkPoint(point);
          boundary = true;
        }
        else
        {
          const auto& incident = mesh.getConnectivity().getIncidence(
            {dimension - 1, dimension}, face->getIndex());
          for (size_t axis = 0; axis < mesh.getSpaceDimension(); ++axis)
          {
            EXPECT_THROW(Derivative(axis, field).getValue(point), Alert::Exception);
            EXPECT_THROW(
              Derivative(axis, field).traceOf(0).getValue(point), Alert::Exception);
            for (const auto cell : incident)
            {
              const auto attribute = mesh.getPolytope(dimension, cell)->getAttribute();
              ASSERT_TRUE(attribute);
              auto partial = Derivative(axis, field);
              auto gradient = Grad(field);
              partial.traceOf(*attribute);
              gradient.traceOf(*attribute);
              const auto expected = gradient.getValue(point.getPoint());
              EXPECT_EQ(partial.getValue(point.getPoint()), expected(axis));
              EXPECT_EQ(partial.getValue(point), expected(axis));
            }
          }
          interior = true;
        }
      }
      EXPECT_TRUE(boundary);
      EXPECT_TRUE(interior);
    }
  }

  class P1DerivativeTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(P1DerivativeTest, ScalarPartialsAndTraces)
  {
    auto mesh = Convergence::UniformGrid(GetParam()).makeMesh(3);
    checkPartials<Real>(mesh);
    checkPartials<Complex>(mesh);
  }

  TEST(P1DerivativeEmbedded, SegmentInThreeDimensions)
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(3)
                  .nodes(3)
                  .vertex({0, 0, 0})
                  .vertex({1, 2, 2})
                  .vertex({2, 4, 4})
                  .polytope(Polytope::Type::Segment, {{0, 1}})
                  .polytope(Polytope::Type::Segment, {{1, 2}})
                  .finalize();
    checkPartials<Real>(mesh);
    checkPartials<Complex>(mesh);
  }

  TEST(P1DerivativeEmbedded, PointInThreeDimensions)
  {
    auto mesh = LocalMesh::Builder().initialize(3).nodes(1).vertex({1, 2, 3}).finalize();
    checkPartials<Real>(mesh);
    checkPartials<Complex>(mesh);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, P1DerivativeTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge));
}
