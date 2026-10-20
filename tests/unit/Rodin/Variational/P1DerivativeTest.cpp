/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Affine scalar derivative and order-dispatch oracles.
 */
#include <gtest/gtest.h>
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "../../../convergence/Convergence.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <class Scalar, bool Quadratic = false>
  void checkAffineDerivative()
  {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 3));
      mesh.getConnectivity().compute(dimension - 1, dimension);
      if constexpr (Quadratic)
      {
        for (size_t d = 1; d <= dimension; ++d)
        {
          mesh.getConnectivity().compute(dimension, d);
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
      }
      for (auto cell = mesh.getPolytope(dimension); cell; ++cell)
        mesh.setAttribute({dimension, cell->getIndex()}, cell->getIndex() + 1);
      auto space = [&] {
        if constexpr (Quadratic)
          return H1<2, Scalar>(std::integral_constant<size_t, 2>{}, mesh);
        else
          return P1<Scalar>(mesh);
      }();
      GridFunction field(space);
      const Scalar factor = [] {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(2, -3);
        else
          return Scalar(2);
      }();
      field = [factor, dimension](const Point& p) {
        Scalar value = 1;
        for (size_t axis = 0; axis < dimension; ++axis)
          value += factor * Scalar(axis + 1) * p(axis);
        return value;
      };
      for (size_t axis = 0; axis < dimension; ++axis)
      {
        Derivative derivative(axis, field);
        for (auto cell = mesh.getPolytope(dimension); cell; ++cell)
        {
          const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
          using Base = typename decltype(derivative)::Parent;
          using ScalarBase = typename Base::Parent;
          const size_t degreeFactor = geometry == Polytope::Type::Hexahedron ? 3
            : (geometry == Polytope::Type::Quadrilateral ||
                geometry == Polytope::Type::Wedge || geometry == Polytope::Type::Pyramid)
            ? 2
            : 1;
          const Optional<size_t> expectedOrder(degreeFactor * (Quadratic ? 2 : 1) - 1);
          EXPECT_EQ(derivative.getOrder(*cell), expectedOrder);
          EXPECT_EQ(static_cast<const Base&>(derivative).getOrder(*cell), expectedOrder);
          EXPECT_EQ(
            static_cast<const ScalarBase&>(derivative).getOrder(*cell), expectedOrder);
          EXPECT_LT(
            std::abs(derivative.getValue(point) - factor * Scalar(axis + 1)), 1e-10);
        }
        for (auto face = mesh.getPolytope(dimension - 1); face; ++face)
        {
          const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
          const auto& cells = mesh.getConnectivity().getIncidence(
            {dimension - 1, dimension}, face->getIndex());
          auto traced = derivative;
          traced.traceOf(cells[0] + 1);
          EXPECT_LT(std::abs(traced.getValue(point) - factor * Scalar(axis + 1)), 1e-10);
          std::unique_ptr<decltype(traced)> clone(traced.copy());
          auto moved = std::move(*clone);
          EXPECT_LT(std::abs(moved.getValue(point) - factor * Scalar(axis + 1)), 1e-10);
        }
      }
    }
  }
}

TEST(Rodin_Variational_P1Derivative, AffineRealAcrossCellsAndTraces)
{
  checkAffineDerivative<Real>();
}

TEST(Rodin_Variational_P1Derivative, AffineComplexAcrossCellsAndTraces)
{
  checkAffineDerivative<Complex>();
}

TEST(Rodin_Variational_H1Derivative, AffineRealAcrossCellsAndTraces)
{
  checkAffineDerivative<Real, true>();
}

TEST(Rodin_Variational_H1Derivative, AffineComplexAcrossCellsAndTraces)
{
  checkAffineDerivative<Complex, true>();
}

namespace
{
  /// @brief Parent/submesh transport preserves physical derivative values and integration-point dispatch.
  template <class Scalar, bool Quadratic>
  void checkTransportedDerivatives()
  {
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -3);
      else
        return Scalar(2);
    }();
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto parent =
        LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 3));
      for (auto cell = parent.getCell(); cell; ++cell)
      {
        const Point center(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
        parent.setAttribute({dimension, cell->getIndex()}, center(0) < 1 ? 17 : 23);
      }
      auto retained = parent.keep(17);
      auto exercise = [&](Mesh<Context::Local>& domain, bool useParentPoint) {
        for (size_t d = 1; d <= dimension; ++d)
        {
          domain.getConnectivity().compute(dimension, d);
          for (size_t lower = 0; lower < d; ++lower)
            domain.getConnectivity().compute(d, lower);
        }
        auto scalarSpace = [&] {
          if constexpr (Quadratic)
            return H1<2, Scalar>(std::integral_constant<size_t, 2>{}, domain);
          else
            return P1<Scalar>(domain);
        }();
        auto vectorSpace = [&] {
          if constexpr (Quadratic)
          {
            return H1<2, Math::SpatialVector<Scalar>>(
              std::integral_constant<size_t, 2>{}, domain, dimension);
          }
          else
            return P1<Math::SpatialVector<Scalar>>(domain, dimension);
        }();
        GridFunction scalarField(scalarSpace);
        GridFunction vectorField(vectorSpace);
        scalarField = [factor, dimension](const Point& p) {
          Scalar value = 1;
          for (size_t axis = 0; axis < dimension; ++axis)
            value += Scalar(axis + 1) * (Quadratic ? p(axis) * p(axis) : p(axis));
          return Scalar(factor * value);
        };
        vectorField = [factor, dimension](const Point& p) {
          Scalar value = 2;
          for (size_t axis = 0; axis < dimension; ++axis)
            value += Scalar(0.5 + axis) * (Quadratic ? p(axis) * p(axis) : p(axis));
          Math::SpatialVector<Scalar> result(dimension);
          for (size_t row = 0; row < dimension; ++row)
            result[row] = factor * Scalar(row + 1) * value;
          return result;
        };
        auto gradient = Grad(scalarField);
        auto jacobian = Jacobian(vectorField);
        auto divergence = Div(vectorField);
        const QF::PolytopeQuadratureFormula formula(4, geometry);
        size_t checked = 0;
        for (auto cell = retained.getCell(); cell; ++cell)
        {
          for (const size_t sample : {size_t(0), formula.getSize() - 1})
          {
            const Point childPoint(*cell, formula.getPoint(sample));
            const auto included = parent.inclusion(childPoint);
            ASSERT_TRUE(included);
            const Point& point = useParentPoint ? *included : childPoint;
            ASSERT_FALSE(domain.isLocalPoint(point));
            auto check = [&](const auto& context) {
              const auto g = gradient.getValue(context);
              const auto j = jacobian.getValue(context);
              Scalar expectedDivergence = 0;
              for (size_t axis = 0; axis < dimension; ++axis)
              {
                const Real derivative = Quadratic ? 2 * point(axis) : 1;
                const Scalar expectedPartial = factor * Scalar(axis + 1) * derivative;
                EXPECT_LT(std::abs(g[axis] - expectedPartial), 1e-10);
                Derivative partial(axis, scalarField);
                EXPECT_LT(std::abs(partial.getValue(context) - expectedPartial), 1e-10);
                for (size_t row = 0; row < dimension; ++row)
                {
                  const Scalar expected =
                    factor * Scalar(row + 1) * Scalar(0.5 + axis) * derivative;
                  EXPECT_LT(std::abs(j(row, axis) - expected), 1e-10);
                  if (row == axis)
                    expectedDivergence += expected;
                }
              }
              EXPECT_LT(
                std::abs(divergence.getValue(context) - expectedDivergence), 1e-10);
            };
            check(point);
            check(IntegrationPoint(point));
            check(IntegrationPoint(point, &formula, sample));
            const Point& localPoint = useParentPoint ? childPoint : *included;
            ASSERT_TRUE(domain.isLocalPoint(localPoint));
            check(localPoint);
            check(IntegrationPoint(localPoint));
            check(IntegrationPoint(localPoint, &formula, sample));
            ++checked;
          }
        }
        EXPECT_GT(checked, 0);
      };
      exercise(parent, false);
      exercise(retained, true);
    }
  }
}

TEST(Rodin_Variational_P1Derivative, TransportedRealDerivatives)
{
  checkTransportedDerivatives<Real, false>();
}

TEST(Rodin_Variational_P1Derivative, TransportedComplexDerivatives)
{
  checkTransportedDerivatives<Complex, false>();
}

TEST(Rodin_Variational_H1Derivative, TransportedRealDerivatives)
{
  checkTransportedDerivatives<Real, true>();
}

TEST(Rodin_Variational_H1Derivative, TransportedComplexDerivatives)
{
  checkTransportedDerivatives<Complex, true>();
}

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
