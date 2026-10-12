/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Independent one-sided trace oracles for P1 and quadratic H1 operators.
 */
#include <gtest/gtest.h>
#include <cmath>
#include <type_traits>
#include <vector>
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief A continuous kink distinguishes the two one-sided derivative traces.
  template <class Scalar, bool Quadratic = false>
  void checkPiecewiseAffineTraces()
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
      for (auto cell = mesh.getCell(); cell; ++cell)
        mesh.setAttribute({dimension, cell->getIndex()}, cell->getIndex() + 1);
      auto scalarSpace = [&] {
        if constexpr (Quadratic)
          return H1<2, Scalar>(std::integral_constant<size_t, 2>{}, mesh);
        else
          return P1<Scalar>(mesh);
      }();
      auto vectorSpace = [&] {
        if constexpr (Quadratic)
          return H1<2, Math::SpatialVector<Scalar>>(
            std::integral_constant<size_t, 2>{}, mesh, dimension);
        else
          return P1<Math::SpatialVector<Scalar>>(mesh, dimension);
      }();
      GridFunction scalarField(scalarSpace);
      GridFunction vectorField(vectorSpace);
      scalarField = [factor, dimension](const Point& p) {
        Scalar value = p(0) <= 1 ? p(0) : 3 * p(0) - 2;
        for (size_t axis = 1; axis < dimension; ++axis)
          value += Scalar(axis + 1) * p(axis);
        return Scalar(factor * value);
      };
      vectorField = [factor, dimension](const Point& p) {
        const Real kink = p(0) <= 1 ? p(0) : 3 * p(0) - 2;
        Math::SpatialVector<Scalar> value(dimension);
        for (size_t axis = 0; axis < dimension; ++axis)
          value(axis) = factor * (Scalar(axis + 1) * kink + Scalar(0.5) * p(axis));
        return value;
      };
      size_t kinkFaces = 0;
      bool checkedMissingTrace = false;
      for (auto face = mesh.getFace(); face; ++face)
      {
        SCOPED_TRACE(face->getIndex());
        const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
        const auto& cells = mesh.getConnectivity().getIncidence(
          {dimension - 1, dimension}, face->getIndex());
        std::vector<Real> slopes;
        for (const auto index : cells)
        {
          const auto cell = mesh.getCell(index);
          const Point center(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
          const Real slope = center(0) < 1 ? 1 : 3;
          slopes.push_back(slope);
          auto gradient = Grad(scalarField);
          auto jacobian = Jacobian(vectorField);
          auto divergence = Div(vectorField);
          gradient.traceOf(index + 1);
          jacobian.traceOf(index + 1);
          divergence.traceOf(index + 1);
          const auto g = gradient.getValue(point);
          const auto j = jacobian.getValue(point);
          for (size_t row = 0; row < dimension; ++row)
          {
            EXPECT_LT(
              std::abs(g(row) - factor * Scalar(row == 0 ? slope : row + 1)), 1e-10);
            for (size_t column = 0; column < dimension; ++column)
            {
              const Scalar expected = factor *
                Scalar((column == 0 ? (row + 1) * slope : 0) + (row == column ? 0.5 : 0));
              EXPECT_LT(std::abs(j(row, column) - expected), 1e-10);
            }
          }
          EXPECT_LT(std::abs(divergence.getValue(point) -
                      factor * Scalar(slope + 0.5 * dimension)),
            1e-10);
        }
        if (slopes.size() == 2)
        {
          if (slopes[0] != slopes[1])
            ++kinkFaces;
          if (checkedMissingTrace)
            continue;
          auto gradient = Grad(scalarField);
          auto jacobian = Jacobian(vectorField);
          auto divergence = Div(vectorField);
          const Attribute absent = mesh.getCellCount() + 1;
          gradient.traceOf(absent);
          jacobian.traceOf(absent);
          divergence.traceOf(absent);
          EXPECT_THROW(gradient.getValue(point), Alert::Exception);
          EXPECT_THROW(jacobian.getValue(point), Alert::Exception);
          EXPECT_THROW(divergence.getValue(point), Alert::Exception);
          checkedMissingTrace = true;
        }
      }
      EXPECT_GT(kinkFaces, 0u);
      EXPECT_TRUE(checkedMissingTrace);
    }
  }
}

TEST(Rodin_Variational_DerivativeTrace, P1RealGradientJacobianAndDivergence)
{
  checkPiecewiseAffineTraces<Real>();
}

TEST(Rodin_Variational_DerivativeTrace, P1ComplexGradientJacobianAndDivergence)
{
  checkPiecewiseAffineTraces<Complex>();
}

TEST(Rodin_Variational_DerivativeTrace, H1QuadraticRealGradientJacobianAndDivergence)
{
  checkPiecewiseAffineTraces<Real, true>();
}

TEST(Rodin_Variational_DerivativeTrace, H1QuadraticComplexGradientJacobianAndDivergence)
{
  checkPiecewiseAffineTraces<Complex, true>();
}
