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
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"

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
