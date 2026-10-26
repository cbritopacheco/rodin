/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Analytic oblique characteristics through every supported cell geometry. */
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational/Flow.h"
#include "Rodin/Variational/RealFunction.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  class FlowCrossing : public ::testing::TestWithParam<Polytope::Type>
  {};

  /// @brief Constant physical velocity gives an exact line across internal faces and corners.
  TEST_P(FlowCrossing, ObliqueCharacteristicPreservesPositionAndAffineField)
  {
    const auto geometry = GetParam();
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    constexpr size_t PointsPerAxis = 5;
    auto mesh =
      LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, PointsPerAxis));
    mesh.scale(Real(1) / Real(PointsPerAxis - 1));
    mesh.getConnectivity().discover(dimension, dimension - 1);
    mesh.getConnectivity().restrict(dimension - 1, dimension);
    const auto reference = Polytope::Traits(geometry).getCentroid();
    Math::SpatialVector<Real> referenceVelocity(static_cast<std::uint8_t>(dimension));
    referenceVelocity.setConstant(1);
    const auto affine = [](const Math::SpatialPoint& point) {
      Real result = 1;
      for (size_t coordinate = 0; coordinate < point.size(); ++coordinate)
        result += Real(coordinate + 2) * point[coordinate];
      return result;
    };
    const auto operand = RealFunction(
      [affine](const Point& point) { return affine(point.getPhysicalCoordinates()); });
    // Unit reference velocity reaches a face from every centroid before this time.
    constexpr Real Duration = 0.6;
    constexpr Real CoordinateTolerance = 1e-6;
    // Refine face hits below the coordinate oracle before exhausting the work budget.
    constexpr size_t BisectionLimit = 48;
    constexpr Real ReferenceRoundoffMultiplier = 4;
    for (const Real sign : {-1.0, 1.0})
    {
      SCOPED_TRACE(sign);
      bool found = false;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const Point start(*cell, reference);
        const Math::SpatialVector<Real> direction =
          start.getJacobian() * referenceVelocity;
        const Math::SpatialPoint expected =
          start.getPhysicalCoordinates() + sign * Duration * direction;
        bool interior = true;
        for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
        {
          interior &= start(coordinate) > 0.3 && start(coordinate) < 0.7;
          interior &= expected[coordinate] > 0.05 && expected[coordinate] < 0.95;
        }
        if (!interior)
          continue;
        const auto velocity = [direction](const Point&) { return direction; };
        Flow flow(sign * Duration, operand, velocity);
        flow.setMaxBisections(BisectionLimit).setTolFactor(ReferenceRoundoffMultiplier);
        const auto trace = flow.trace(start);
        ASSERT_FALSE(trace.exited());
        EXPECT_NEAR(trace.getTime(), 0, 1e-12);
        EXPECT_LT((trace.getPoint().getPhysicalCoordinates() - expected).norm(),
          CoordinateTolerance);
        EXPECT_NE(trace.getPoint().getPolytope().getIndex(), cell->getIndex());
        EXPECT_NEAR(flow(start), affine(expected), CoordinateTolerance);
        found = true;
        break;
      }
      EXPECT_TRUE(found);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, FlowCrossing,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge));
}
