/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Advection/Lagrangian.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Advection;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Unit tests for Flow class basic functionality
   */
  class FlowTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      Mesh<Context::Local> getMesh()
      {
        Mesh mesh;
        mesh = mesh.UniformGrid(GetParam(), { 8, 8 });
        mesh.scale(1.0 / 7.0);
        mesh.getConnectivity().compute(1, 2);
        return mesh;
      }
  };

  /**
   * @brief Test basic Flow class construction and coordinate handling.
   *
   * This test verifies that Flow objects can be constructed with valid
   * velocity fields and that basic coordinate system handling works.
   */
  TEST_P(FlowTest, BasicConstruction)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    // Define constant velocity field: v = (0.1, 0.2)
    auto velocity = VectorFunction{
      RealFunction([](const Point&) { return 0.1; }),
      RealFunction([](const Point&) { return 0.2; })
    };

    // Define the time step for advection
    Real dt = 0.1;
    (void)dt; // Suppress unused variable warning

    // Create test function
    TestFunction v(vh);

    // Test that we can create points with the correct coordinate format
    Math::SpatialVector<Real> coords{{ 0.5, 0.5 }};
    EXPECT_EQ(coords.size(), 2);
    EXPECT_NEAR(coords[0], 0.5, 1e-10);
    EXPECT_NEAR(coords[1], 0.5, 1e-10);

    // Test that we can create a point on the mesh
    auto polytope = mesh.getPolytope(2, 0);
    EXPECT_TRUE(polytope != mesh.getPolytope(2).end());

    // Create the test point 
    Geometry::Point testPoint(*polytope, coords);

    // Test velocity evaluation
    auto vel_value = velocity(testPoint);
    EXPECT_NEAR(vel_value[0], 0.1, 1e-10);
    EXPECT_NEAR(vel_value[1], 0.2, 1e-10);
  }

  /**
   * @brief Test different velocity field types and their evaluation.
   *
   * This test verifies that different types of velocity fields can be
   * properly constructed and evaluated at specific points.
   */
  TEST_P(FlowTest, VelocityFieldEvaluation)
  {
    Mesh mesh = this->getMesh();

    // Test constant velocity field
    auto constant_velocity = VectorFunction{
      RealFunction([](const Point&) { return 1.0; }),
      RealFunction([](const Point&) { return 0.5; })
    };

    // Test rotational velocity field: v = (-y, x)
    auto rotational_velocity = VectorFunction{
      -F::y,
      F::x
    };

    // Test linear shear velocity field: v = (y, 0)
    auto shear_velocity = VectorFunction{
      F::y,
      Zero()
    };

    // Create a test point
    Math::SpatialVector<Real> coords{{ 0.3, 0.4 }};
    auto polytope_iter = mesh.getPolytope(2, 0);
    Geometry::Point testPoint(*polytope_iter, coords);

    // Get the actual physical coordinates of the point
    auto physCoords = testPoint.getPhysicalCoordinates();
    Real actual_x = physCoords[0];
    Real actual_y = physCoords[1];

    // Evaluate constant velocity - should be (1.0, 0.5) everywhere
    auto const_vel_value = constant_velocity(testPoint);
    EXPECT_NEAR(const_vel_value[0], 1.0, 1e-10);
    EXPECT_NEAR(const_vel_value[1], 0.5, 1e-10);

    // Evaluate rotational velocity at the actual physical coordinates
    auto rot_vel_value = rotational_velocity(testPoint);
    EXPECT_NEAR(rot_vel_value[0], -actual_y, 1e-3);
    EXPECT_NEAR(rot_vel_value[1], actual_x, 1e-3);

    // Evaluate shear velocity at the actual physical coordinates
    auto shear_vel_value = shear_velocity(testPoint);
    EXPECT_NEAR(shear_vel_value[0], actual_y, 1e-3);
    EXPECT_NEAR(shear_vel_value[1], 0.0, 1e-10);
  }

  /**
   * @brief Test mesh polytope access and validation.
   *
   * This test verifies basic mesh functionality needed for Flow operations.
   */
  TEST_P(FlowTest, MeshAccess)
  {
    Mesh mesh = this->getMesh();

    // Check that mesh has elements
    EXPECT_GT(mesh.getCellCount(), 0);

    // Test polytope iteration
    auto polytope = mesh.getPolytope(2, 0);
    EXPECT_TRUE(polytope != mesh.getPolytope(2).end());

    // Test coordinate creation and point construction
    Math::SpatialVector<Real> coords{{ 0.25, 0.75 }};
    Geometry::Point testPoint(*polytope, coords);

    // Verify physical coordinates are reasonable (within [0,1]x[0,1])
    auto physCoords = testPoint.getPhysicalCoordinates();
    EXPECT_GE(physCoords[0], 0.0);
    EXPECT_LE(physCoords[0], 1.0);
    EXPECT_GE(physCoords[1], 0.0);
    EXPECT_LE(physCoords[1], 1.0);
  }

  /// @brief Instantiates Flow Test over the Flow Mesh Params parameter coverage.
  INSTANTIATE_TEST_SUITE_P(
    FlowMeshParams,
    FlowTest,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}

TEST(Rodin_Variational_Flow, ConstantCharacteristicCrossesCellsInBothDirections)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
  mesh.getConnectivity().compute(0, 1);
  const auto cell = mesh.getPolytope(1, 1);
  const Point start(*cell, Math::SpatialPoint{0.4});
  const auto operand = RealFunction([](const Point& p) { return 2 + 3 * p(0); });
  const auto velocity = [](const Point&) { return Math::SpatialVector<Real>{1}; };
  for (Real time : {-1.1, 0.0, 1.2})
  {
    Flow flow(time, operand, velocity);
    flow.setMaxZeroHops(64)
      .setMaxBisections(32)
      .setMaxSubdivisions(12)
      .setMaxOuterIterations(10000)
      .setMaxStagnations(256)
      .setFaceTrials(8)
      .setBracketGrowth(1.5)
      .setMaxBracketExpansions(5)
      .setStepSafety(0.4)
      .setMinStepFactor(8)
      .setTolFactor(40)
      .setClampFactor(8);
    const auto trace = flow.trace(start);
    ASSERT_FALSE(trace.exited());
    EXPECT_NEAR(trace.getTime(), 0, 1e-12);
    EXPECT_NEAR(trace.getCorrection(), 0, 1e-12);
    EXPECT_NEAR(trace.getPoint()(0), start(0) + time, 1e-7);
    EXPECT_NEAR(flow(start), 2 + 3 * (start(0) + time), 1e-6);
    EXPECT_FALSE(flow.getOrder(*cell).has_value());
    EXPECT_DOUBLE_EQ(flow.getOperand()(start), operand(start));
    EXPECT_DOUBLE_EQ(flow.getVelocity()(start)(0), 1);
    std::unique_ptr<decltype(flow)> clone(flow.copy());
    auto moved = std::move(*clone);
    EXPECT_NEAR(moved(start), flow(start), 1e-12);
  }
}

TEST(Rodin_Variational_Flow, BoundaryStopCorrectionAndWorkLimit)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
  mesh.getConnectivity().compute(0, 1);
  const auto cell = mesh.getPolytope(1, 0);
  const Point start(*cell, Math::SpatialPoint{0.4});
  const auto operand = RealFunction([](const Point& p) { return 2 + 3 * p(0); });
  const auto velocity = [](const Point&) { return Math::SpatialVector<Real>{-1}; };
  Flow stopped(1.0, operand, velocity);
  auto trace = stopped.trace(start);
  EXPECT_TRUE(trace.exited());
  EXPECT_NEAR(trace.getPoint()(0), 0, 1e-7);
  EXPECT_NEAR(trace.getTime(), 0.6, 1e-7);
  EXPECT_DOUBLE_EQ(stopped(start), 0);
  const auto correction = [](const BoundaryHit& hit) {
    hit.correction += 2.5;
    return false;
  };
  Flow corrected(1.0, operand, velocity, Math::RungeKutta::RK4{}, correction);
  const auto correctedTrace = corrected.trace(start);
  EXPECT_TRUE(correctedTrace.exited());
  EXPECT_NEAR(correctedTrace.getCorrection(), 2.5, 1e-12);
  EXPECT_DOUBLE_EQ(corrected(start), 0);
  const auto consumeRemainder = [](const BoundaryHit& hit) {
    hit.tau = 0;
    hit.correction += 2.5;
    return true;
  };
  Flow completed(1.0, operand, velocity, Math::RungeKutta::RK4{}, consumeRemainder);
  const auto completedTrace = completed.trace(start);
  EXPECT_FALSE(completedTrace.exited());
  EXPECT_DOUBLE_EQ(completedTrace.getTime(), 0);
  EXPECT_NEAR(completed(start), 4.5, 1e-6);
  corrected.setMaxOuterIterations(0);
  const auto limited = corrected.trace(start);
  EXPECT_TRUE(limited.exited());
  EXPECT_DOUBLE_EQ(limited.getTime(), 1);
  EXPECT_NEAR(limited.getPoint()(0), start(0), 1e-12);
}

TEST(Rodin_Variational_Flow, StationaryCharacteristicsPreserveEveryCellGeometry)
{
  for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    mesh.getConnectivity().compute(dimension - 1, dimension);
    const auto cell = mesh.getPolytope(dimension, 0);
    const Point start(*cell, Polytope::Traits(geometry).getCentroid());
    const auto operand = RealFunction([](const Point& p) { return 1 + p(0); });
    // Matrix-valued velocity exercises the column extraction specialization.
    const auto velocity = [dimension](const Point&) {
      Math::SpatialMatrix<Real> value(dimension, 1);
      value.setZero();
      return value;
    };
    Flow stationary(0.3, operand, velocity);
    const auto trace = stationary.trace(start);
    EXPECT_FALSE(trace.exited());
    EXPECT_NEAR(trace.getTime(), 0, 1e-12);
    EXPECT_LT(
      (trace.getPoint().getPhysicalCoordinates() - start.getPhysicalCoordinates()).norm(),
      1e-12);
    EXPECT_NEAR(stationary(start), operand(start), 1e-12);
  }
}

TEST(Rodin_Variational_Flow, InteriorFaceStartChoosesTheUpwindCell)
{
  for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 3));
    mesh.getConnectivity().compute(dimension - 1, dimension);
    const auto& adjacent = mesh.getConnectivity().getIncidence(dimension - 1, dimension);
    bool checked = false;
    for (auto face = mesh.getPolytope(dimension - 1); face; ++face)
    {
      const auto& cells = adjacent.at(face->getIndex());
      if (cells.size() != 2)
        continue;
      const Point start(*face, Polytope::Traits(face->getGeometry()).getCentroid());
      const auto cell = mesh.getPolytope(dimension, cells[0]);
      const Point center(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      const Math::SpatialVector<Real> direction =
        center.getPhysicalCoordinates() - start.getPhysicalCoordinates();
      const auto operand = RealFunction([](const Point& p) { return 1 + p(0); });
      const auto velocity = [direction](const Point&) { return direction; };
      for (const Real time : {-0.1, 0.1})
      {
        Flow flow(time, operand, velocity);
        const auto trace = flow.trace(start);
        ASSERT_FALSE(trace.exited());
        EXPECT_LT((trace.getPoint().getPhysicalCoordinates() -
                    (start.getPhysicalCoordinates() + time * direction))
                    .norm(),
          1e-6);
        EXPECT_NEAR(flow(start), 1 + start(0) + time * direction(0), 1e-6);
      }
      checked = true;
      break;
    }
    EXPECT_TRUE(checked);
  }
}
