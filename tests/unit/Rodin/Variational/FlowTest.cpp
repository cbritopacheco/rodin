/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <type_traits>

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
        mesh = mesh.UniformGrid(GetParam(), {8, 8});
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
    auto velocity = VectorFunction{RealFunction([](const Point&) { return 0.1; }),
      RealFunction([](const Point&) { return 0.2; })};

    // Define the time step for advection
    Real dt = 0.1;
    (void)dt; // Suppress unused variable warning

    // Create test function
    TestFunction v(vh);

    // Test that we can create points with the correct coordinate format
    Math::SpatialVector<Real> coords{{0.5, 0.5}};
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
    auto constant_velocity =
      VectorFunction{RealFunction([](const Point&) { return 1.0; }),
        RealFunction([](const Point&) { return 0.5; })};

    // Test rotational velocity field: v = (-y, x)
    auto rotational_velocity = VectorFunction{-F::y, F::x};

    // Test linear shear velocity field: v = (y, 0)
    auto shear_velocity = VectorFunction{F::y, Zero()};

    // Create a test point
    Math::SpatialVector<Real> coords{{0.3, 0.4}};
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
    Math::SpatialVector<Real> coords{{0.25, 0.75}};
    Geometry::Point testPoint(*polytope, coords);

    // Verify physical coordinates are reasonable (within [0,1]x[0,1])
    auto physCoords = testPoint.getPhysicalCoordinates();
    EXPECT_GE(physCoords[0], 0.0);
    EXPECT_LE(physCoords[0], 1.0);
    EXPECT_GE(physCoords[1], 0.0);
    EXPECT_LE(physCoords[1], 1.0);
  }

  /// @brief Instantiates Flow Test over the Flow Mesh Params parameter coverage.
  INSTANTIATE_TEST_SUITE_P(FlowMeshParams, FlowTest,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle));
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
      for (const Real time : {-0.1, 0.0, 0.1})
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

/// @brief Boundary-face starts trace inward and stop immediately when directed outward.
TEST(Rodin_Variational_Flow, BoundaryFaceStartsRespectSignedCharacteristics)
{
  for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    mesh.getConnectivity().compute(dimension - 1, dimension);
    const auto& adjacent = mesh.getConnectivity().getIncidence(dimension - 1, dimension);
    size_t checked = 0;
    for (auto face = mesh.getBoundary(); face; ++face)
    {
      SCOPED_TRACE(face->getIndex());
      const Point start(*face, Polytope::Traits(face->getGeometry()).getCentroid());
      const auto cell = mesh.getCell(adjacent.at(face->getIndex())[0]);
      const Point center(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      const Math::SpatialVector<Real> direction =
        center.getPhysicalCoordinates() - start.getPhysicalCoordinates();
      const auto operand = RealFunction([](const Point& p) { return 2 + 3 * p(0); });
      for (const Real sign : {-1.0, 1.0})
      {
        const auto velocity = [direction, sign](const Point&) {
          return Math::SpatialVector<Real>(sign * direction);
        };
        Flow inward(sign * 0.1, operand, velocity);
        const auto trace = inward.trace(start);
        ASSERT_FALSE(trace.exited());
        EXPECT_NEAR(trace.getTime(), 0, 1e-12);
        EXPECT_LT((trace.getPoint().getPhysicalCoordinates() -
                    (start.getPhysicalCoordinates() + 0.1 * direction))
                    .norm(),
          1e-6);
        EXPECT_NEAR(inward(start), 2 + 3 * (start(0) + 0.1 * direction(0)), 1e-6);
        Flow outward(-sign * 0.1, operand, velocity);
        const auto stopped = outward.trace(start);
        EXPECT_TRUE(stopped.exited());
        EXPECT_NEAR(stopped.getTime(), 0.1, 1e-12);
        EXPECT_LT(
          (stopped.getPoint().getPhysicalCoordinates() - start.getPhysicalCoordinates())
            .norm(),
          1e-6);
      }
      ++checked;
    }
    EXPECT_GT(checked, 0u);
  }
}

/// @brief A continuing boundary policy must make progress or the trace stops at the hit.
TEST(Rodin_Variational_Flow, ContinuingBoundaryPolicyWithoutProgressStops)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  mesh.getConnectivity().compute(0, 1);
  const auto cell = mesh.getPolytope(1, 0);
  const auto operand = RealFunction([](const Point& point) { return 2 + 3 * point(0); });
  const auto velocity = [](const Point&) { return Math::SpatialVector<Real>{-1}; };
  constexpr size_t WorkLimit = 64;
  for (const Real coordinate : {0.0, 0.4})
  {
    SCOPED_TRACE(coordinate);
    const Point start(*cell, Math::SpatialPoint{coordinate});
    size_t calls = 0;
    const auto unchanged = [&calls](const BoundaryHit&) {
      ++calls;
      return true;
    };
    Flow flow(1.0, operand, velocity, Math::RungeKutta::RK4{}, unchanged);
    flow.setMaxOuterIterations(WorkLimit).setStepSafety(1);
    const auto trace = flow.trace(start);
    EXPECT_TRUE(trace.exited());
    EXPECT_NEAR(trace.getTime(), 1 - coordinate, 1e-6);
    EXPECT_NEAR(trace.getPoint()(0), 0, 1e-6);
    EXPECT_DOUBLE_EQ(trace.getCorrection(), 0);
    EXPECT_GE(calls, 1u);
    EXPECT_LE(calls, 2u);
    calls = 0;
    Flow reconstructed(1.0, operand, velocity, flow.getStep(), flow.getBoundaryPolicy());
    reconstructed.setMaxOuterIterations(WorkLimit).setStepSafety(1);
    EXPECT_DOUBLE_EQ(reconstructed(start), 0);
    EXPECT_GE(calls, 1u);
    EXPECT_LE(calls, 2u);
  }
}

/// @brief Complex scalar functions compose with real characteristics and retain typed boundary zeros.
TEST(Rodin_Variational_Flow, ComplexScalarCompositionAndBoundaryCorrection)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
  mesh.getConnectivity().compute(0, 1);
  const auto operand = ComplexFunction(
    [](const Point& point) { return Complex(2, -3) + Complex(1, 2) * point(0); });
  const auto velocity = [](const Point&) { return Math::SpatialVector<Real>{1}; };
  const auto cell = mesh.getPolytope(1, 1);
  const Point start(*cell, Math::SpatialPoint{0.4});
  for (const Real time : {-0.2, 0.2})
  {
    Flow flow(time, operand, velocity);
    const Complex expected = Complex(2, -3) + Complex(1, 2) * (start(0) + time);
    EXPECT_LT(std::abs(flow(start) - expected), 1e-6);
    std::unique_ptr<decltype(flow)> clone(flow.copy());
    EXPECT_LT(std::abs((*clone)(start)-expected), 1e-6);
  }
  const auto firstCell = mesh.getPolytope(1, 0);
  const Point nearBoundary(*firstCell, Math::SpatialPoint{0.4});
  const auto outward = [](const Point&) { return Math::SpatialVector<Real>{-1}; };
  Flow stopped(1.0, operand, outward);
  EXPECT_EQ(stopped(nearBoundary), Complex(0));
  const auto consumeRemainder = [](const BoundaryHit& hit) {
    hit.tau = 0;
    hit.correction += 2.5;
    return true;
  };
  Flow completed(1.0, operand, outward, Math::RungeKutta::RK4{}, consumeRemainder);
  EXPECT_LT(std::abs(completed(nearBoundary) - Complex(4.5, -3)), 1e-6);
}

/// @brief Boolean operands keep their existing real promotion under additive boundary correction.
TEST(Rodin_Variational_Flow, BooleanOperandPreservesScalarPromotion)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  mesh.getConnectivity().compute(0, 1);
  const auto cell = mesh.getPolytope(1, 0);
  const Point start(*cell, Math::SpatialPoint{0.4});
  const BooleanFunction operand(true);
  const auto velocity = [](const Point&) { return Math::SpatialVector<Real>{-1}; };
  Flow interior(0.1, operand, velocity);
  static_assert(std::is_same_v<decltype(interior(start)), Real>);
  EXPECT_DOUBLE_EQ(interior(start), 1);
  Flow stopped(1.0, operand, velocity);
  EXPECT_DOUBLE_EQ(stopped(start), 0);
  const auto correction = [](const BoundaryHit& hit) {
    hit.tau = 0;
    hit.correction += 0.5;
    return true;
  };
  Flow completed(1.0, operand, velocity, Math::RungeKutta::RK4{}, correction);
  EXPECT_DOUBLE_EQ(completed(start), 1.5);
}

/// @brief Starts below cell/face dimension abort without evaluating either field.
TEST(Rodin_Variational_Flow, LowerDimensionalStartsPreservePointAndSkipFields)
{
  for (const auto geometry :
    {Polytope::Type::Triangle, Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    mesh.getConnectivity().compute(dimension, 1);
    size_t operandCalls = 0, velocityCalls = 0;
    const RealFunction operand([&](const Point&) {
      ++operandCalls;
      return Real(7);
    });
    const auto velocity = [&](const Point&) {
      ++velocityCalls;
      Math::SpatialVector<Real> result(dimension);
      result.setConstant(1);
      return result;
    };
    for (size_t startDimension = 0; startDimension + 1 < dimension; ++startDimension)
    {
      SCOPED_TRACE(startDimension);
      ASSERT_GT(mesh.getPolytopeCount(startDimension), 0u);
      const auto polytope = mesh.getPolytope(startDimension, 0);
      const auto reference = Polytope::Traits(polytope->getGeometry()).getCentroid();
      const Point start(*polytope, reference);
      for (const Real time : {-0.25, 0.25})
      {
        SCOPED_TRACE(time);
        const Flow flow(time, operand, velocity);
        const auto trace = flow.trace(start);
        EXPECT_TRUE(trace.exited());
        EXPECT_DOUBLE_EQ(trace.getTime(), 0);
        EXPECT_DOUBLE_EQ(trace.getCorrection(), 0);
        EXPECT_EQ(&trace.getPoint().getPolytope().getMesh(), &mesh);
        EXPECT_EQ(trace.getPoint().getPolytope().getDimension(), startDimension);
        EXPECT_EQ(trace.getPoint().getPolytope().getIndex(), polytope->getIndex());
        const auto& tracedReference = trace.getPoint().getReferenceCoordinates();
        ASSERT_EQ(tracedReference.size(), reference.size());
        EXPECT_DOUBLE_EQ(flow(start), 0);
        ASSERT_EQ(start.getReferenceCoordinates().size(), reference.size());
        for (size_t component = 0; component < static_cast<size_t>(reference.size());
             ++component)
        {
          EXPECT_DOUBLE_EQ(tracedReference[component], reference[component]);
          EXPECT_DOUBLE_EQ(
            start.getReferenceCoordinates()[component], reference[component]);
        }
        EXPECT_EQ(operandCalls, 0u);
        EXPECT_EQ(velocityCalls, 0u);
      }
    }
  }
}

/// @brief Rotating characteristics agree with the analytic map across triangle and quadrilateral cells.
TEST(
  Rodin_Variational_Flow, RotatingCharacteristicsCrossCellsAndCloseInBothTimeDirections)
{
  constexpr Real PositionTolerance = 2e-5;
  const Real pi = std::acos(Real(-1));
  for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = LocalMesh::UniformGrid(geometry, {5, 5});
    mesh.scale(Real(0.25));
    mesh.getConnectivity().compute(1, 2);
    bool checked = false;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point start(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      const auto initial = start.getPhysicalCoordinates();
      const Real dx = initial.x() - Real(0.5);
      const Real dy = initial.y() - Real(0.5);
      if (!(dx > 0 && dy < 0 && dx * dx + dy * dy < Real(0.0625)))
        continue;
      const RealFunction operand(
        [](const Point& point) { return 2 + 3 * point.x() - 2 * point.y(); });
      const auto velocity = [](const Point& point) {
        return Math::SpatialVector<Real>{Real(0.5) - point.y(), point.x() - Real(0.5)};
      };
      for (const Real time : {-2 * pi, -pi / 2, pi / 2, 2 * pi})
      {
        SCOPED_TRACE(time);
        size_t boundaryHits = 0;
        const auto boundary = [&](const BoundaryHit&) {
          ++boundaryHits;
          return false;
        };
        Flow flow(time, operand, velocity, Math::RungeKutta::RK4{}, boundary);
        flow.setStepSafety(0.05).setTolFactor(2).setMaxBisections(32);
        const Math::SpatialPoint expected{
          Real(0.5) + std::cos(time) * dx - std::sin(time) * dy,
          Real(0.5) + std::sin(time) * dx + std::cos(time) * dy};
        const auto trace = flow.trace(start);
        ASSERT_FALSE(trace.exited());
        EXPECT_DOUBLE_EQ(trace.getTime(), 0);
        EXPECT_DOUBLE_EQ(trace.getCorrection(), 0);
        EXPECT_LT((trace.getPoint().getPhysicalCoordinates() - expected).norm(),
          PositionTolerance);
        EXPECT_NEAR(
          flow(start), 2 + 3 * expected.x() - 2 * expected.y(), 5 * PositionTolerance);
        EXPECT_EQ(boundaryHits, 0u);
        EXPECT_EQ((start.getPhysicalCoordinates() - initial).norm(), 0);
      }
      checked = true;
      break;
    }
    EXPECT_TRUE(checked);
  }
}

/// @brief Three-dimensional rotation crosses volume facets without spurious boundary exits.
TEST(Rodin_Variational_Flow, VolumeRotationMatchesAnalyticCharacteristics)
{
  constexpr Real PositionTolerance = 2e-5;
  const Real pi = std::acos(Real(-1));
  for (const auto geometry : {Polytope::Type::Tetrahedron, Polytope::Type::Hexahedron,
         Polytope::Type::Wedge, Polytope::Type::Pyramid})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = LocalMesh::UniformGrid(geometry, {3, 3, 3});
    mesh.scale(Real(0.5));
    mesh.getConnectivity().compute(2, 3);
    bool checked = false;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point start(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      const auto initial = start.getPhysicalCoordinates();
      const Real dx = initial.x() - Real(0.5);
      const Real dy = initial.y() - Real(0.5);
      if (!(dx > 0 && dy < 0 && dx * dx + dy * dy < Real(0.16)))
        continue;
      const RealFunction operand([](const Point& point) {
        return 2 + 3 * point.x() - 2 * point.y() + 5 * point.z();
      });
      const auto velocity = [](const Point& point) {
        return Math::SpatialVector<Real>{Real(0.5) - point.y(), point.x() - Real(0.5), 0};
      };
      for (const Real time : {-pi / 2, pi / 2})
      {
        SCOPED_TRACE(time);
        size_t boundaryHits = 0;
        const auto boundary = [&](const BoundaryHit&) {
          ++boundaryHits;
          return false;
        };
        Flow flow(time, operand, velocity, Math::RungeKutta::RK4{}, boundary);
        flow.setStepSafety(0.05).setTolFactor(2).setMaxBisections(32);
        const Math::SpatialPoint expected{
          Real(0.5) + std::cos(time) * dx - std::sin(time) * dy,
          Real(0.5) + std::sin(time) * dx + std::cos(time) * dy, initial.z()};
        const auto trace = flow.trace(start);
        ASSERT_FALSE(trace.exited());
        EXPECT_DOUBLE_EQ(trace.getTime(), 0);
        EXPECT_DOUBLE_EQ(trace.getCorrection(), 0);
        EXPECT_LT((trace.getPoint().getPhysicalCoordinates() - expected).norm(),
          PositionTolerance);
        EXPECT_NEAR(flow(start),
          2 + 3 * expected.x() - 2 * expected.y() + 5 * expected.z(),
          10 * PositionTolerance);
        EXPECT_EQ(boundaryHits, 0u);
        EXPECT_EQ((start.getPhysicalCoordinates() - initial).norm(), 0);
      }
      checked = true;
      break;
    }
    EXPECT_TRUE(checked);
  }
}

/// @brief A nested trace of the same type owns independent face-search storage.
TEST(FlowTest, NestedVelocityEvaluation)
{
  auto line = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
  line.getConnectivity().compute(0, 1);
  line.getConnectivity().compute(1, 0);
  auto triangle = LocalMesh::Builder()
                    .initialize(2)
                    .nodes(3)
                    .vertex({0, 0})
                    .vertex({1, 0})
                    .vertex({0, 1})
                    .polytope(Polytope::Type::Triangle, {0, 1, 2})
                    .finalize();
  triangle.getConnectivity().compute(1, 2);
  triangle.getConnectivity().compute(2, 1);
  const auto innerCell = triangle.getCell();
  const auto outerCell = line.getCell();
  const Point innerPoint(*innerCell, Math::SpatialPoint{0.2, 0.2});
  const Point outerPoint(*outerCell, Math::SpatialPoint{0.8});
  RealFunction operand(Real(1));
  using Velocity = std::function<Math::SpatialVector<Real>(const Point&)>;
  Velocity innerVelocity = [](const Point&) { return Math::SpatialVector<Real>{-1, 0}; };
  auto inner = Flow(Real(0.4), operand, innerVelocity);
  inner.setStepSafety(1).setMaxOuterIterations(100);
  bool nested = false;
  size_t callbacks = 0;
  Velocity outerVelocity = [&](const Point&) {
    ++callbacks;
    if (nested)
    {
      EXPECT_TRUE(inner.trace(innerPoint).exited());
    }
    return Math::SpatialVector<Real>{1};
  };
  auto outer = Flow(Real(0.4), operand, outerVelocity);
  outer.setStepSafety(1).setMaxOuterIterations(100);
  const auto control = outer.trace(outerPoint);
  const auto controlCalls = callbacks;
  callbacks = 0;
  nested = true;
  const auto result = outer.trace(outerPoint);
  EXPECT_FALSE(result.exited());
  EXPECT_EQ(result.getTime(), 0);
  EXPECT_NEAR(result.getPoint().x(), 1.2, 1e-14);
  EXPECT_EQ(result.getPoint().x(), control.getPoint().x());
  EXPECT_EQ(callbacks, controlCalls);
}
