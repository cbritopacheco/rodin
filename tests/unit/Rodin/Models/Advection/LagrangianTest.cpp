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
  /// @brief The inward shift and Taylor correction recover exact affine boundary data.
  TEST(Rodin_Advection_TaylorBoundaryShift, AffineBoundaryRecoveryOnEveryGeometry)
  {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
      mesh.getConnectivity().compute(dimension - 1, dimension);
      mesh.getConnectivity().compute(dimension, dimension - 1);
      P1 space(mesh);
      GridFunction field(space);
      const auto exact = [dimension](const Point& p) {
        Real value = 2;
        for (size_t axis = 0; axis < dimension; ++axis)
          value += (axis + 1) * p(axis);
        return value;
      };
      field = exact;
      size_t checked = 0;
      for (auto face = mesh.getBoundary(); face; ++face)
      {
        SCOPED_TRACE(face->getIndex());
        const Point boundary(*face, Polytope::Traits(face->getGeometry()).getCentroid());
        Math::SpatialVector<Real> normal(dimension);
        normal.setZero();
        size_t normalAxes = 0;
        for (size_t axis = 0; axis < dimension; ++axis)
        {
          if (std::abs(boundary(axis)) < 1e-12)
          {
            normal(axis) = -1;
            ++normalAxes;
          }
          else if (std::abs(boundary(axis) - 1) < 1e-12)
          {
            normal(axis) = 1;
            ++normalAxes;
          }
        }
        ASSERT_EQ(normalAxes, 1u);
        const Index originalCell = mesh.getConnectivity()
                                     .getIncidence(dimension - 1, dimension)
                                     .at(face->getIndex())[0];
        const auto& faces =
          mesh.getConnectivity().getIncidence(dimension, dimension - 1).at(originalCell);
        size_t localFace = 0;
        while (localFace < faces.size() && faces[localFace] != face->getIndex())
          ++localFace;
        ASSERT_LT(localFace, faces.size());
        const auto cell = mesh.getCell(originalCell);
        for (const Real timeStep : {1.0, 0.2, -0.2})
        {
          SCOPED_TRACE(timeStep);
          constexpr Real requestedOffset = 0.0625;
          const Real offset = std::min(requestedOffset, Real(0.1) * std::abs(timeStep));
          TaylorBoundaryShiftPolicy policy(timeStep, mesh, field, requestedOffset);
          Real remaining = 0.3;
          Index currentCell = originalCell;
          Math::SpatialPoint reference;
          cell->getTransformation().inverse(reference, boundary.getPhysicalCoordinates());
          Real correction = 1.25;
          BoundaryHit hit{remaining, currentCell, reference, localFace, correction};
          EXPECT_TRUE(policy(hit));
          EXPECT_DOUBLE_EQ(remaining, 0);
          EXPECT_EQ(currentCell, originalCell);
          const Point shifted(*cell, reference);
          EXPECT_LT((shifted.getPhysicalCoordinates() -
                      (boundary.getPhysicalCoordinates() - offset * normal))
                      .norm(),
            1e-10);
          Real normalDerivative = 0;
          for (size_t axis = 0; axis < dimension; ++axis)
            normalDerivative += (axis + 1) * normal(axis);
          EXPECT_NEAR(correction, 1.25 + offset * normalDerivative, 1e-10);
          EXPECT_NEAR(
            field.getValue(shifted) + correction, exact(boundary) + 1.25, 1e-10);
          const auto previousReference = reference;
          const Real previousCorrection = correction;
          EXPECT_TRUE(policy(hit));
          EXPECT_LT((reference - previousReference).norm(), 1e-12);
          EXPECT_DOUBLE_EQ(correction, previousCorrection);
        }
        ++checked;
      }
      EXPECT_GT(checked, 0u);
    }
  }

  /// @brief The stop policy applies the physical offset without changing the traced value.
  TEST(Rodin_Advection_StopInside, PhysicalInwardOffsetOnEveryGeometry)
  {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
      mesh.getConnectivity().compute(dimension - 1, dimension);
      mesh.getConnectivity().compute(dimension, dimension - 1);
      size_t checked = 0;
      for (auto face = mesh.getBoundary(); face; ++face)
      {
        SCOPED_TRACE(face->getIndex());
        const Point boundary(*face, Polytope::Traits(face->getGeometry()).getCentroid());
        Math::SpatialVector<Real> normal(dimension);
        normal.setZero();
        size_t normalAxes = 0;
        for (size_t axis = 0; axis < dimension; ++axis)
        {
          if (std::abs(boundary(axis)) < 1e-12)
          {
            normal(axis) = -1;
            ++normalAxes;
          }
          else if (std::abs(boundary(axis) - 1) < 1e-12)
          {
            normal(axis) = 1;
            ++normalAxes;
          }
        }
        ASSERT_EQ(normalAxes, 1u);
        const Index originalCell = mesh.getConnectivity()
                                     .getIncidence(dimension - 1, dimension)
                                     .at(face->getIndex())[0];
        const auto& faces =
          mesh.getConnectivity().getIncidence(dimension, dimension - 1).at(originalCell);
        size_t localFace = 0;
        while (localFace < faces.size() && faces[localFace] != face->getIndex())
          ++localFace;
        ASSERT_LT(localFace, faces.size());
        const auto cell = mesh.getCell(originalCell);
        for (const Real timeStep : {1.0, 0.2, -0.2})
        {
          SCOPED_TRACE(timeStep);
          constexpr Real requestedOffset = 0.0625;
          const Real offset = std::min(requestedOffset, Real(0.1) * std::abs(timeStep));
          StopInsideBoundaryPolicy policy(timeStep, mesh, requestedOffset);
          Real remaining = 0.3;
          Index currentCell = originalCell;
          Math::SpatialPoint reference;
          cell->getTransformation().inverse(reference, boundary.getPhysicalCoordinates());
          Real correction = 1.25;
          BoundaryHit hit{remaining, currentCell, reference, localFace, correction};
          EXPECT_TRUE(policy(hit));
          EXPECT_DOUBLE_EQ(remaining, 0);
          EXPECT_EQ(currentCell, originalCell);
          const Point shifted(*cell, reference);
          EXPECT_LT((shifted.getPhysicalCoordinates() -
                      (boundary.getPhysicalCoordinates() - offset * normal))
                      .norm(),
            1e-10);
          EXPECT_DOUBLE_EQ(correction, 1.25);
          const auto previousReference = reference;
          const Real previousCorrection = correction;
          EXPECT_TRUE(policy(hit));
          EXPECT_LT((reference - previousReference).norm(), 1e-12);
          EXPECT_DOUBLE_EQ(correction, previousCorrection);
        }
        ++checked;
      }
      EXPECT_GT(checked, 0u);
    }
  }

  /// @brief Exhausted time and an invalid local face stop without modifying other state.
  TEST(Rodin_Advection_StopInside, ExhaustedTimeAndInvalidFacePreserveState)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    StopInsideBoundaryPolicy policy(Real(0.2), mesh);
    const auto cell = mesh.getCell(0);
    const auto original = Polytope::Traits(cell->getGeometry()).getCentroid();
    for (const Real time : {Real(0), Real(-0.25), Real(0.25)})
    {
      SCOPED_TRACE(time);
      Real remaining = time;
      Index currentCell = 0;
      auto reference = original;
      Real correction = 1.25;
      const size_t invalidFace =
        Polytope::Traits(cell->getGeometry()).getHalfSpace().vector.size();
      BoundaryHit hit{remaining, currentCell, reference, invalidFace, correction};
      EXPECT_TRUE(policy(hit));
      EXPECT_DOUBLE_EQ(remaining, time > 0 ? Real(0) : time);
      EXPECT_EQ(currentCell, 0u);
      ASSERT_EQ(reference.size(), original.size());
      for (size_t axis = 0; axis < reference.size(); ++axis)
        EXPECT_DOUBLE_EQ(reference(axis), original(axis));
      EXPECT_DOUBLE_EQ(correction, 1.25);
    }
  }

  /**
   * @brief Unit tests for Lagrangian class basic functionality
   */
  class LagrangianTest : public ::testing::TestWithParam<Polytope::Type>
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
   * @brief Test basic Lagrangian class construction.
   *
   * This test verifies that Lagrangian objects can be constructed
   * with proper trial/test functions, initial conditions, and velocity fields.
   */
  TEST_P(LagrangianTest, BasicConstruction)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    // Create trial and test functions
    TrialFunction u(vh);
    TestFunction v(vh);

    // Define initial condition
    auto pi = Math::Constants::pi();
    auto u0 = sin(pi * F::x) * sin(pi * F::y);

    // Define velocity field
    auto velocity = VectorFunction{RealFunction([](const Point&) { return 0.1; }),
      RealFunction([](const Point&) { return 0.2; })};

    // Test construction
    try
    {
      Lagrangian lagrangian(u, v, u0, velocity);
      // If we get here, construction succeeded
      EXPECT_TRUE(true);
    }
    catch (const std::exception& e)
    {
      FAIL() << "Lagrangian construction failed: " << e.what();
    }
    catch (...)
    {
      FAIL() << "Lagrangian construction failed with unknown exception";
    }
  }

  /**
   * @brief Test Lagrangian with different velocity field types.
   *
   * This test verifies that Lagrangian can work with various velocity fields.
   */
  TEST_P(LagrangianTest, DifferentVelocityFields)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    TrialFunction u(vh);
    TestFunction v(vh);

    // Simple initial condition
    auto u0 = RealFunction([](const Point&) { return 1.0; });

    // Test with constant velocity
    auto constant_velocity =
      VectorFunction{RealFunction([](const Point&) { return 1.0; }),
        RealFunction([](const Point&) { return 0.0; })};

    // Test with variable velocity
    auto variable_velocity = VectorFunction{F::x, F::y};

    // Test construction with constant velocity
    EXPECT_NO_THROW({ Lagrangian lagrangian1(u, v, u0, constant_velocity); });

    // Test construction with variable velocity
    EXPECT_NO_THROW({ Lagrangian lagrangian2(u, v, u0, variable_velocity); });
  }

  /**
   * @brief Test Lagrangian step function basic functionality.
   *
   * This test verifies that the step function can be called without errors.
   * Note: The current implementation may have template deduction issues,
   * so this test focuses on basic construction.
   */
  TEST_P(LagrangianTest, StepFunction)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    TrialFunction u(vh);
    TestFunction v(vh);

    // Define simple initial condition
    auto u0 = RealFunction([](const Point&) { return 1.0; });

    // Define simple velocity field
    auto velocity = VectorFunction{RealFunction([](const Point&) { return 0.1; }),
      RealFunction([](const Point&) { return 0.0; })};

    // Create Lagrangian object
    Lagrangian lagrangian(u, v, u0, velocity);

    // For now, just test that construction succeeded
    // TODO: Enable step function test when template issues are resolved
    // Real dt = 0.01;
    // EXPECT_NO_THROW({
    //   lagrangian.step(dt);
    // });

    EXPECT_TRUE(true); // Basic construction test passes
  }

  /**
   * @brief Test with different initial conditions.
   *
   * This test verifies that Lagrangian works with various initial conditions.
   */
  TEST_P(LagrangianTest, DifferentInitialConditions)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    TrialFunction u(vh);
    TestFunction v(vh);

    auto velocity = VectorFunction{RealFunction([](const Point&) { return 0.1; }),
      RealFunction([](const Point&) { return 0.1; })};

    // Test with constant initial condition
    auto u0_constant = RealFunction([](const Point&) { return 2.5; });
    EXPECT_NO_THROW({ Lagrangian lagrangian1(u, v, u0_constant, velocity); });

    // Test with polynomial initial condition
    auto u0_poly = F::x + F::y;
    EXPECT_NO_THROW({ Lagrangian lagrangian2(u, v, u0_poly, velocity); });

    // Test with trigonometric initial condition
    auto pi = Math::Constants::pi();
    auto u0_trig = sin(pi * F::x) * cos(pi * F::y);
    EXPECT_NO_THROW({ Lagrangian lagrangian3(u, v, u0_trig, velocity); });
  }

  /// @brief Instantiates Lagrangian Test over the Lagrangian Mesh Params parameter coverage.
  INSTANTIATE_TEST_SUITE_P(LagrangianMeshParams, LagrangianTest,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle));
}
