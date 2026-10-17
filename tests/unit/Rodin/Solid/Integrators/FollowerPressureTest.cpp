/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Exact affine loads and finite-difference follower-pressure tangents. */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Solid/Integrators/FollowerPressure.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Solid_FollowerPressure, AffineLoadAndEveryTangentColumn)
{
  // These integrators claim P1 segment/triangle faces, so test both supported
  // dimensions; higher-order and quadrilateral faces are outside that contract.
  for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    mesh.getConnectivity().compute(dimension - 1, dimension);
    P1 space(mesh, dimension);
    TrialFunction u(space);
    TestFunction v(space);
    GridFunction displacement(space);
    Math::SpatialMatrix<Real> gradient(dimension, dimension);
    for (size_t row = 0; row < dimension; ++row)
    {
      for (size_t column = 0; column < dimension; ++column)
        gradient(row, column) = 0.03 * (row + column + 1);
    }
    displacement = [gradient](
                     const Point& p) { return gradient * p.getPhysicalCoordinates(); };
    Real pressure = 2.3;
    Solid::FollowerPressureForce force(pressure, v, displacement);
    Solid::FollowerPressureTangent tangent(pressure, u, v, displacement);
    EXPECT_EQ(force.getRegion(), Region::Boundary);
    EXPECT_EQ(tangent.getRegion(), Region::Boundary);
    for (auto face = mesh.getBoundary(); face; ++face)
    {
      SCOPED_TRACE(face->getIndex());
      force.setPolytope(*face);
      tangent.setPolytope(*face);
      EXPECT_EQ(force.getPolytope().getIndex(), face->getIndex());
      EXPECT_EQ(tangent.getPolytope().getIndex(), face->getIndex());
      const auto& vertices = face->getVertices();
      const size_t count = vertices.size() * dimension;
      Math::Vector<Real> original(count);
      Math::Matrix<Real> jacobian(count, count);
      for (size_t test = 0; test < count; ++test)
      {
        original(test) = force.integrate(test);
        for (size_t trial = 0; trial < count; ++trial)
          jacobian(test, trial) = tangent.integrate(trial, test);
      }
      // Independent geometric oracle: transform the face vertices by the
      // prescribed affine displacement and form their oriented area vector.
      Eigen::Vector3d x[3], reference[3];
      for (size_t a = 0; a < vertices.size(); ++a)
      {
        reference[a].setZero();
        x[a].setZero();
        const auto& coordinates = mesh.getVertexCoordinates(vertices[a]);
        const auto shifted = coordinates + gradient * coordinates;
        for (size_t axis = 0; axis < dimension; ++axis)
        {
          reference[a](axis) = coordinates(axis);
          x[a](axis) = shifted(axis);
        }
      }
      Eigen::Vector3d area, referenceArea;
      if (dimension == 3)
      {
        area = (x[1] - x[0]).cross(x[2] - x[0]) / 2;
        referenceArea = (reference[1] - reference[0]).cross(reference[2] - reference[0]);
      }
      else
      {
        area << x[1](1) - x[0](1), x[0](0) - x[1](0), 0;
        referenceArea << reference[1](1) - reference[0](1),
          reference[0](0) - reference[1](0), 0;
      }
      const auto& adjacent =
        mesh.getConnectivity().getIncidence({dimension - 1, dimension}, face->getIndex());
      ASSERT_EQ(adjacent.size(), 1);
      const auto cell = mesh.getPolytope(dimension, adjacent[0]);
      const Point center(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      Eigen::Vector3d inward = -reference[0];
      for (size_t axis = 0; axis < dimension; ++axis)
        inward(axis) += center(axis);
      if (referenceArea.dot(inward) > 0)
        area *= -1;
      for (size_t a = 0; a < vertices.size(); ++a)
      {
        for (size_t axis = 0; axis < dimension; ++axis)
        {
          EXPECT_NEAR(original(a * dimension + axis),
            pressure * area(axis) / vertices.size(), 1e-12);
        }
      }
      // Perturb every independent nodal displacement, not only an affine
      // direction, to verify the entire generally nonsymmetric tangent.
      const Real step = 1e-6;
      for (size_t trial = 0; trial < count; ++trial)
      {
        const Index dof =
          space.getGlobalIndex({0, vertices[trial / dimension]}, trial % dimension);
        const Real saved = displacement[dof];
        displacement[dof] = saved + step;
        force.setPolytope(*face);
        Math::Vector<Real> plus(count);
        for (size_t test = 0; test < count; ++test)
          plus(test) = force.integrate(test);
        displacement[dof] = saved - step;
        force.setPolytope(*face);
        for (size_t test = 0; test < count; ++test)
        {
          EXPECT_NEAR(jacobian(test, trial),
            (plus(test) - force.integrate(test)) / (2 * step), 1e-8);
        }
        displacement[dof] = saved;
      }
      std::unique_ptr<decltype(force)> forceCopy(force.copy());
      std::unique_ptr<decltype(tangent)> tangentCopy(tangent.copy());
      pressure *= 2;
      forceCopy->setPolytope(*face);
      tangentCopy->setPolytope(*face);
      for (size_t test = 0; test < count; ++test)
      {
        EXPECT_NEAR(forceCopy->integrate(test), 2 * original(test), 1e-12);
        for (size_t trial = 0; trial < count; ++trial)
        {
          EXPECT_NEAR(
            tangentCopy->integrate(trial, test), 2 * jacobian(test, trial), 1e-12);
        }
      }
      pressure /= 2;
    }
  }
}
