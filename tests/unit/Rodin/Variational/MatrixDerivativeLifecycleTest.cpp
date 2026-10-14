/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Matrix derivative and Jacobian ownership and direction contracts. */
#include <gtest/gtest.h>
#include <memory>
#include <limits>
#include <string>
#include <type_traits>
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <class Scalar>
  Scalar complexFactor()
  {
    if constexpr (std::is_same_v<Scalar, Complex>)
      return Scalar(1, -0.5);
    else
      return Scalar(1);
  }

  template <class FES>
  void checkMatrixDerivatives(const FES& space, bool constant)
  {
    using Scalar = typename FormLanguage::Traits<FES>::ScalarType;
    const size_t dimension = space.getMesh().getSpaceDimension();
    const Scalar factor = complexFactor<Scalar>();
    GridFunction field(space);
    const size_t rows = space.getRows(), columns = space.getColumns();
    field = MatrixFunction(rows, columns, [=](const Geometry::Point& point) {
      Math::SpatialMatrix<Scalar> result(rows, columns);
      for (size_t row = 0; row < rows; ++row)
        for (size_t col = 0; col < columns; ++col)
        {
          result(row, col) = factor * Real(1 + row + col);
          if (!constant)
            for (size_t k = 0; k < dimension; ++k)
              result(row, col) += factor * Real((row + 1) * (col + 2) * (k + 1)) *
                point.getPhysicalCoordinates()[k];
        }
      return result;
    });
    using Expression = decltype(Jacobian(field));
    std::unique_ptr<Expression> owned;
    {
      Expression original(field);
      Expression copied(original);
      Expression moved(std::move(copied));
      owned.reset(moved.copy());
    }
    EXPECT_EQ(&owned->getOperand(), &field);
    TrialFunction trial(space);
    using ShapeExpression = decltype(Derivative(size_t(0), trial));
    std::unique_ptr<ShapeExpression> shape;
    {
      ShapeExpression original(size_t(0), trial);
      ShapeExpression copied(original);
      ShapeExpression moved(std::move(copied));
      shape.reset(moved.copy());
    }
    EXPECT_EQ(shape->getOperand().getUUID(), trial.getUUID());
    EXPECT_EQ(&shape->getOperand().getFiniteElementSpace(), &space);
    for (auto cell = space.getMesh().getCell(); cell; ++cell)
    {
      const Geometry::Point point(
        *cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      const IntegrationPoint ip(point);
      shape->setIntegrationPoint(ip);
      const auto evaluated = owned->getValue(ip);
      EXPECT_EQ(owned->getOrder(*cell), Grad(field).getOrder(*cell));
      for (size_t direction = 0; direction < dimension; ++direction)
      {
        using Partial = decltype(Derivative(direction, field));
        std::unique_ptr<Partial> partial;
        {
          Partial original(direction, field);
          Partial copied(original);
          Partial moved(std::move(copied));
          partial.reset(moved.copy());
        }
        EXPECT_EQ(&partial->getOperand(), &field);
        EXPECT_EQ(partial->getOrder(*cell), Grad(field).getOrder(*cell));
        const auto derivative = partial->getValue(ip);
        std::unique_ptr<ShapeExpression> basis;
        {
          ShapeExpression expression(direction, trial);
          ShapeExpression copied(expression);
          ShapeExpression moved(std::move(copied));
          basis.reset(moved.copy());
        }
        basis->setIntegrationPoint(ip);
        Math::SpatialMatrix<Scalar> reconstructed(rows, columns);
        reconstructed.setZero();
        const auto& dofs = space.getDOFs(cell->getDimension(), cell->getIndex());
        ASSERT_EQ(basis->getDOFs(*cell), static_cast<size_t>(dofs.size()));
        for (size_t local = 0; local < static_cast<size_t>(dofs.size()); ++local)
          reconstructed += field[dofs[local]] * basis->getBasis(local);
        for (size_t row = 0; row < rows; ++row)
          for (size_t col = 0; col < columns; ++col)
          {
            const Scalar expected = constant
              ? Scalar(0)
              : factor * Real((row + 1) * (col + 2) * (direction + 1));
            EXPECT_NEAR(Math::abs(evaluated(row, col, direction) - expected), 0, 1e-11);
            EXPECT_NEAR(Math::abs(derivative(row, col) - expected), 0, 1e-11);
            EXPECT_NEAR(Math::abs(reconstructed(row, col) - expected), 0, 1e-11);
          }
      }
    }
  }

  template <class Scalar>
  void checkSpaces()
  {
    using Range = Math::SpatialMatrix<Scalar>;
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(int(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {2})
        : dimension == 2         ? LocalMesh::UniformGrid(geometry, {2, 2})
                                 : LocalMesh::UniformGrid(geometry, {2, 2, 2});
      for (size_t d = 1; d <= dimension; ++d)
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      checkMatrixDerivatives(P0<Range>(mesh, 2, dimension), true);
      checkMatrixDerivatives(P0g<Range, LocalMesh>(mesh, 2, dimension), true);
      checkMatrixDerivatives(P1<Range>(mesh, 2, dimension), false);
      checkMatrixDerivatives(
        H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh, 2, dimension), false);
    }
  }

  template <class Action>
  void expectDimensionDiagnostic(Action action, const std::string& member)
  {
    try
    {
      action();
      FAIL() << "Out-of-range derivative direction was accepted";
    }
    catch (const Alert::Exception& exception)
    {
      const std::string message = exception.what();
      EXPECT_NE(
        message.find("In member function " + member + " of class "), std::string::npos);
      EXPECT_NE(message.find("Derivative<"), std::string::npos);
      EXPECT_NE(
        message.find("Partial derivative direction exceeds the spatial dimension."),
        std::string::npos);
    }
  }

  template <class Scalar>
  void checkInvalidDirections()
  {
    using Range = Math::SpatialMatrix<Scalar>;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    // Both the first excluded direction and an arbitrarily large direction must reject evaluation.
    P1<Range> space(mesh, 2, 3);
    GridFunction field(space);
    Range constant(2, 3);
    constant.setConstant(Scalar(1));
    field = constant;
    TrialFunction trial(space);

    const auto cell = mesh.getCell();
    const Geometry::Point point(
      *cell, Polytope::Traits(cell->getGeometry()).getCentroid());
    const IntegrationPoint ip(point);
    for (const size_t direction : {size_t(2), std::numeric_limits<size_t>::max()})
    {
      auto partial = Derivative(direction, field);
      auto shape = Derivative(direction, trial);
      shape.setIntegrationPoint(ip);
      expectDimensionDiagnostic([&] { (void)partial.getValue(ip); }, "getValue");
      expectDimensionDiagnostic([&] { (void)shape.getBasis(0); }, "getBasis");
    }
    EXPECT_EQ(field.getRows(), 2);
    EXPECT_EQ(field.getColumns(), 3);
    EXPECT_EQ((field.getValue(ip) - constant).norm(), 0);
  }
}

/// @brief Real and complex copies retain analytic matrix Jacobians and partial derivatives.
TEST(MatrixDerivativeLifecycle, RealFieldsAndClones)
{
  checkSpaces<Real>();
}
TEST(MatrixDerivativeLifecycle, ComplexFieldsAndClones)
{
  checkSpaces<Complex>();
}

/// @brief Out-of-range directions report the evaluating member and preserve the field.
TEST(MatrixDerivativeLifecycle, RealDirectionDiagnostics)
{
  checkInvalidDirections<Real>();
}
TEST(MatrixDerivativeLifecycle, ComplexDirectionDiagnostics)
{
  checkInvalidDirections<Complex>();
}
