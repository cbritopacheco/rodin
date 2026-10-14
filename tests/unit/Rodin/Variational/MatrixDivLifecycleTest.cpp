/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Matrix divergence ownership and dimension contracts. */
#include <gtest/gtest.h>
#include <memory>
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
  void checkMatrixDivergence(const FES& space, bool constant)
  {
    using Scalar = typename FormLanguage::Traits<FES>::ScalarType;
    const size_t dimension = space.getMesh().getSpaceDimension();
    const Scalar factor = complexFactor<Scalar>();
    GridFunction field(space);
    const size_t rows = space.getRows(), columns = space.getColumns();
    field = MatrixFunction(rows, columns, [=](const Geometry::Point& point) {
      Math::SpatialMatrix<Scalar> result(rows, columns);
      for (size_t row = 0; row < rows; ++row)
      {
        for (size_t col = 0; col < columns; ++col)
        {
          result(row, col) = factor * Real(1 + row + col);
          if (!constant)
            for (size_t k = 0; k < dimension; ++k)
            {
              result(row, col) += factor * Real((row + 1) * (col + 2) * (k + 1)) *
                point.getPhysicalCoordinates()[k];
            }
        }
      }
      return result;
    });
    using Expression = decltype(Div(field));
    std::unique_ptr<Expression> owned;
    {
      Expression original(field);
      Expression copied(original);
      Expression moved(std::move(copied));
      owned.reset(moved.copy());
    }
    EXPECT_EQ(&owned->getOperand(), &field);
    TrialFunction trial(space);
    using ShapeExpression = decltype(Div(trial));
    std::unique_ptr<ShapeExpression> shape;
    {
      ShapeExpression original(trial);
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
      Math::SpatialVector<Scalar> reconstructed(space.getRows());
      reconstructed.setZero();
      const auto& dofs = space.getDOFs(cell->getDimension(), cell->getIndex());
      ASSERT_EQ(shape->getDOFs(*cell), dofs.size());
      for (size_t local = 0; local < static_cast<size_t>(dofs.size()); ++local)
        reconstructed += field[dofs[local]] * shape->getBasis(local);
      for (size_t row = 0; row < rows; ++row)
      {
        Real derivativeSum = 0;
        if (!constant)
          for (size_t k = 0; k < dimension; ++k)
            derivativeSum += Real((row + 1) * (k + 2) * (k + 1));
        const Scalar expected = factor * derivativeSum;
        EXPECT_NEAR(Math::abs(evaluated[row] - expected), 0, 1e-11);
        EXPECT_NEAR(Math::abs(reconstructed[row] - expected), 0, 1e-11);
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
      {
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      }
      checkMatrixDivergence(P0<Range>(mesh, 2, dimension), true);
      checkMatrixDivergence(P0g<Range, LocalMesh>(mesh, 2, dimension), true);
      checkMatrixDivergence(P1<Range>(mesh, 2, dimension), false);
      checkMatrixDivergence(
        H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh, 2, dimension), false);
    }
  }

  template <class Action>
  void expectDimensionDiagnostic(Action action, const std::string& member)
  {
    try
    {
      action();
      FAIL() << "Incompatible matrix column count was accepted";
    }
    catch (const Alert::Exception& exception)
    {
      const std::string message = exception.what();
      EXPECT_NE(
        message.find("In member function " + member + " of class "), std::string::npos);
      EXPECT_NE(message.find("Div<"), std::string::npos);
      EXPECT_NE(message.find(
                  "Matrix divergence requires columns equal to the spatial dimension."),
        std::string::npos);
    }
  }

  template <class Scalar>
  void checkInvalidColumns()
  {
    using Range = Math::SpatialMatrix<Scalar>;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    // A three-column matrix field is valid; its divergence in two dimensions is undefined.
    P1<Range> space(mesh, 2, 3);
    GridFunction field(space);
    Range constant(2, 3);
    constant.setConstant(Scalar(1));
    field = constant;
    TrialFunction trial(space);
    auto divergence = Div(field);
    auto shape = Div(trial);
    const auto cell = mesh.getCell();
    const Geometry::Point point(
      *cell, Polytope::Traits(cell->getGeometry()).getCentroid());
    const IntegrationPoint ip(point);
    shape.setIntegrationPoint(ip);
    expectDimensionDiagnostic([&] { (void)divergence.getValue(ip); }, "getValue");
    expectDimensionDiagnostic([&] { (void)shape.getBasis(0); }, "getBasis");
    EXPECT_EQ(field.getRows(), 2);
    EXPECT_EQ(field.getColumns(), 3);
    EXPECT_EQ((field.getValue(ip) - constant).norm(), 0);
  }
}

/// @brief Real and complex copies retain analytic matrix divergence and basis reconstruction.
TEST(MatrixDivLifecycle, RealFieldsAndClones)
{
  checkSpaces<Real>();
}
TEST(MatrixDivLifecycle, ComplexFieldsAndClones)
{
  checkSpaces<Complex>();
}

/// @brief Undefined matrix divergences report the evaluating member and preserve the field.
TEST(MatrixDivLifecycle, RealDimensionDiagnostics)
{
  checkInvalidColumns<Real>();
}
TEST(MatrixDivLifecycle, ComplexDimensionDiagnostics)
{
  checkInvalidColumns<Complex>();
}
