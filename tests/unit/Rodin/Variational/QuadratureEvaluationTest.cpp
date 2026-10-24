/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Reference association in coefficient and field derivative evaluation.
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_QuadratureEvaluation, MatrixCallablePreservesPointContract)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = mesh.getCell(0);
    const auto& formula = QF::PolytopeQuadratureFormula::get(3, cell->getGeometry());
    const PolytopeQuadrature quadrature(*cell, formula);
    size_t points = 0;
    const auto coefficient = MatrixFunction(2, 2, [&](const auto& point) {
      ++points;
      auto value = Math::SpatialMatrix<Real>::Identity(2, 2);
      value(0, 0) = point.x();
      return value;
    });
    const auto expected = coefficient.getValue(quadrature.getPoint(0));
    const auto actual = coefficient.getValue(IntegrationPoint(quadrature, 0));
    EXPECT_EQ(points, 2u);
    EXPECT_EQ((actual - expected).norm(), Real(0));
    const auto invalid = MatrixFunction(
      2, 2, [](const auto&) { return Math::SpatialMatrix<Real>::Identity(1, 1); });
    EXPECT_THROW(invalid.getValue(IntegrationPoint(quadrature, 0)), Alert::Exception);
    struct PointOnly
    {
      Math::SpatialMatrix<Real> operator()(const Point&) const
      {
        return Math::SpatialMatrix<Real>::Identity(2, 2);
      }
      Math::SpatialMatrix<Real> operator()(const IntegrationPoint&) &
      {
        ADD_FAILURE() << "A const coefficient cannot invoke a nonconst sample overload";
        return Math::SpatialMatrix<Real>::Identity(2, 2);
      }
    };
    const auto pointOnly = MatrixFunction(2, 2, PointOnly{});
    EXPECT_EQ((pointOnly.getValue(IntegrationPoint(quadrature, 0)) -
                Math::SpatialMatrix<Real>::Identity(2, 2))
                .norm(),
      Real(0));
  }

  TEST(Rodin_Variational_QuadratureEvaluation, MassCoefficientPreservesAssociation)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    H1<2, Real> space(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction trial(space);
    TestFunction test(space);
    size_t samples = 0;
    const auto coefficient = RealFunction([&](const auto& point) -> Real {
      if constexpr (std::is_same_v<std::decay_t<decltype(point)>, IntegrationPoint>)
      {
        EXPECT_NE(point.getReferenceQuadratureFormula(), nullptr);
        ++samples;
      }
      else
        ADD_FAILURE() << "Mass coefficient lost quadrature association";
      return 2;
    });
    MassForm actual(coefficient, trial, test);
    MassForm expected(Real(2), trial, test);
    EXPECT_GT(samples, 0u);
    EXPECT_NEAR((actual.getOperator() - expected.getOperator()).norm(), 0, 1e-13);
  }

  template <class Space>
  void checkMatrixQuadratureDerivatives(const Space& space)
  {
    using Scalar = typename FormLanguage::Traits<Space>::ScalarType;
    GridFunction field(space);
    for (Index i = 0; i < static_cast<Index>(field.getSize()); ++i)
    {
      if constexpr (std::is_same_v<Scalar, Complex>)
        field.getData()[i] = Scalar(Real(i + 1) / 17, Real(i + 2) / 23);
      else
        field.getData()[i] = Real(i + 1) / 17;
    }
    const auto cell = space.getMesh().getCell(0);
    const auto& formula = QF::PolytopeQuadratureFormula::get(3, cell->getGeometry());
    const PolytopeQuadrature quadrature(*cell, formula);
    Grad gradient(field);
    TrialFunction trial(space);
    Grad shape(trial);
    shape.setPolytope(*cell);
    const auto& element = space.getFiniteElement(cell->getDimension(), cell->getIndex());
    const auto& shapeFormula = QF::PolytopeQuadratureFormula::get(4, cell->getGeometry());
    const PolytopeQuadrature shapeQuadrature(*cell, shapeFormula);
    const IntegrationPoint shapeSample(shapeQuadrature, 0);
    shape.setIntegrationPoint(shapeSample);
    EXPECT_TRUE(element.getScalarElement().findTabulation(shapeFormula).has_value());
    for (size_t q = 0; q < quadrature.getSize(); ++q)
    {
      const auto& point = quadrature.getPoint(q);
      const IntegrationPoint sample(quadrature, q);
      const auto expected = gradient.getValue(point);
      const auto actual = gradient.getValue(sample);
      for (size_t row = 0; row < space.getRows(); ++row)
      {
        for (size_t column = 0; column < space.getColumns(); ++column)
        {
          for (size_t axis = 0; axis < space.getMesh().getSpaceDimension(); ++axis)
          {
            EXPECT_NEAR(std::abs(actual(row, column, axis) - expected(row, column, axis)),
              0, 1e-10);
          }
        }
      }
      EXPECT_TRUE(element.getScalarElement().findTabulation(formula).has_value());
      // A changed reference coordinate cannot inherit the original sample index.
      auto coordinates = point.getReferenceCoordinates();
      coordinates *= Real(0.75);
      const Point shifted(*cell, coordinates);
      const IntegrationPoint provenance(shifted, &formula, q);
      const auto shiftedExpected = gradient.getValue(shifted);
      const auto shiftedActual = gradient.getValue(provenance);
      const IntegrationPoint unproven(point, &formula, q);
      shape.setIntegrationPoint(sample);
      for (size_t local = 0; local < element.getCount(); ++local)
      {
        const auto cached = shape.getBasis(local);
        shape.setIntegrationPoint(unproven);
        const auto direct = shape.getBasis(local);
        for (size_t row = 0; row < space.getRows(); ++row)
        {
          for (size_t column = 0; column < space.getColumns(); ++column)
          {
            for (size_t axis = 0; axis < space.getMesh().getSpaceDimension(); ++axis)
            {
              EXPECT_NEAR(std::abs(cached(row, column, axis) - direct(row, column, axis)),
                0, 1e-11);
              EXPECT_NEAR(std::abs(shiftedActual(row, column, axis) -
                            shiftedExpected(row, column, axis)),
                0, 1e-11);
            }
          }
        }
        shape.setIntegrationPoint(sample);
      }
    }
  }

  TEST(Rodin_Variational_QuadratureEvaluation, MatrixDerivativesUseReferenceTables)
  {
    const auto check = []<class Scalar>() {
      for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
             Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
             Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
      {
        SCOPED_TRACE(static_cast<int>(geometry));
        Array<size_t> shape(Polytope::Traits(geometry).getDimension());
        shape.setConstant(2);
        auto mesh = LocalMesh::UniformGrid(geometry, shape);
        for (size_t d = mesh.getDimension(); d > 0; --d)
          mesh.getConnectivity().compute(d, d - 1);
        H1<3, Math::SpatialMatrix<Scalar>> high(
          std::integral_constant<size_t, 3>{}, mesh, 2, 2);
        checkMatrixQuadratureDerivatives(high);
        P1<Math::SpatialMatrix<Scalar>> low(mesh, 2, 2);
        checkMatrixQuadratureDerivatives(low);
      }
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
  }

  TEST(Rodin_Variational_QuadratureEvaluation, P1FieldDerivativesUseReferenceTables)
  {
    const auto check = []<class Scalar>() {
      for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
             Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
             Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
      {
        Array<size_t> shape(Polytope::Traits(geometry).getDimension());
        shape.setConstant(2);
        auto mesh = LocalMesh::UniformGrid(geometry, shape);
        for (size_t d = mesh.getDimension(); d > 0; --d)
          mesh.getConnectivity().compute(d, d - 1);
        P1<Scalar> scalar(mesh);
        P1<Math::SpatialVector<Scalar>> vector(mesh, mesh.getSpaceDimension());
        GridFunction scalarField(scalar);
        GridFunction vectorField(vector);
        for (Index i = 0; i < static_cast<Index>(scalarField.getSize()); ++i)
        {
          if constexpr (std::is_same_v<Scalar, Complex>)
            scalarField.getData()[i] = Scalar(Real(i + 1) / 11, Real(i + 3) / 13);
          else
            scalarField.getData()[i] = Real(i + 1) / 11;
        }
        for (Index i = 0; i < static_cast<Index>(vectorField.getSize()); ++i)
        {
          if constexpr (std::is_same_v<Scalar, Complex>)
            vectorField.getData()[i] = Scalar(Real(i + 2) / 17, Real(i + 5) / 19);
          else
            vectorField.getData()[i] = Real(i + 2) / 17;
        }
        const auto cell = mesh.getCell(0);
        const auto& formula = QF::PolytopeQuadratureFormula::get(3, geometry);
        const PolytopeQuadrature quadrature(*cell, formula);
        for (size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& point = quadrature.getPoint(q);
          const IntegrationPoint sample(quadrature, q);
          EXPECT_NEAR(
            (Grad(scalarField).getValue(sample) - Grad(scalarField).getValue(point))
              .norm(),
            0, 1e-11);
          EXPECT_NEAR((Jacobian(vectorField).getValue(sample) -
                        Jacobian(vectorField).getValue(point))
                        .norm(),
            0, 1e-11);
          EXPECT_NEAR(std::abs(Div(vectorField).getValue(sample) -
                        Div(vectorField).getValue(point)),
            0, 1e-11);
          auto coordinates = point.getReferenceCoordinates();
          coordinates *= Real(0.75);
          const Point shifted(*cell, coordinates);
          const IntegrationPoint provenance(shifted, &formula, q);
          EXPECT_NEAR(
            (Grad(scalarField).getValue(provenance) - Grad(scalarField).getValue(shifted))
              .norm(),
            0, 1e-11);
          EXPECT_NEAR((Jacobian(vectorField).getValue(provenance) -
                        Jacobian(vectorField).getValue(shifted))
                        .norm(),
            0, 1e-11);
          EXPECT_NEAR(std::abs(Div(vectorField).getValue(provenance) -
                        Div(vectorField).getValue(shifted)),
            0, 1e-11);
          for (size_t axis = 0; axis < mesh.getSpaceDimension(); ++axis)
          {
            EXPECT_NEAR(std::abs(Derivative(axis, scalarField).getValue(sample) -
                          Derivative(axis, scalarField).getValue(point)),
              0, 1e-11);
            EXPECT_NEAR(std::abs(Derivative(axis, scalarField).getValue(provenance) -
                          Derivative(axis, scalarField).getValue(shifted)),
              0, 1e-11);
          }
        }
        EXPECT_TRUE(scalar.getFiniteElement(cell->getDimension(), cell->getIndex())
            .findTabulation(formula)
            .has_value());
      }
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
  }
  TEST(Rodin_Variational_QuadratureEvaluation, P1DivergenceRejectsInvalidVectorDimension)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = mesh.getCell(0);
    const auto& formula = QF::PolytopeQuadratureFormula::get(3, cell->getGeometry());
    const PolytopeQuadrature quadrature(*cell, formula);
    for (const size_t dimension : {size_t(1), size_t(3)})
    {
      P1<Math::SpatialVector<Real>> space(mesh, dimension);
      GridFunction field(space);
      field.getData().setOnes();
      const Div divergence(field);
      EXPECT_THROW(divergence.getValue(quadrature.getPoint(0)), Alert::Exception);
      EXPECT_THROW(
        divergence.getValue(IntegrationPoint(quadrature, 0)), Alert::Exception);
    }
  }

  TEST(Rodin_Variational_QuadratureEvaluation, ReferenceVectorBasisJacobians)
  {
    const auto check = []<class Scalar>() {
      for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
             Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
             Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
      {
        const auto& formula = QF::PolytopeQuadratureFormula::get(3, geometry);
        const auto compare = [&](const auto& element) {
          for (size_t q = 0; q < formula.getSize(); ++q)
          {
            for (size_t local = 0; local < element.getCount(); ++local)
            {
              const auto derivative = element.getBasis(local).getJacobian();
              EXPECT_NEAR(
                (derivative(formula, q) - derivative(formula.getPoint(q))).norm(), 0,
                1e-11);
            }
          }
        };
        compare(H1Element<2, Math::SpatialVector<Scalar>>(geometry, 3));
        compare(P1Element<Math::SpatialVector<Scalar>>(geometry, 3));
      }
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
  }

}
