/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "../../../QuadratureReference.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief Reference Jacobians are evaluated at each nonsimplex quadrature sample.
  TEST(RankOneJacobianQuadrature, NonsimplexGradients)
  {
    for (const auto geometry : {Polytope::Type::Quadrilateral, Polytope::Type::Hexahedron,
           Polytope::Type::Wedge, Polytope::Type::Pyramid})
    {
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
      P1<Math::SpatialVector<Real>> space(mesh, 1);
      TrialFunction u(space);
      TestFunction v(space);
      Math::SpatialMatrix<Real> matrix(1, dimension);
      for (size_t axis = 0; axis < dimension; ++axis)
        matrix(0, axis) = Real(axis + 1);
      MatrixFunction coefficient(matrix);
      auto integral =
        Integral(Dot(coefficient, Jacobian(u)), Dot(coefficient, Jacobian(v)));
      integral.setOrder(6);
      const auto cell = mesh.getCell();
      QuadratureReference<decltype(integral)> reference(integral);
      reference.assemble(*cell, 6);
      integral.setPolytope(*cell);
      for (size_t test = 0; test < v.getDOFs(*cell); ++test)
      {
        for (size_t trial = 0; trial < u.getDOFs(*cell); ++trial)
        {
          EXPECT_NEAR(
            integral.integrate(trial, test), reference.getOperator()(test, trial), 1e-11);
        }
      }
    }
  }

  /// @brief An explicitly polynomial matrix coefficient used to check inferred order.
  class AffineRankOneCoefficient final
    : public MatrixFunctionBase<Real, AffineRankOneCoefficient>
  {
    public:
      size_t getRows() const
      {
        return 1;
      }
      size_t getColumns() const
      {
        return 2;
      }
      Math::SpatialMatrix<Real> getValue(const Point& point) const
      {
        Math::SpatialMatrix<Real> value(1, 2);
        value(0, 0) = point.x();
        value(0, 1) = 0;
        return value;
      }
      Optional<size_t> getOrder([[maybe_unused]] const Polytope& cell) const
      {
        return 1;
      }
      AffineRankOneCoefficient* copy() const noexcept override
      {
        return new AffineRankOneCoefficient(*this);
      }
  };

  /// @brief Both coefficient degrees contribute to default P1 and H1 quadrature.
  TEST(RankOneJacobianQuadrature, CoefficientOrder)
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(2)
                  .nodes(3)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({0, 1})
                  .polytope(Polytope::Type::Triangle, {0, 1, 2})
                  .finalize();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    const auto cell = mesh.getCell();
    const auto check = [&](const auto& space) {
      TrialFunction u(space);
      TestFunction v(space);
      AffineRankOneCoefficient coefficient;
      auto integral =
        Integral(Dot(coefficient, Jacobian(u)), Dot(coefficient, Jacobian(v)));
      integral.setPolytope(*cell);
      EXPECT_NEAR(integral.integrate(0, 0), Real(1) / 12, 1e-13);
    };
    P1<Math::SpatialVector<Real>> linear(mesh, 1);
    H1<1, Math::SpatialVector<Real>> higher(std::integral_constant<size_t, 1>{}, mesh, 1);
    check(linear);
    check(higher);
  }

  /// @brief Shared and distinct coefficients retain high-order inference and overrides.
  TEST(RankOneJacobianQuadrature, HigherOrderCoefficients)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    H1<2, Math::SpatialVector<Real>> space(std::integral_constant<size_t, 2>{}, mesh, 1);
    TrialFunction u(space);
    TestFunction v(space);
    const auto cell = mesh.getCell();
    AffineRankOneCoefficient first, second;
    for (bool shared : {false, true})
    {
      const auto& rhs = shared ? first : second;
      const auto trial = Dot(first, Jacobian(u));
      const auto test = Dot(rhs, Jacobian(v));
      auto integral = Integral(trial, test);
      QuadratureReference<decltype(integral)> reference(Dot(trial, test));
      reference.assemble(*cell, 6);
      for (bool explicitOrder : {false, true})
      {
        if (explicitOrder)
          integral.setOrder(6);
        integral.setPolytope(*cell);
        for (size_t test = 0; test < v.getDOFs(*cell); ++test)
        {
          for (size_t trial = 0; trial < u.getDOFs(*cell); ++trial)
          {
            EXPECT_NEAR(integral.integrate(trial, test),
              reference.getOperator()(test, trial), 1e-12);
          }
        }
      }
    }
  }

  /// @brief Separate contractions permit different trial and test component counts.
  TEST(RankOneJacobianQuadrature, RectangularActions)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    H1<2, Math::SpatialVector<Real>> trialSpace(
      std::integral_constant<size_t, 2>{}, mesh, 3);
    H1<2, Math::SpatialVector<Real>> testSpace(
      std::integral_constant<size_t, 2>{}, mesh, 2);
    TrialFunction u(trialSpace);
    TestFunction v(testSpace);
    Math::SpatialMatrix<Real> first(3, 2), second(2, 2);
    for (size_t row = 0; row < first.rows(); ++row)
    {
      for (size_t column = 0; column < first.cols(); ++column)
        first(row, column) = 1 + row + 2 * column;
    }
    for (size_t row = 0; row < second.rows(); ++row)
    {
      for (size_t column = 0; column < second.cols(); ++column)
        second(row, column) = 2 + 2 * row + column;
    }
    const auto trial = Dot(MatrixFunction(first), Jacobian(u));
    const auto test = Dot(MatrixFunction(second), Jacobian(v));
    auto integral = Integral(trial, test);
    QuadratureReference<decltype(integral)> reference(Dot(trial, test));
    const auto cell = mesh.getCell();
    reference.assemble(*cell, 6);
    integral.setPolytope(*cell);
    for (size_t row = 0; row < v.getDOFs(*cell); ++row)
    {
      for (size_t column = 0; column < u.getDOFs(*cell); ++column)
      {
        EXPECT_NEAR(
          integral.integrate(column, row), reference.getOperator()(row, column), 1e-12);
      }
    }
  }
}
