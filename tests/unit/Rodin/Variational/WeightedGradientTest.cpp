/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <functional>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "../../../QuadratureReference.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief Complex nonsymmetric tensors act linearly on the trial gradient.
  TEST(WeightedGradient, ComplexTensor)
  {
    for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral})
    {
      auto mesh = LocalMesh::UniformGrid(geometry, {2, 2});
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1<2, Complex> trialSpace(std::integral_constant<size_t, 2>{}, mesh);
      H1<1, Complex> testSpace(std::integral_constant<size_t, 1>{}, mesh);
      TrialFunction u(trialSpace);
      TestFunction v(testSpace);
      Math::SpatialMatrix<Complex> tensor(2, 2);
      tensor(0, 0) = Complex(0, 1);
      tensor(0, 1) = Complex(2, 3);
      tensor(1, 0) = Complex(-1, 2);
      tensor(1, 1) = Complex(1, -1);
      MatrixFunction coefficient(tensor);
      auto integral = Integral(coefficient * Grad(u), Grad(v));
      integral.setOrder(5);
      const auto cell = mesh.getCell();
      QuadratureReference<decltype(integral)> reference(integral);
      reference.assemble(*cell, 5);
      integral.setPolytope(*cell);
      for (size_t test = 0; test < v.getDOFs(*cell); ++test)
      {
        for (size_t trial = 0; trial < u.getDOFs(*cell); ++trial)
        {
          EXPECT_NEAR(std::abs(integral.integrate(trial, test) -
                        reference.getOperator()(test, trial)),
            0, 1e-12);
        }
      }
    }
  }

  /// @brief Nested independent coefficients cannot change outer scratch dimensions.
  TEST(WeightedGradient, NestedCoefficientEvaluation)
  {
    auto segment = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
    auto hex = LocalMesh::UniformGrid(Polytope::Type::Hexahedron, {2, 2, 2});
    for (auto* mesh : {&segment, &hex})
    {
      for (size_t dimension = 1; dimension <= mesh->getDimension(); ++dimension)
      {
        for (size_t lower = 0; lower < dimension; ++lower)
          mesh->getConnectivity().compute(dimension, lower);
      }
    }
    for (const bool reverse : {false, true})
    {
      const auto& outerMesh = reverse ? hex : segment;
      const auto& innerMesh = reverse ? segment : hex;
      H1<1, Real> outerSpace(std::integral_constant<size_t, 1>{}, outerMesh);
      H1<1, Real> innerSpace(std::integral_constant<size_t, 1>{}, innerMesh);
      TrialFunction ou(outerSpace), iu(innerSpace);
      TestFunction ov(outerSpace), iv(innerSpace);
      using Callable = std::function<Real(const Point&)>;
      Callable constant = [](const Point&) { return Real(2); };
      auto inner = Integral(RealFunction(constant) * Grad(iu), Grad(iv));
      inner.setOrder(2);
      const auto innerCell = innerMesh.getCell();
      size_t calls = 0;
      Callable nested = [&](const Point&) {
        ++calls;
        inner.setPolytope(*innerCell);
        return Real(2);
      };
      auto actual = Integral(RealFunction(nested) * Grad(ou), Grad(ov));
      auto expected = Integral(RealFunction(constant) * Grad(ou), Grad(ov));
      static_assert(std::is_same_v<decltype(inner), decltype(actual)>);
      const auto cell = outerMesh.getCell();
      actual.setOrder(2);
      expected.setOrder(2);
      expected.setPolytope(*cell);
      actual.setPolytope(*cell);
      ASSERT_GT(calls, 0u);
      for (size_t test = 0; test < ov.getDOFs(*cell); ++test)
      {
        for (size_t trial = 0; trial < ou.getDOFs(*cell); ++trial)
          EXPECT_EQ(actual.integrate(trial, test), expected.integrate(trial, test));
      }
    }
  }
}

/// @brief Real tensor coefficients promote naturally on complex H1 carriers.
TEST(WeightedGradient, RealTensorOnComplexSpace)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  H1<1, Complex> scalar(std::integral_constant<size_t, 1>{}, mesh);
  H1<1, Math::SpatialVector<Complex>> vector(
    std::integral_constant<size_t, 1>{}, mesh, 2);
  TrialFunction u(scalar);
  TrialFunction w(vector);
  TestFunction v(scalar);
  TestFunction z(vector);
  Math::SpatialMatrix<Real> tensor(2, 2);
  tensor(0, 0) = 2;
  tensor(0, 1) = 3;
  tensor(1, 0) = -1;
  tensor(1, 1) = 4;
  MatrixFunction coefficient(tensor);
  auto gradient = Integral(coefficient * Grad(u), Grad(v));
  auto jacobian = Integral(coefficient * Jacobian(w), Jacobian(z));
  auto mass = Integral(coefficient * w, z);
  const auto cell = mesh.getCell();
  const auto check = [&](auto& integral, const auto& trial, const auto& test) {
    using IntegralType = std::remove_cvref_t<decltype(integral)>;
    Rodin::Tests::QuadratureReference<IntegralType> reference(integral);
    reference.assemble(*cell, 5);
    integral.setOrder(5);
    integral.setPolytope(*cell);
    for (size_t te = 0; te < test.getDOFs(*cell); ++te)
    {
      for (size_t tr = 0; tr < trial.getDOFs(*cell); ++tr)
      {
        EXPECT_LT(
          std::abs(integral.integrate(tr, te) - reference.getOperator()(te, tr)), 1e-12);
      }
    }
  };
  check(gradient, u, v);
  check(jacobian, w, z);
  check(mass, w, z);
}

/// @brief Matrix weights map different vector and matrix component counts in mass forms.
TEST(WeightedGradient, RectangularMassCarriers)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  const auto cell = mesh.getCell();
  const auto check = [&]<bool Matrix>() {
    using Range = std::conditional_t<Matrix, Math::SpatialMatrix<Complex>,
      Math::SpatialVector<Complex>>;
    auto trialSpace = [&] {
      if constexpr (Matrix)
        return H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh, 2, 3);
      else
        return H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh, 3);
    }();
    auto testSpace = [&] {
      if constexpr (Matrix)
        return H1<1, Range>(std::integral_constant<size_t, 1>{}, mesh, 1, 3);
      else
        return H1<1, Range>(std::integral_constant<size_t, 1>{}, mesh, 2);
    }();
    TrialFunction u(trialSpace);
    TestFunction v(testSpace);
    Math::SpatialMatrix<Real> tensor(Matrix ? 1 : 2, Matrix ? 2 : 3);
    for (size_t row = 0; row < tensor.rows(); ++row)
    {
      for (size_t col = 0; col < tensor.cols(); ++col)
        tensor(row, col) = 1 + row + 2 * col;
    }
    auto integral = Integral(MatrixFunction(tensor) * u, v);
    Rodin::Tests::QuadratureReference<decltype(integral)> reference(integral);
    reference.assemble(*cell, 5);
    integral.setOrder(5);
    integral.setPolytope(*cell);
    for (size_t te = 0; te < v.getDOFs(*cell); ++te)
    {
      for (size_t tr = 0; tr < u.getDOFs(*cell); ++tr)
      {
        EXPECT_LT(
          std::abs(integral.integrate(tr, te) - reference.getOperator()(te, tr)), 1e-12);
      }
    }
    auto scalarWeight = Integral(RealFunction(2) * u, v);
    EXPECT_THROW(scalarWeight.setPolytope(*cell), Alert::Exception);
    auto wrongShape =
      Integral(MatrixFunction(Math::SpatialMatrix<Real>::Identity(3, 3)) * u, v);
    EXPECT_THROW(wrongShape.setPolytope(*cell), Alert::Exception);
  };
  check.template operator()<false>();
  check.template operator()<true>();
}
