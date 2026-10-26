/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Distinct trial/test spaces agree with independent generic quadrature.
 *
 * Real/complex scalar and vector P1/H1 mass and derivative kernels are checked
 * entrywise on every positive-dimensional geometry. Complex coefficients use
 * complex-valued spaces. Distinct identities exercise independent test-basis
 * caches; quadratic trials with linear tests exercise rectangular operators.
 */
#include <algorithm>
#include <complex>
#include <type_traits>

#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include "../../../QuadratureReference.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Compares every local matrix entry with independent quadrature evaluation.
   * @param trial Trial expression, including any coefficient.
   * @param test Test expression.
   * @param cell Physical integration cell.
   */
  template <class Trial, class Test>
  void checkDistinctRule(const Trial& trial, const Test& test, const Polytope& cell)
  {
    constexpr size_t order = 8;
    // Equivalent contractions and summation orders differ by roundoff.
    constexpr Real relativeTolerance = 1e-12;
    auto integral = Integral(trial, test);
    QuadratureReference<decltype(integral)> reference(Dot(trial, test));
    reference.assemble(cell, order);
    const auto& expected = reference.getOperator();
    ASSERT_GT(expected.norm(), 0);
    integral.setOrder(order);
    integral.setPolytope(cell);
    EXPECT_EQ(&integral.getPolytope(), &cell);
    for (Index row = 0; row < static_cast<Index>(expected.rows()); ++row)
    {
      for (Index column = 0; column < static_cast<Index>(expected.cols()); ++column)
      {
        const auto actual = integral.integrate(column, row);
        const Real tolerance =
          relativeTolerance * std::max(Real(1), std::abs(expected(row, column)));
        EXPECT_LE(std::abs(actual - expected(row, column)), tolerance)
          << "row " << row << ", column " << column;
      }
    }
  }

  /**
   * @brief Checks mass and derivative forms with compatible scalar coefficients.
   * @param trialSpace Independently constructed trial finite element space.
   * @param testSpace Independently constructed test finite element space.
   * @param cell Physical integration cell shared by both spaces.
   */
  template <class TrialFES, class TestFES>
  void checkDistinctSpaces(
    const TrialFES& trialSpace, const TestFES& testSpace, const Polytope& cell)
  {
    ASSERT_NE(
      static_cast<const void*>(&trialSpace), static_cast<const void*>(&testSpace));
    TrialFunction u(trialSpace);
    TestFunction v(testSpace);
    const auto check = [&](const auto& trial, const auto& test) {
      checkDistinctRule(trial, test, cell);
      checkDistinctRule(0.25 * trial, test, cell);
      if constexpr (std::is_same_v<typename FormLanguage::Traits<TrialFES>::ScalarType,
                      Complex>)
        checkDistinctRule(Complex(0.25, 0.5) * trial, test, cell);
    };
    check(u, v);
    if constexpr (FormLanguage::IsVectorRange<
                    typename FormLanguage::Traits<TrialFES>::RangeType>::Value)
      check(Jacobian(u), Jacobian(v));
    else
      check(Grad(u), Grad(v));
  }

  /// @brief Exercises trial/test identities and polynomial orders on each cell geometry.
  class DistinctSpaceQuadrature : public ::testing::TestWithParam<Polytope::Type>
  {};

  /// @brief Distinct P1 instances retain every scalar/vector trial-test entry.
  TEST_P(DistinctSpaceQuadrature, P1RealComplexScalarVector)
  {
    const size_t dimension = Polytope::Traits(GetParam()).getDimension();
    auto mesh = LocalMesh::UniformGrid(GetParam(), Array<size_t>::Constant(dimension, 2));
    const auto cell = mesh.getCell();
    const auto check = [&]<class Scalar>() {
      const P1<Scalar> scalarTrial(mesh), scalarTest(mesh);
      const P1<Math::SpatialVector<Scalar>> vectorTrial(mesh, dimension),
        vectorTest(mesh, dimension);
      checkDistinctSpaces(scalarTrial, scalarTest, *cell);
      checkDistinctSpaces(vectorTrial, vectorTest, *cell);
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
  }

  /// @brief Distinct quadratic H1 instances retain independent derivative caches.
  TEST_P(DistinctSpaceQuadrature, H1RealComplexScalarVector)
  {
    const size_t dimension = Polytope::Traits(GetParam()).getDimension();
    auto mesh = LocalMesh::UniformGrid(GetParam(), Array<size_t>::Constant(dimension, 2));
    for (size_t d = 1; d <= dimension; ++d)
    {
      mesh.getConnectivity().compute(d, 0);
      mesh.getConnectivity().compute(dimension, d);
      for (size_t lower = 1; lower < d; ++lower)
        mesh.getConnectivity().compute(d, lower);
    }
    const auto cell = mesh.getCell();
    const auto check = [&]<class Scalar>() {
      const H1<2, Scalar> scalarTrial(std::integral_constant<size_t, 2>{}, mesh),
        scalarTest(std::integral_constant<size_t, 2>{}, mesh);
      const H1<2, Math::SpatialVector<Scalar>> vectorTrial(
        std::integral_constant<size_t, 2>{}, mesh, dimension),
        vectorTest(std::integral_constant<size_t, 2>{}, mesh, dimension);
      checkDistinctSpaces(scalarTrial, scalarTest, *cell);
      checkDistinctSpaces(vectorTrial, vectorTest, *cell);
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
  }

  /// @brief Quadratic trials and linear tests produce the full rectangular H1 operator.
  TEST_P(DistinctSpaceQuadrature, H1MixedOrdersRealComplexScalarVector)
  {
    const size_t dimension = Polytope::Traits(GetParam()).getDimension();
    auto mesh = LocalMesh::UniformGrid(GetParam(), Array<size_t>::Constant(dimension, 2));
    for (size_t d = 1; d <= dimension; ++d)
    {
      mesh.getConnectivity().compute(d, 0);
      mesh.getConnectivity().compute(dimension, d);
      for (size_t lower = 1; lower < d; ++lower)
        mesh.getConnectivity().compute(d, lower);
    }
    const auto cell = mesh.getCell();
    const auto check = [&]<class Scalar>() {
      const H1<2, Scalar> scalarTrial(std::integral_constant<size_t, 2>{}, mesh);
      const H1<1, Scalar> scalarTest(std::integral_constant<size_t, 1>{}, mesh);
      const H1<2, Math::SpatialVector<Scalar>> vectorTrial(
        std::integral_constant<size_t, 2>{}, mesh, dimension);
      const H1<1, Math::SpatialVector<Scalar>> vectorTest(
        std::integral_constant<size_t, 1>{}, mesh, dimension);
      checkDistinctSpaces(scalarTrial, scalarTest, *cell);
      checkDistinctSpaces(vectorTrial, vectorTest, *cell);
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, DistinctSpaceQuadrature,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge));
}
