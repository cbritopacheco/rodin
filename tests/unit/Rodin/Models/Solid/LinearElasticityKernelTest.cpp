/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Exact local elasticity weighting and kernel-state regressions.
 */
#include <gtest/gtest.h>
#include "Rodin/Solid.h"
#include "Rodin/Variational.h"
#include "../../../../convergence/CurvedGeometry.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    /**
     * @brief Evaluates the original pairwise weighting expression independently.
     *
     * Reference derivatives are evaluated directly from the element, without
     * the production kernel's retained derivative or weighted-trial buffers.
     * Coefficient multiplication remains inside the pair loop and quadrature
     * accumulation order is unchanged. This is an exact implemented-rule
     * oracle, not a polynomial-exactness claim for curved geometry.
     *
     * @param space Vector-valued finite-element space.
     * @param cell Entity on which the local operator is assembled.
     * @param lambda First Lamé parameter.
     * @param mu Shear modulus.
     * @returns Original local matrix, with test rows and trial columns.
     */
    auto pairwiseOperator(
      const auto& space, const Polytope& cell, const auto& lambda, const auto& mu)
    {
      using Scalar =
        typename FormLanguage::Traits<std::decay_t<decltype(space)>>::ScalarType;
      const size_t dimension = cell.getDimension();
      const auto& element = space.getFiniteElement(dimension, cell.getIndex());
      const size_t count = element.getCount();
      const size_t order = 2 * element.getOrder() +
        std::max(lambda.getOrder(cell).value_or(size_t(0)),
          mu.getOrder(cell).value_or(size_t(0)));
      const auto& formula = QF::PolytopeQuadratureFormula::get(order, cell.getGeometry());
      const PolytopeQuadrature quadrature(cell, formula);
      Math::Matrix<Scalar> result(count, count);
      result.setZero();
      std::vector<Math::SpatialMatrix<Scalar>> symmetric(count);
      std::vector<Scalar> divergence(count);
      for (size_t q = 0; q < formula.getSize(); ++q)
      {
        const auto& point = quadrature.getPoint(q);
        const Scalar weight = Scalar(formula.getWeight(q) * point.getDistortion());
        const auto first = lambda.getValue(point);
        const auto shear = mu.getValue(point);
        for (size_t a = 0; a < count; ++a)
        {
          Math::SpatialMatrix<Scalar> reference(dimension, dimension);
          for (size_t row = 0; row < dimension; ++row)
          {
            for (size_t col = 0; col < dimension; ++col)
            {
              reference(row, col) = element.getBasis(a).template getDerivative<1>(
                row, col)(formula.getPoint(q));
            }
          }
          const Math::SpatialMatrix<Scalar> jacobian =
            reference * point.getJacobianInverse();
          symmetric[a] = jacobian + jacobian.adjoint();
          divergence[a] = jacobian.trace();
        }
        for (size_t te = 0; te < count; ++te)
        {
          for (size_t tr = 0; tr < count; ++tr)
          {
            result(te, tr) += weight *
              (Math::dot(first * divergence[tr], divergence[te]) +
                Scalar(0.5) * Math::dot(shear * symmetric[tr], symmetric[te]));
          }
        }
      }
      return result;
    }

    /**
     * @brief Certifies varying real/complex weights, cell changes and kernel copies.
     * @tparam RealCoefficients Whether coefficients stay real for complex fields.
     * @tparam DistinctSpaces Whether test functions borrow a separate copied space.
     * @param space Vector space on native or prescribed quadratic geometry.
     */
    template <bool RealCoefficients = false, bool DistinctSpaces = false>
    void checkWeighting(auto& space)
    {
      using Scalar =
        typename FormLanguage::Traits<std::decay_t<decltype(space)>>::ScalarType;
      using CoefficientScalar = std::conditional_t<RealCoefficients, Real, Scalar>;
      auto& mesh = space.getMesh();
      Real scale = 1;
      size_t lambdaCalls = 0, muCalls = 0;
      const auto phase = [](bool shear) {
        if constexpr (std::is_same_v<CoefficientScalar, Complex>)
          return CoefficientScalar(1, shear ? -0.25 : 0.5);
        else
          return CoefficientScalar(shear ? 1.5 : 1);
      };
      const auto first = [&](const Point& point) {
        ++lambdaCalls;
        Real value = 2;
        for (size_t j = 0; j < point.getDimension(); ++j)
          value += Real(0.25) * point(j);
        return phase(false) * (scale * value);
      };
      const auto shear = [&](const Point& point) {
        ++muCalls;
        Real value = 3;
        for (size_t j = 0; j < point.getDimension(); ++j)
          value += Real(0.125) * point(j) * point(j);
        return phase(true) * (scale * value);
      };
      const auto function = [](auto callable) {
        if constexpr (std::is_same_v<CoefficientScalar, Complex>)
          return ComplexFunction(callable);
        else
          return RealFunction(callable);
      };
      const auto lambda = function(first);
      const auto mu = function(shear);
      Optional<std::decay_t<decltype(space)>> separate;
      if constexpr (DistinctSpaces)
        separate.emplace(space);
      auto& testSpace = [&]() -> auto& {
        if constexpr (DistinctSpaces)
          return *separate;
        else
          return space;
      }();
      TrialFunction trial(space);
      TestFunction test(testSpace);
      using Form = decltype(LinearElasticityForm(lambda, mu, trial, test));
      typename Form::Kernel kernel(space, testSpace, lambda, mu);
      EXPECT_EQ(std::addressof(space) == std::addressof(testSpace), !DistinctSpaces);
      const auto check = [&](auto& selected, const Polytope& cell) {
        lambdaCalls = muCalls = 0;
        typename Form::Kernel::MatrixType actual;
        selected.compute(actual, cell);
        if (space.getFiniteElement(cell.getDimension(), cell.getIndex()).getOrder() == 0)
          EXPECT_TRUE((actual.array() == Scalar(0)).all());
        const size_t order =
          2 * space.getFiniteElement(cell.getDimension(), cell.getIndex()).getOrder();
        const auto& formula =
          QF::PolytopeQuadratureFormula::get(order, cell.getGeometry());
        EXPECT_EQ(lambdaCalls, formula.getSize());
        EXPECT_EQ(muCalls, formula.getSize());
        const auto expected = pairwiseOperator(space, cell, lambda, mu);
        ASSERT_EQ(actual.rows(), expected.rows());
        ASSERT_EQ(actual.cols(), expected.cols());
        for (Eigen::Index row = 0; row < actual.rows(); ++row)
        {
          for (Eigen::Index col = 0; col < actual.cols(); ++col)
          {
            ASSERT_EQ(actual(row, col), expected(row, col))
              << "entry=" << row << ',' << col;
          }
        }
      };
      for (const Real value : {Real(1), Real(2)})
      {
        scale = value;
        for (const Index cell : {Index(0), Index(mesh.getCellCount() - 1), Index(0)})
          check(kernel, *mesh.getCell(cell));
      }
      auto copied(kernel);
      scale = 3;
      check(copied, *mesh.getCell());
    }
  }

  class LinearElasticityKernelTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(LinearElasticityKernelTest, PairwiseWeightingAndState)
  {
    for (const bool curved : {false, true})
    {
      auto mesh = Convergence::UniformGrid(GetParam()).makeMesh(3);
      Convergence::CurvedGeometry geometry(mesh);
      if (curved)
        geometry.template install<2>();
      const auto check = [&]<class Scalar>() {
        P1<Math::SpatialVector<Scalar>> linear(mesh, mesh.getSpaceDimension());
        H1<2, Math::SpatialVector<Scalar>> quadratic(
          std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
        checkWeighting(linear);
        checkWeighting(quadratic);
        if constexpr (std::is_same_v<Scalar, Complex>)
        {
          checkWeighting<true>(linear);
          checkWeighting<true>(quadratic);
        }
      };
      check.template operator()<Real>();
      check.template operator()<Complex>();
    }
  }

  /**
   * @brief Distinct space objects preserve the original pairwise result.
   *
   * A copied space has the same basis definition but different object identity.
   * The original pairwise oracle therefore remains applicable, providing a
   * separate-object control for shared-space reuse experiments. Native and
   * quadratic geometry, real/complex weights, coefficient changes, cell changes and copies are
   * checked on all seven positive-dimensional cell families.
   */
  TEST_P(LinearElasticityKernelTest, DistinctSpacesPairwiseWeightingAndState)
  {
    for (const bool curved : {false, true})
    {
      auto mesh = Convergence::UniformGrid(GetParam()).makeMesh(3);
      Convergence::CurvedGeometry geometry(mesh);
      if (curved)
        geometry.template install<2>();
      const auto check = [&]<class Scalar>() {
        P1<Math::SpatialVector<Scalar>> linear(mesh, mesh.getSpaceDimension());
        H1<2, Math::SpatialVector<Scalar>> quadratic(
          std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
        checkWeighting<false, true>(linear);
        checkWeighting<false, true>(quadratic);
        if constexpr (std::is_same_v<Scalar, Complex>)
        {
          checkWeighting<true, true>(linear);
          checkWeighting<true, true>(quadratic);
        }
      };
      check.template operator()<Real>();
      check.template operator()<Complex>();
    }
  }

  /**
   * @brief Constant vector bases retain exactly zero elasticity kernels.
   *
   * Broken and global constants exercise small pair counts independently of
   * the physical dimension. The original derivative-based oracle checks native
   * and quadratic maps, real/complex coefficients, state changes and copies.
   * This is a local zero-derivative contract, not a conforming approximation
   * or convergence claim for the broken constant space.
   */
  TEST_P(LinearElasticityKernelTest, ConstantVectorSpaces)
  {
    for (const bool curved : {false, true})
    {
      auto mesh = Convergence::UniformGrid(GetParam()).makeMesh(3);
      Convergence::CurvedGeometry geometry(mesh);
      if (curved)
        geometry.template install<2>();
      const auto check = [&]<class Scalar>() {
        P0<Math::SpatialVector<Scalar>> broken(mesh, mesh.getSpaceDimension());
        P0g<Math::SpatialVector<Scalar>, LocalMesh> global(
          mesh, mesh.getSpaceDimension());
        checkWeighting(broken);
        checkWeighting(global);
        if constexpr (std::is_same_v<Scalar, Complex>)
        {
          checkWeighting<true>(broken);
          checkWeighting<true>(global);
        }
      };
      check.template operator()<Real>();
      check.template operator()<Complex>();
    }
  }

  /**
   * @brief Complex integrators preserve material values and scalar promotion.
   *
   * P1 and H1 degree two are compared with the independent pairwise oracle
   * on native and quadratic geometry. All real/complex Lamé parameter pairs
   * act on complex fields. Nonreal parameters produce a non-Hermitian operator;
   * coincident spaces alone must not trigger adjoint reconstruction. This test
   * does not require a complex PETSc installation.
   */
  TEST_P(LinearElasticityKernelTest, ComplexIntegralParameterPromotion)
  {
    /// @brief Relative Frobenius roundoff budget on dimensionless reference meshes.
    constexpr Real IntegralTolerance = 1e-11;
    const RealFunction realLambda(Real(2)), realMu(Real(3));
    const ComplexFunction complexLambda(Complex(2, 0.75)), complexMu(Complex(3, -0.25));
    for (const bool curved : {false, true})
    {
      auto mesh = Convergence::UniformGrid(GetParam()).makeMesh(2);
      Convergence::CurvedGeometry geometry(mesh);
      if (curved)
        geometry.template install<2>();
      const auto check = [&](auto& space) {
        TrialFunction trial(space);
        TestFunction test(space);
        const auto verify = [&](const auto& lambda, const auto& mu) {
          auto integrator = LinearElasticityIntegral(trial, test)(lambda, mu);
          const auto cell = *mesh.getCell();
          integrator.setPolytope(cell);
          const auto expected = pairwiseOperator(space, cell, lambda, mu);
          Math::Matrix<Complex> actual(expected.rows(), expected.cols());
          for (Eigen::Index row = 0; row < actual.rows(); ++row)
          {
            for (Eigen::Index col = 0; col < actual.cols(); ++col)
              actual(row, col) = integrator.integrate(col, row);
          }
          EXPECT_LE((actual - expected).norm(),
            IntegralTolerance * std::max(Real(1), expected.norm()));
          using FirstRange = typename std::decay_t<decltype(lambda)>::ScalarType;
          using ShearRange = typename std::decay_t<decltype(mu)>::ScalarType;
          if constexpr (std::is_same_v<FirstRange, Real> &&
            std::is_same_v<ShearRange, Real>)
            EXPECT_EQ((actual - actual.adjoint()).norm(), 0);
          else
            EXPECT_GT((actual - actual.adjoint()).norm(), 0);
        };
        verify(realLambda, realMu);
        verify(realLambda, complexMu);
        verify(complexLambda, realMu);
        verify(complexLambda, complexMu);
      };
      P1<Math::SpatialVector<Complex>> linear(mesh, mesh.getSpaceDimension());
      H1<2, Math::SpatialVector<Complex>> quadratic(
        std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
      check(linear);
      check(quadratic);
    }
  }

  /// @brief A retained kernel resizes and rebinds when the reference cell family changes.
  TEST(LinearElasticityKernelState, MixedCellFamilies)
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(2)
                  .nodes(5)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({1, 1})
                  .vertex({2, 0})
                  .vertex({2, 1})
                  .polytope(Polytope::Type::Triangle, {{0, 1, 2}})
                  .polytope(Polytope::Type::Quadrilateral, {{1, 3, 4, 2}})
                  .finalize();
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    const auto check = [&]<class Scalar>() {
      P1<Math::SpatialVector<Scalar>> linear(mesh, 2);
      H1<2, Math::SpatialVector<Scalar>> quadratic(
        std::integral_constant<size_t, 2>{}, mesh, 2);
      checkWeighting(linear);
      checkWeighting(quadratic);
      if constexpr (std::is_same_v<Scalar, Complex>)
      {
        checkWeighting<true>(linear);
        checkWeighting<true>(quadratic);
      }
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
  }

  /**
   * @brief Higher H1 degrees retain the same weighting and coefficient-read contract.
   *
   * Degrees three through six are checked on quadratic segments. The
   * all-geometry matrix above separately exercises degree two in 2D and 3D;
   * this test does not claim high-degree 3D recertification.
   */
  TEST(LinearElasticityKernelState, HigherDegreeSegments)
  {
    auto mesh = Convergence::UniformGrid(Polytope::Type::Segment).makeMesh(3);
    Convergence::CurvedGeometry geometry(mesh);
    geometry.template install<2>();
    const auto check = [&]<size_t K>() {
      H1<K, Math::SpatialVector<Real>> real(std::integral_constant<size_t, K>{}, mesh, 1);
      H1<K, Math::SpatialVector<Complex>> complex(
        std::integral_constant<size_t, K>{}, mesh, 1);
      checkWeighting(real);
      checkWeighting(complex);
      checkWeighting<true>(complex);
    };
    check.template operator()<3>();
    check.template operator()<4>();
    check.template operator()<5>();
    check.template operator()<6>();
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, LinearElasticityKernelTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge));
}
