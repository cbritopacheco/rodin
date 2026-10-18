/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <array>
#include <type_traits>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "QuadratureReference.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Certifies logical component addresses of H1 mass quadrature.
   *
   * ## Mathematical contract
   *
   * A scalar weight gives @f$A_{bm+e,am+c}=\delta_{ec}S_{ba}@f$, where
   * @f$S@f$ is the scalar quadrature operator. The scalar space has one
   * component and does not use the intercomponent stride; its entries qualify
   * every vector/matrix component,
   * including structural zeros. Equality is exact, without a tolerance or
   * coordinate-based identification. For a row-major local matrix with
   * @f$n_r@f$ trial DOFs, successive diagonal components differ by
   * @f$n_r+1@f$ scalar indices. Rectangular mixed-degree operators use
   * @f$n_r@f$, not the test-space count, as their row stride.
   *
   * One cell of every reference family exercises real/complex, scalar-weighted
   * and unit forms, empty vectors, one through three vector components and
   * rectangular matrix ranges. Orders zero, two, four and zero again test
   * formula replacement and return to the initial rule while the coefficient
   * scale changes. Point mass is nonzero counting-measure algebra. These are
   * discrete-operator identities, not integration-exactness or PDE-rate claims.
   * All degrees exercise two-component matrices; six-component matrices are
   * exercised through degree three to isolate the component traversal without
   * constructing unnecessarily large dense unit-test operators. The separate
   * opt-in benchmark retains the higher-degree, wider-range cases.
   */
  class H1MassQuadratureTest : public ::testing::Test
  {
    protected:
      /**
       * @brief Preserves coefficient promotion in right-weighted mass entries.
       *
       * Complex fields with real coefficients are compared with the original
       * entrywise algorithm. Complex coefficients on real trial/test fields
       * are excluded by the existing equal-range contract of the dot product.
       * Mixed degrees make the local operator rectangular; scalar, vector and
       * matrix fields retain the same trial-first conjugation convention.
       *
       * @tparam Scalar Finite-element scalar.
       * @tparam Weight Scalar coefficient range.
       */
      template <class Scalar, class Weight>
      static void checkRightPromotion()
      {
        using namespace Geometry;
        using namespace Variational;
        for (const auto geometry :
          {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
        {
          SCOPED_TRACE(static_cast<int>(geometry));
          SCOPED_TRACE((std::is_same_v<Scalar, Real> ? "Real field" : "Complex field"));
          const size_t dimension = Polytope::Traits(geometry).getDimension();
          auto mesh = [&]() {
            if (dimension == 0)
              return LocalMesh::Builder()
                .initialize(3)
                .nodes(1)
                .vertex({0, 0, 0})
                .finalize();
            return LocalMesh::UniformGrid(
              geometry, Array<size_t>::Constant(dimension, 2));
          }();
          for (size_t d = 1; d <= dimension; ++d)
          {
            for (size_t lower = 0; lower < d; ++lower)
              mesh.getConnectivity().compute(d, lower);
          }
          const auto cell = *mesh.getCell();
          const auto coefficient = [&]() {
            const auto value = [](const Point& point) -> Weight {
              const Real x = point.getDimension() == 0 ? Real(0) : point.x();
              if constexpr (std::is_same_v<Weight, Complex>)
                return Complex(Real(2) + x, Real(0.5));
              else
                return Real(2) + x;
            };
            if constexpr (std::is_same_v<Weight, Complex>)
              return ComplexFunction(value);
            else
              return RealFunction(value);
          }();
          const auto checkPair = [&](auto& trial, auto& test, bool specialized) {
            TrialFunction u(trial);
            TestFunction v(test);
            auto actual = Integral(u * coefficient, v);
            ReferenceIntegral reference(actual);
            EXPECT_EQ(decltype(actual)::Specialized, specialized);
            for (const size_t order : {size_t(0), size_t(4)})
            {
              actual.setOrder(order);
              reference.setOrder(order);
              reference.setPolytope(cell);
              actual.setPolytope(cell);
              for (size_t tr = 0; tr < u.getDOFs(cell); ++tr)
              {
                for (size_t te = 0; te < v.getDOFs(cell); ++te)
                {
                  SCOPED_TRACE(::testing::Message() << "trial=" << tr << ", test=" << te);
                  ASSERT_EQ(actual.integrate(tr, te), reference.integrate(tr, te));
                }
              }
            }
          };
          H1<2, Scalar> scalarTrial(std::integral_constant<size_t, 2>{}, mesh);
          H1<1, Scalar> scalarTest(std::integral_constant<size_t, 1>{}, mesh);
          checkPair(scalarTrial, scalarTest, false);
          H1<2, Math::SpatialVector<Scalar>> vectorTrial(
            std::integral_constant<size_t, 2>{}, mesh, 2);
          H1<1, Math::SpatialVector<Scalar>> vectorTest(
            std::integral_constant<size_t, 1>{}, mesh, 2);
          checkPair(vectorTrial, vectorTest, true);
          H1<2, Math::SpatialMatrix<Scalar>> matrixTrial(
            std::integral_constant<size_t, 2>{}, mesh, 2, 3);
          H1<1, Math::SpatialMatrix<Scalar>> matrixTest(
            std::integral_constant<size_t, 1>{}, mesh, 2, 3);
          checkPair(matrixTrial, matrixTest, true);
        }
      }

      /**
       * @brief Requires component-aware right scaling without commuting other expressions.
       *
       * Scalar fields retain the generic cached-trial route. Vector and matrix
       * carriers use component separation rather than dense range contraction.
       *
       * @tparam Scalar Real or complex finite-element scalar.
       */
      template <class Scalar>
      static void checkRightDispatch()
      {
        using namespace Geometry;
        using namespace Variational;
        auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
        mesh.getConnectivity().compute(1, 0);
        const auto coefficient = [&]() {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction([](const Point&) { return Complex(2, 0.5); });
          else
            return RealFunction([](const Point&) { return Real(2); });
        }();
        const auto checkPair = [&](auto& trial, auto& test, bool specialized) {
          TrialFunction u(trial);
          TestFunction v(test);
          auto right = Integral(u * coefficient, v);
          EXPECT_EQ(decltype(right)::Specialized, specialized);
          auto nested = Integral((u * coefficient) * coefficient, v);
          EXPECT_FALSE(decltype(nested)::Specialized);
        };
        H1<2, Scalar> scalarTrial(std::integral_constant<size_t, 2>{}, mesh);
        H1<1, Scalar> scalarTest(std::integral_constant<size_t, 1>{}, mesh);
        checkPair(scalarTrial, scalarTest, false);
        H1<2, Math::SpatialVector<Scalar>> vectorTrial(
          std::integral_constant<size_t, 2>{}, mesh, 2);
        H1<1, Math::SpatialVector<Scalar>> vectorTest(
          std::integral_constant<size_t, 1>{}, mesh, 2);
        checkPair(vectorTrial, vectorTest, true);
        H1<2, Math::SpatialMatrix<Scalar>> matrixTrial(
          std::integral_constant<size_t, 2>{}, mesh, 2, 3);
        H1<1, Math::SpatialMatrix<Scalar>> matrixTest(
          std::integral_constant<size_t, 1>{}, mesh, 2, 3);
        checkPair(matrixTrial, matrixTest, true);
        Math::SpatialMatrix<Scalar> matrix(3, 3);
        matrix.setZero();
        for (size_t i = 0; i < 3; ++i)
          matrix(i, i) = Scalar(1);
        matrix(0, 1) = Scalar(2);
        TrialFunction u(matrixTrial);
        TestFunction v(matrixTest);
        auto coupled = Integral(u * MatrixFunction(matrix), v);
        EXPECT_FALSE(decltype(coupled)::Specialized);
      }

      /**
       * @brief Checks equal or mixed degrees against their scalar mass operator.
       * @tparam TrialDegree Trial polynomial degree.
       * @tparam TestDegree Test polynomial degree.
       * @tparam Scalar Real or complex field scalar.
       */
      template <size_t TrialDegree, size_t TestDegree, class Scalar>
      static void check()
      {
        using namespace Geometry;
        using namespace Variational;
        for (const auto geometry :
          {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
        {
          SCOPED_TRACE(static_cast<int>(geometry));
          SCOPED_TRACE(TrialDegree);
          SCOPED_TRACE(TestDegree);
          const size_t dimension = Polytope::Traits(geometry).getDimension();
          auto mesh = [&]() {
            if (dimension == 0)
              return LocalMesh::Builder()
                .initialize(3)
                .nodes(1)
                .vertex({0, 0, 0})
                .finalize();
            return LocalMesh::UniformGrid(
              geometry, Array<size_t>::Constant(dimension, 2));
          }();
          for (size_t d = 1; d <= dimension; ++d)
          {
            for (size_t lower = 0; lower < d; ++lower)
              mesh.getConnectivity().compute(d, lower);
          }
          const auto cell = *mesh.getCell();
          H1<TrialDegree, Scalar> scalarTrial(
            std::integral_constant<size_t, TrialDegree>{}, mesh);
          H1<TestDegree, Scalar> scalarTest(
            std::integral_constant<size_t, TestDegree>{}, mesh);
          TrialFunction scalarU(scalarTrial);
          TestFunction scalarV(scalarTest);
          Real scale = 1;
          const auto coefficient = [&]() {
            if constexpr (std::is_same_v<Scalar, Complex>)
              return ComplexFunction([&](const Point& point) {
                const Real x = point.getDimension() == 0 ? Real(0) : point.x();
                return scale * Complex(Real(2) + x, Real(0.5));
              });
            else
              return RealFunction([&](const Point& point) {
                const Real x = point.getDimension() == 0 ? Real(0) : point.x();
                return scale * (Real(2) + x);
              });
          }();
          const auto checkRange = [&](auto& trial, auto& test, size_t components) {
            SCOPED_TRACE(components);
            TrialFunction u(trial);
            TestFunction v(test);
            const size_t nr = scalarU.getDOFs(cell), nt = scalarV.getDOFs(cell);
            ASSERT_EQ(u.getDOFs(cell), components * nr);
            ASSERT_EQ(v.getDOFs(cell), components * nt);
            const auto checkIntegral = [&](auto& actual, auto& reference) {
              constexpr std::array<size_t, 5> Orders{
                0, 2, 4, TrialDegree + TestDegree + 1, 0};
              for (size_t step = 0; step < Orders.size(); ++step)
              {
                SCOPED_TRACE(step);
                scale = Real(step + 1);
                actual.setOrder(Orders[step]);
                reference.setOrder(Orders[step]);
                reference.setPolytope(cell);
                actual.setPolytope(cell);
                bool nonzero = false;
                for (size_t b = 0; b < nt; ++b)
                {
                  for (size_t a = 0; a < nr; ++a)
                  {
                    const Scalar expected = reference.integrate(a, b);
                    nonzero = nonzero || expected != Scalar(0);
                    for (size_t e = 0; e < components; ++e)
                    {
                      for (size_t c = 0; c < components; ++c)
                      {
                        ASSERT_EQ(
                          actual.integrate(a * components + c, b * components + e),
                          c == e ? expected : Scalar(0));
                      }
                    }
                  }
                }
                ASSERT_TRUE(nonzero);
                actual.setPolytope(cell);
                for (size_t b = 0; b < nt; ++b)
                {
                  for (size_t a = 0; a < nr; ++a)
                  {
                    const Scalar expected = reference.integrate(a, b);
                    for (size_t c = 0; c < components; ++c)
                    {
                      ASSERT_EQ(actual.integrate(a * components + c, b * components + c),
                        expected);
                    }
                  }
                }
              }
            };
            {
              SCOPED_TRACE("unit mass");
              auto actual = Integral(u, v);
              auto reference = Integral(scalarU, scalarV);
              checkIntegral(actual, reference);
            }
            {
              SCOPED_TRACE("left scalar weight");
              auto actual = Integral(coefficient * u, v);
              auto reference = Integral(coefficient * scalarU, scalarV);
              checkIntegral(actual, reference);
            }
            {
              SCOPED_TRACE("right scalar weight");
              auto actual = Integral(u * coefficient, v);
              const auto original = Integral(scalarU * coefficient, scalarV);
              ReferenceIntegral reference(original);
              checkIntegral(actual, reference);
            }
          };
          checkRange(scalarTrial, scalarTest, 1);
          for (const size_t components : {size_t(0), size_t(1), size_t(2), size_t(3)})
          {
            H1<TrialDegree, Math::SpatialVector<Scalar>> trial(
              std::integral_constant<size_t, TrialDegree>{}, mesh, components);
            H1<TestDegree, Math::SpatialVector<Scalar>> test(
              std::integral_constant<size_t, TestDegree>{}, mesh, components);
            checkRange(trial, test, components);
          }
          for (const auto shape : {std::array<size_t, 2>{1, 1}, {1, 2}, {2, 3}})
          {
            if constexpr (TrialDegree > 3 || TestDegree > 3)
            {
              if (shape[0] == 2)
                continue;
            }
            SCOPED_TRACE(
              ::testing::Message() << "matrix " << shape[0] << 'x' << shape[1]);
            H1<TrialDegree, Math::SpatialMatrix<Scalar>> trial(
              std::integral_constant<size_t, TrialDegree>{}, mesh, shape[0], shape[1]);
            H1<TestDegree, Math::SpatialMatrix<Scalar>> test(
              std::integral_constant<size_t, TestDegree>{}, mesh, shape[0], shape[1]);
            checkRange(trial, test, shape[0] * shape[1]);
          }
        }
      }
  };

  TEST_F(H1MassQuadratureTest, RightScalarDispatch)
  {
    checkRightDispatch<Real>();
    checkRightDispatch<Complex>();
  }

  TEST_F(H1MassQuadratureTest, RightScalarPromotion)
  {
    checkRightPromotion<Complex, Real>();
  }

  TEST_F(H1MassQuadratureTest, EqualDegrees)
  {
    check<1, 1, Real>();
    check<2, 2, Real>();
    check<3, 3, Real>();
    check<4, 4, Real>();
    check<5, 5, Real>();
    check<6, 6, Real>();
    check<1, 1, Complex>();
    check<2, 2, Complex>();
    check<3, 3, Complex>();
    check<4, 4, Complex>();
    check<5, 5, Complex>();
    check<6, 6, Complex>();
  }

  TEST_F(H1MassQuadratureTest, MixedDegrees)
  {
    check<1, 2, Real>();
    check<2, 1, Real>();
    check<2, 3, Real>();
    check<3, 2, Real>();
    check<3, 6, Real>();
    check<6, 3, Real>();
    check<1, 2, Complex>();
    check<2, 1, Complex>();
    check<2, 3, Complex>();
    check<3, 2, Complex>();
    check<3, 6, Complex>();
    check<6, 3, Complex>();
  }
}
