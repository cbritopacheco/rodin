/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <memory>
#include <type_traits>
#include <vector>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Certifies empty component ranges without accessing nonexistent entries.
   *
   * For @f$V_h=\{0\}@f$, mass and Jacobian forms have no local entries,
   * independently of the scalar polynomial degrees. Binding must not construct
   * scalar reference tables or contract scalar basis pairs. In particular,
   * equality of two zero DOF counts does not imply equal scalar basis counts.
   *
   * Nine distinct owned formulas evict the documented eight-entry reference
   * store before each check. Absence is checked by logical formula identity;
   * no floating-point tolerance or invalid matrix index is used.
   */
  class H1EmptyBlocksTest : public ::testing::Test
  {
    protected:
      /**
       * @brief Checks equal or mixed scalar degrees on all reference-cell families.
       * @tparam TrialDegree Scalar degree of the empty trial space.
       * @tparam TestDegree Scalar degree of the empty test space.
       * @tparam Scalar Real or complex scalar field.
       */
      template <size_t TrialDegree, size_t TestDegree, class Scalar>
      static void check()
      {
        using namespace Geometry;
        using namespace Variational;
        constexpr size_t EvictionRules = 9;
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
            {
              auto builder = LocalMesh::Builder();
              builder.initialize(3).nodes(1).vertex({0, 0, 0});
              return builder.finalize();
            }
            const Array<size_t> shape = Array<size_t>::Constant(dimension, 2);
            return LocalMesh::UniformGrid(geometry, shape);
          }();
          for (size_t d = 1; d <= dimension; ++d)
          {
            for (size_t lower = 0; lower < d; ++lower)
              mesh.getConnectivity().compute(d, lower);
          }
          const auto cell = *mesh.getCell();
          H1<TrialDegree, Math::SpatialVector<Scalar>> trial(
            std::integral_constant<size_t, TrialDegree>{}, mesh, 0);
          H1<TestDegree, Math::SpatialVector<Scalar>> test(
            std::integral_constant<size_t, TestDegree>{}, mesh, 0);
          TrialFunction u(trial);
          TestFunction v(test);
          ASSERT_EQ(u.getDOFs(cell), 0);
          ASSERT_EQ(v.getDOFs(cell), 0);
          const H1Element<TrialDegree, Scalar> trialElement(geometry);
          const H1Element<TestDegree, Scalar> testElement(geometry);
          using TrialTables =
            ReferenceTabulationCache<H1Element<TrialDegree, Scalar>, Scalar>;
          using TestTables =
            ReferenceTabulationCache<H1Element<TestDegree, Scalar>, Scalar>;
          const auto& bound = QF::PolytopeQuadratureFormula::get(0, geometry);
          std::vector<const QF::QuadratureFormula*> eviction;
          for (size_t i = 0; i < EvictionRules; ++i)
          {
            eviction.emplace_back(&QF::PolytopeQuadratureFormula::get(i + 1, geometry));
            trialElement.getTabulation(*eviction.back());
            testElement.getTabulation(*eviction.back());
          }
          ASSERT_FALSE(TrialTables::find(trialElement, bound));
          ASSERT_FALSE(TestTables::find(testElement, bound));
          auto mass = Integral(u, v);
          mass.setOrder(0);
          mass.setPolytope(cell);
          ASSERT_FALSE(TrialTables::find(trialElement, bound));
          ASSERT_FALSE(TestTables::find(testElement, bound));
          auto stiffness = Integral(Jacobian(u), Jacobian(v));
          stiffness.setOrder(0);
          stiffness.setPolytope(cell);
          ASSERT_FALSE(TrialTables::find(trialElement, bound));
          ASSERT_FALSE(TestTables::find(testElement, bound));
        }
      }
  };

  TEST_F(H1EmptyBlocksTest, EqualDegrees)
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

  TEST_F(H1EmptyBlocksTest, MixedDegrees)
  {
    check<2, 3, Real>();
    check<3, 2, Real>();
    check<3, 6, Real>();
    check<6, 3, Real>();
    check<2, 3, Complex>();
    check<3, 2, Complex>();
    check<3, 6, Complex>();
    check<6, 3, Complex>();
  }

  /**
   * @brief Certifies complete component blocks independently of triangular copying.
   *
   * For the scalar quadrature block @f$S@f$, every local component entry must
   * satisfy @f$A_{bm+e,am+c}=\delta_{ec}S_{ba}@f$. The one-component rule
   * supplies @f$S@f$ without invoking component replication. Both triangles,
   * all cross-component zeros and rectangular mixed-degree operators are
   * checked exactly after quadrature-order replacement and repeated binding.
   * The indices are the H1 element's local component convention, not an
   * assumption about global assembly numbering.
   *
   * Vector and rectangular matrix ranges exercise real and complex fields on
   * every reference-cell family. The existing embedded-centroid test retains
   * its independent public-basis derivative oracle. Order zero denotes the
   * discrete one-point operator, not higher-degree integration exactness.
   */
  class H1JacobianComponentBlocksTest : public ::testing::Test
  {
    protected:
      /**
       * @brief Checks equal or mixed degrees against their scalar block.
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
            {
              auto builder = LocalMesh::Builder();
              builder.initialize(3).nodes(1).vertex({0, 0, 0});
              return builder.finalize();
            }
            return LocalMesh::UniformGrid(
              geometry, Array<size_t>::Constant(dimension, 2));
          }();
          for (size_t d = 1; d <= dimension; ++d)
          {
            for (size_t lower = 0; lower < d; ++lower)
              mesh.getConnectivity().compute(d, lower);
          }
          using Vector = Math::SpatialVector<Scalar>;
          H1<TrialDegree, Vector> oneTrial(
            std::integral_constant<size_t, TrialDegree>{}, mesh, 1);
          H1<TestDegree, Vector> oneTest(
            std::integral_constant<size_t, TestDegree>{}, mesh, 1);
          TrialFunction oneU(oneTrial);
          TestFunction oneV(oneTest);
          auto scalar = Integral(Jacobian(oneU), Jacobian(oneV));
          const auto checkRange = [&](auto& trial, auto& test, size_t components) {
            SCOPED_TRACE(components);
            TrialFunction u(trial);
            TestFunction v(test);
            auto integral = Integral(Jacobian(u), Jacobian(v));
            for (const size_t order : {size_t(0), size_t(4), size_t(0)})
            {
              SCOPED_TRACE(order);
              scalar.setOrder(order);
              integral.setOrder(order);
              for (auto cell = mesh.getCell(); cell; ++cell)
              {
                scalar.setPolytope(*cell);
                integral.setPolytope(*cell);
                const size_t nr = oneU.getDOFs(*cell), nt = oneV.getDOFs(*cell);
                ASSERT_EQ(u.getDOFs(*cell), components * nr);
                ASSERT_EQ(v.getDOFs(*cell), components * nt);
                bool nonzero = false;
                for (size_t b = 0; b < nt; ++b)
                {
                  for (size_t a = 0; a < nr; ++a)
                  {
                    const Scalar value = scalar.integrate(a, b);
                    nonzero = nonzero || value != Scalar(0);
                    for (size_t e = 0; e < components; ++e)
                    {
                      for (size_t c = 0; c < components; ++c)
                      {
                        ASSERT_EQ(
                          integral.integrate(a * components + c, b * components + e),
                          c == e ? value : Scalar(0));
                      }
                    }
                  }
                }
                if (dimension > 0)
                {
                  ASSERT_TRUE(nonzero);
                }
              }
            }
          };
          for (const size_t components : {size_t(0), size_t(1), size_t(2), size_t(3)})
          {
            H1<TrialDegree, Vector> trial(
              std::integral_constant<size_t, TrialDegree>{}, mesh, components);
            H1<TestDegree, Vector> test(
              std::integral_constant<size_t, TestDegree>{}, mesh, components);
            checkRange(trial, test, components);
          }
          for (const size_t rows : {size_t(1), size_t(2)})
          {
            const size_t cols = 3 - rows;
            SCOPED_TRACE(::testing::Message() << "matrix " << rows << 'x' << cols);
            H1<TrialDegree, Math::SpatialMatrix<Scalar>> trial(
              std::integral_constant<size_t, TrialDegree>{}, mesh, rows, cols);
            H1<TestDegree, Math::SpatialMatrix<Scalar>> test(
              std::integral_constant<size_t, TestDegree>{}, mesh, rows, cols);
            checkRange(trial, test, rows * cols);
          }
        }
      }
  };

  TEST_F(H1JacobianComponentBlocksTest, EqualDegrees)
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

  TEST_F(H1JacobianComponentBlocksTest, MixedDegrees)
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

  /**
   * @brief Checks every component block of an embedded one-point Jacobian form.
   *
   * ## Mathematical foundation
   *
   * The segment map is @f$x(r)=(r,r,0)@f$, with measure factor @f$\sqrt2@f$
   * and inverse row @f$(1/2,1/2,0)@f$. At the centroid, physical basis
   * Jacobians follow directly from the public reference-basis Jacobians.
   * Every trial/test entry is checked exactly, including cross-component
   * zeros, for equal and mixed degrees and real/complex scalar fields.
   * Component counts one through three are independent of the embedding.
   * No assumed local component numbering or coordinate-based identity is used.
   * This is the implemented one-point quadrature, not a claim that one point
   * integrates the complete higher-degree stiffness operator exactly.
   * Scalar degrees one through six and both mixed-degree orientations are
   * included; the same scalar contraction must populate every component block.
   */
  TEST(H1JacobianBlocks, EmbeddedCentroidMatchesReferenceBasisExactly)
  {
    using namespace Geometry;
    using namespace Variational;
    auto builder = LocalMesh::Builder();
    builder.initialize(3).nodes(2).vertex({0, 0, 0}).vertex({1, 1, 0});
    auto mesh = builder.polytope(Polytope::Type::Segment, {0, 1}).finalize();
    const auto cell = *mesh.getCell();
    const Math::SpatialPoint centroid{Real(0.5)};
    const Math::SpatialVector<Real> inverseRow{Real(0.5), Real(0.5), Real(0)};
    const Real measure = std::sqrt(Real(2));
    const auto check = [&]<size_t TrialDegree, size_t TestDegree, class Scalar>() {
      for (const size_t components : {size_t(1), size_t(2), size_t(3)})
      {
        SCOPED_TRACE(TrialDegree);
        SCOPED_TRACE(TestDegree);
        SCOPED_TRACE(components);
        H1<TrialDegree, Math::SpatialVector<Scalar>> trial(
          std::integral_constant<size_t, TrialDegree>{}, mesh, components);
        H1<TestDegree, Math::SpatialVector<Scalar>> test(
          std::integral_constant<size_t, TestDegree>{}, mesh, components);
        TrialFunction u(trial);
        TestFunction v(test);
        auto integral = Integral(Jacobian(u), Jacobian(v));
        integral.setOrder(0);
        integral.setPolytope(cell);
        const auto& trialElement = trial.getFiniteElement(1, cell.getIndex());
        const auto& testElement = test.getFiniteElement(1, cell.getIndex());
        for (size_t te = 0; te < testElement.getCount(); ++te)
        {
          const auto testReference = testElement.getBasis(te).getJacobian()(centroid);
          ASSERT_EQ(testReference.rows(), components);
          ASSERT_EQ(testReference.cols(), 1);
          for (size_t tr = 0; tr < trialElement.getCount(); ++tr)
          {
            SCOPED_TRACE(::testing::Message() << "test=" << te << ", trial=" << tr);
            const auto trialReference = trialElement.getBasis(tr).getJacobian()(centroid);
            ASSERT_EQ(trialReference.rows(), components);
            ASSERT_EQ(trialReference.cols(), 1);
            Scalar expected = 0;
            for (size_t component = 0; component < components; ++component)
            {
              for (size_t axis = 0; axis < inverseRow.size(); ++axis)
              {
                const Scalar trialValue = trialReference(component, 0) * inverseRow(axis);
                const Scalar testValue = testReference(component, 0) * inverseRow(axis);
                if constexpr (std::is_same_v<Scalar, Complex>)
                  expected += trialValue * std::conj(testValue);
                else
                  expected += trialValue * testValue;
              }
            }
            expected *= measure;
            EXPECT_EQ(integral.integrate(tr, te), expected);
          }
        }
      }
    };
    check.template operator()<1, 1, Real>();
    check.template operator()<2, 2, Real>();
    check.template operator()<1, 2, Real>();
    check.template operator()<2, 1, Real>();
    check.template operator()<3, 3, Real>();
    check.template operator()<4, 4, Real>();
    check.template operator()<5, 5, Real>();
    check.template operator()<6, 6, Real>();
    check.template operator()<2, 3, Real>();
    check.template operator()<3, 2, Real>();
    check.template operator()<3, 6, Real>();
    check.template operator()<6, 3, Real>();
    check.template operator()<1, 1, Complex>();
    check.template operator()<2, 2, Complex>();
    check.template operator()<1, 2, Complex>();
    check.template operator()<2, 1, Complex>();
    check.template operator()<3, 3, Complex>();
    check.template operator()<4, 4, Complex>();
    check.template operator()<5, 5, Complex>();
    check.template operator()<6, 6, Complex>();
    check.template operator()<2, 3, Complex>();
    check.template operator()<3, 2, Complex>();
    check.template operator()<3, 6, Complex>();
    check.template operator()<6, 3, Complex>();
  }
}
