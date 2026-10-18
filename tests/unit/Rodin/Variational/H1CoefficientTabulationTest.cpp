/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <memory>
#include <type_traits>
#include <vector>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

namespace Rodin::Tests::Unit
{
  enum class H1CoefficientForm
  {
    Mass,
    Gradient,
    Jacobian,
    RankOneShared,
    RankOneDistinct,
    Convection,
    Load,
    MassRight
  };

  /**
   * @brief Certifies reference-table borrowing across coefficient callbacks.
   *
   * ## Mathematical foundation
   *
   * A coefficient that returns a fixed value must define the same discrete
   * form whether or not its evaluation queries unrelated reference tables.
   * The expected and actual forms have identical quadrature and arithmetic;
   * every entry is therefore compared exactly, without a numerical tolerance.
   * Shared and distinct coefficient identities use matching kernel branches
   * in the two forms; triangular mirroring is not compared with full assembly.
   * The same callable types are used in both forms. A captured flag changes
   * only whether table queries are performed, not coefficient representation.
   * Real and complex quadratic scalar/vector spaces are exercised on a segment.
   * Empty right-scaled vector spaces retain coefficient evaluation at each
   * sample even though the local operator has no entries.
   *
   * ## Architecture
   *
   * Nine owned formulas have distinct logical identities and exceed the
   * documented eight-entry reference store. Their sample counts exceed the
   * bound rule's count, so a broken borrow reads replaced samples rather than
   * relying on an out-of-bounds access. Coefficient values remain unchanged.
   * Callback counts also reject a fallback to per-entry generic evaluation.
   */
  class H1CoefficientTabulationTest : public ::testing::TestWithParam<H1CoefficientForm>
  {
    protected:
      /**
       * @brief Compares one constant-coefficient form with its cache-querying counterpart.
       * @tparam Scalar Real or complex field scalar.
       * @param form Local operator family to exercise.
       */
      template <class Scalar>
      static void check(H1CoefficientForm form)
      {
        SCOPED_TRACE((std::is_same_v<Scalar, Complex> ? "Complex" : "Real"));
        using namespace Geometry;
        using namespace Variational;
        constexpr size_t NestedRules = 9;
        constexpr size_t BoundOrder = 4;
        constexpr size_t NestedOrder = 16;
        auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
        const auto cell = *mesh.getCell();
        const auto& bound =
          QF::PolytopeQuadratureFormula::get(BoundOrder, Polytope::Type::Segment);
        std::vector<std::unique_ptr<const QF::QuadratureFormulaBase>> nested;
        for (size_t i = 0; i < NestedRules; ++i)
        {
          nested.emplace_back(
            QF::PolytopeQuadratureFormula::build(NestedOrder, Polytope::Type::Segment));
          ASSERT_GT(nested.back()->getSize(), bound.getSize());
          for (size_t previous = 0; previous < i; ++previous)
          {
            ASSERT_NE(
              nested[previous]->getCacheIdentity(), nested.back()->getCacheIdentity());
          }
        }
        const H1Element<2, Scalar> element(Polytope::Type::Segment);
        size_t evaluations = 0;
        bool queryTables = false;
        const auto evict = [&]() {
          ++evaluations;
          if (queryTables)
          {
            for (const auto& formula : nested)
              element.getTabulation(*formula);
          }
        };
        const auto scalarCallable = [&](const Geometry::Point&) -> Scalar {
          evict();
          return Scalar(2);
        };
        const auto scalarCoefficient = [&]() {
          if constexpr (std::is_same_v<Scalar, Real>)
            return RealFunction(scalarCallable);
          else
            return ComplexFunction(scalarCallable);
        }();
        Math::Vector<Scalar> direction(1);
        direction[0] = Scalar(2);
        const auto vectorCoefficient =
          VectorFunction(size_t(1), [&](const Geometry::Point&) {
            evict();
            return direction;
          });
        Math::SpatialMatrix<Scalar> matrix(2, 1);
        matrix(0, 0) = Scalar(2);
        matrix(1, 0) = Scalar(-1);
        const auto matrixCallable = [&](const Geometry::Point&) {
          evict();
          return matrix;
        };
        const auto matrixCoefficient = MatrixFunction(2, 1, matrixCallable);
        const auto otherMatrixCoefficient = MatrixFunction(2, 1, matrixCallable);
        H1<2, Scalar> scalarSpace(std::integral_constant<size_t, 2>{}, mesh);
        H1<2, Math::SpatialVector<Scalar>> vectorSpace(
          std::integral_constant<size_t, 2>{}, mesh, 2);
        TrialFunction u(scalarSpace);
        TestFunction v(scalarSpace);
        TrialFunction w(vectorSpace);
        TestFunction z(vectorSpace);
        const auto compare = [&](auto& expected, auto& actual, size_t count) {
          expected.setOrder(BoundOrder);
          actual.setOrder(BoundOrder);
          queryTables = false;
          expected.setPolytope(cell);
          evaluations = 0;
          queryTables = true;
          actual.setPolytope(cell);
          ASSERT_GE(evaluations, bound.getSize());
          ASSERT_LE(evaluations, 2 * bound.getSize());
          for (size_t te = 0; te < count; ++te)
          {
            for (size_t tr = 0; tr < count; ++tr)
            {
              SCOPED_TRACE(::testing::Message() << "test=" << te << ", trial=" << tr);
              EXPECT_EQ(actual.integrate(tr, te), expected.integrate(tr, te));
            }
          }
        };
        const size_t scalarCount =
          scalarSpace.getFiniteElement(1, cell.getIndex()).getCount();
        const size_t vectorCount =
          vectorSpace.getFiniteElement(1, cell.getIndex()).getCount();
        switch (form)
        {
          case H1CoefficientForm::Mass:
          {
            auto expected = Integral(scalarCoefficient * w, z);
            auto actual = Integral(scalarCoefficient * w, z);
            compare(expected, actual, vectorCount);
            break;
          }
          case H1CoefficientForm::MassRight:
          {
            auto expected = Integral(w * scalarCoefficient, z);
            auto actual = Integral(w * scalarCoefficient, z);
            compare(expected, actual, vectorCount);
            H1<2, Math::SpatialVector<Scalar>> emptySpace(
              std::integral_constant<size_t, 2>{}, mesh, 0);
            TrialFunction emptyTrial(emptySpace);
            TestFunction emptyTest(emptySpace);
            auto emptyExpected = Integral(emptyTrial * scalarCoefficient, emptyTest);
            auto emptyActual = Integral(emptyTrial * scalarCoefficient, emptyTest);
            compare(emptyExpected, emptyActual, 0);
            break;
          }
          case H1CoefficientForm::Gradient:
          {
            auto expected = Integral(scalarCoefficient * Grad(u), Grad(v));
            auto actual = Integral(scalarCoefficient * Grad(u), Grad(v));
            compare(expected, actual, scalarCount);
            break;
          }
          case H1CoefficientForm::Jacobian:
          {
            auto expected = Integral(scalarCoefficient * Jacobian(w), Jacobian(z));
            auto actual = Integral(scalarCoefficient * Jacobian(w), Jacobian(z));
            compare(expected, actual, vectorCount);
            break;
          }
          case H1CoefficientForm::RankOneShared:
          {
            auto expected = Integral(
              Dot(matrixCoefficient, Jacobian(w)), Dot(matrixCoefficient, Jacobian(z)));
            auto actual = Integral(
              Dot(matrixCoefficient, Jacobian(w)), Dot(matrixCoefficient, Jacobian(z)));
            compare(expected, actual, vectorCount);
            break;
          }
          case H1CoefficientForm::RankOneDistinct:
          {
            auto expected = Integral(Dot(matrixCoefficient, Jacobian(w)),
              Dot(otherMatrixCoefficient, Jacobian(z)));
            auto actual = Integral(Dot(matrixCoefficient, Jacobian(w)),
              Dot(otherMatrixCoefficient, Jacobian(z)));
            compare(expected, actual, vectorCount);
            break;
          }
          case H1CoefficientForm::Convection:
          {
            auto expected = Integral(Jacobian(w) * vectorCoefficient, z);
            auto actual = Integral(Jacobian(w) * vectorCoefficient, z);
            compare(expected, actual, vectorCount);
            break;
          }
          case H1CoefficientForm::Load:
          {
            auto expected = Integral(scalarCoefficient, v);
            auto actual = Integral(scalarCoefficient, v);
            expected.setOrder(BoundOrder);
            actual.setOrder(BoundOrder);
            queryTables = false;
            expected.setPolytope(cell);
            evaluations = 0;
            queryTables = true;
            actual.setPolytope(cell);
            ASSERT_EQ(evaluations, bound.getSize());
            for (size_t te = 0; te < scalarCount; ++te)
            {
              SCOPED_TRACE(te);
              EXPECT_EQ(actual.integrate(te), expected.integrate(te));
            }
            break;
          }
        }
      }
  };

  TEST_P(H1CoefficientTabulationTest, ConstantValueIsIndependentOfCacheQueries)
  {
    check<Real>(GetParam());
    check<Complex>(GetParam());
  }

  INSTANTIATE_TEST_SUITE_P(Forms, H1CoefficientTabulationTest,
    ::testing::Values(H1CoefficientForm::Mass, H1CoefficientForm::Gradient,
      H1CoefficientForm::Jacobian, H1CoefficientForm::RankOneShared,
      H1CoefficientForm::RankOneDistinct, H1CoefficientForm::Convection,
      H1CoefficientForm::Load, H1CoefficientForm::MassRight));
}
