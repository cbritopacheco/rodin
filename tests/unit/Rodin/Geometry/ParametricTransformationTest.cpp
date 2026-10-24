/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <thread>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Test/Random/RandomPointOnTriangle.h>
#include <Rodin/Variational/H1.h>
#include <Rodin/Variational/P1.h>

using namespace Rodin;
using namespace Rodin::Test;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Verifies sanity test reference triangle for geometry parametric
   * transformation by checking tolerance-based numerical results.
   */
  TEST(Rodin_Geometry_ParametricTransformation, SanityTest_ReferenceTriangle)
  {
    constexpr const size_t sdim = 2;
    constexpr const size_t n = 3;

    PointCloud pm(sdim, n);
    pm(0, 0) = 0;
    pm(0, 1) = 1;
    pm(0, 2) = 0;
    pm(1, 0) = 0;
    pm(1, 1) = 0;
    pm(1, 2) = 1;

    Variational::RealP1Element fe(Polytope::Type::Triangle);
    ParametricTransformation trans(pm, fe);

    Math::SpatialPoint res;
    Math::SpatialPoint inv;

    for (size_t i = 0; i < 3; i++)
    {
      trans.transform(res, pm.col(i));
      trans.inverse(inv, res);
      EXPECT_NEAR((res - pm.col(i)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR((inv - pm.col(i)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  /**
   * @brief Verifies sanity test triangle 1 for geometry parametric
   * transformation by checking tolerance-based numerical results.
   */
  TEST(Rodin_Geometry_ParametricTransformation, SanityTest_Triangle_1)
  {
    constexpr const size_t rdim = 2;
    constexpr const size_t sdim = 2;
    constexpr const size_t n = 3;

    PointCloud pm(sdim, n);
    pm(0, 0) = -1;
    pm(0, 1) = 1;
    pm(0, 2) = 0;
    pm(1, 0) = -1;
    pm(1, 1) = 1;
    pm(1, 2) = 1;

    Variational::RealP1Element fe(Polytope::Type::Triangle);
    ParametricTransformation trans(pm, fe);

    Math::SpatialPoint rc(rdim);
    Math::SpatialPoint res;
    Math::SpatialPoint inv;

    {
      rc[0] = 0;
      rc[1] = 0;
      trans.transform(res, rc);
      trans.inverse(inv, res);
      EXPECT_NEAR((res - pm.col(0)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      rc[0] = 1;
      rc[1] = 0;
      trans.transform(res, rc);
      trans.inverse(inv, res);
      EXPECT_NEAR((res - pm.col(1)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      rc[0] = 0;
      rc[1] = 1;
      trans.transform(res, rc);
      trans.inverse(inv, res);
      EXPECT_NEAR((res - pm.col(2)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::SpatialPoint rc(rdim);

      rc[0] = 1.0 / 3.0;
      rc[1] = 1.0 / 3.0;

      Math::SpatialPoint pc(sdim);
      pc[0] = 0;
      pc[1] = (1.0 / 3.0);

      trans.transform(res, rc);
      EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);

      trans.inverse(inv, res);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::SpatialPoint rc(rdim);
      rc[0] = 0.5;
      rc[1] = 0;

      Math::SpatialPoint pc(sdim);
      pc[0] = 0;
      pc[1] = 0;

      trans.transform(res, rc);
      EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);

      trans.inverse(inv, res);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::SpatialPoint rc(rdim);
      rc[0] = 0.5;
      rc[1] = 0.5;

      Math::SpatialPoint pc(sdim);
      pc[0] = 0.5;
      pc[1] = 1;

      trans.transform(res, rc);
      EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);

      trans.inverse(inv, res);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::SpatialPoint rc(rdim);
      rc[0] = 0.5;
      rc[1] = 0.5;

      Math::SpatialPoint pc(sdim);
      pc[0] = 0.5;
      pc[1] = 1;

      trans.transform(res, rc);
      EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);

      trans.inverse(inv, res);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  template <size_t K>
  static PointCloud makeQuarterCircleTrianglePointCloud()
  {
    constexpr Real Pi = 3.14159265358979323846;
    Variational::RealH1Element<K> fe(Polytope::Type::Triangle);
    PointCloud pm(2, fe.getCount());
    const Math::SpatialPoint v0{1, 0};
    const Math::SpatialPoint v1{0, 1};
    const Math::SpatialPoint v2{0, 0};
    for (size_t i = 0; i < fe.getCount(); ++i)
    {
      const auto& rc = fe.getNode(i);
      Math::SpatialPoint x = (Real(1) - rc[0] - rc[1]) * v0 + rc[0] * v1 + rc[1] * v2;
      if (std::abs(rc[1]) <= Real(1e-12))
      {
        const Real theta = (Pi / Real(2)) * rc[0];
        x = Math::SpatialPoint{std::cos(theta), std::sin(theta)};
      }
      pm(0, i) = x[0];
      pm(1, i) = x[1];
    }
    return pm;
  }

  /**
   * @brief Verifies P 2 curves interface edge for geometry parametric
   * transformation by checking tolerance-based numerical results.
   */
  TEST(Rodin_Geometry_ParametricTransformation, P2CurvesInterfaceEdge)
  {
    Variational::RealH1Element<2> fe(Polytope::Type::Triangle);
    ParametricTransformation trans(makeQuarterCircleTrianglePointCloud<2>(), fe);

    Math::SpatialPoint x;
    trans.transform(x, Math::SpatialPoint{Real(0.5), Real(0)});
    EXPECT_NEAR(x.norm(), Real(1), Real(1e-12));

    const Math::SpatialPoint linearMidpoint{Real(0.5), Real(0.5)};
    EXPECT_GT(std::abs(linearMidpoint.norm() - Real(1)), std::abs(x.norm() - Real(1)));
  }

  /**
   * @brief Verifies P 3 curved edge improves fit for geometry parametric
   * transformation.
   */
  TEST(Rodin_Geometry_ParametricTransformation, P3CurvedEdgeImprovesFit)
  {
    Variational::RealH1Element<3> fe(Polytope::Type::Triangle);
    ParametricTransformation trans(makeQuarterCircleTrianglePointCloud<3>(), fe);

    Real maxError = 0;
    for (Real s : {Real(0.25), Real(0.5), Real(0.75)})
    {
      Math::SpatialPoint x;
      trans.transform(x, Math::SpatialPoint{s, Real(0)});
      maxError = std::max(maxError, std::abs(x.norm() - Real(1)));
    }

    const Math::SpatialPoint linearMidpoint{Real(0.5), Real(0.5)};
    EXPECT_LT(maxError, std::abs(linearMidpoint.norm() - Real(1)));
  }

  /**
   * @brief Verifies curved triangle jacobian stays positive for geometry
   * parametric transformation.
   */
  TEST(Rodin_Geometry_ParametricTransformation, CurvedTriangleJacobianStaysPositive)
  {
    Variational::RealH1Element<2> fe(Polytope::Type::Triangle);
    ParametricTransformation trans(makeQuarterCircleTrianglePointCloud<2>(), fe);

    for (const Math::SpatialPoint& rc :
      {Math::SpatialPoint{Real(1) / Real(3), Real(1) / Real(3)},
        Math::SpatialPoint{Real(0.25), Real(0.25)},
        Math::SpatialPoint{Real(0.5), Real(0.25)},
        Math::SpatialPoint{Real(0.25), Real(0.5)}})
    {
      Math::SpatialMatrix<Real> J;
      trans.jacobian(J, rc);
      EXPECT_GT(J.determinant(), Real(0));
    }
  }
  namespace
  {
    // A pure basis evaluation with an observable call count checks collaboration
    // with generic scalar elements, including embedded transformations.
    template <size_t K>
    struct CountedElement
    {
      using RangeType = Real;
      using Element = Variational::RealH1Element<K>;
      Element element;
      size_t* evaluations;

      auto getGeometry() const
      {
        return element.getGeometry();
      }
      size_t getCount() const
      {
        return element.getCount();
      }
      size_t getOrder() const
      {
        return element.getOrder();
      }

      struct Basis
      {
        typename Element::BasisFunction basis;
        size_t* evaluations;
        Real operator()(const Math::SpatialPoint& r) const
        {
          return basis(r);
        }
        template <size_t Order>
        struct Derivative
        {
          typename Element::BasisFunction::template DerivativeFunction<Order>
            derivative;
          size_t* evaluations;
          Real operator()(const Math::SpatialPoint& r) const
          {
            ++*evaluations;
            return derivative(r);
          }
        };
        template <size_t Order>
        auto getDerivative(size_t axis) const
        {
          return Derivative<Order>{
            basis.template getDerivative<Order>(axis), evaluations};
        }
      };
      Basis getBasis(size_t local) const
      {
        return {element.getBasis(local), evaluations};
      }
    };

    template <size_t K>
    void checkEvaluation(Polytope::Type geometry)
    {
      const Polytope::Traits traits(geometry);
      Variational::RealH1Element<K> fe(geometry);
      for (size_t physicalDimension = std::max(size_t(1), traits.getDimension());
           physicalDimension <= 3; ++physicalDimension)
      {
        PointCloud nodes(physicalDimension, fe.getCount());
        // Alternating coefficients exercise cancellation and every basis term.
        for (size_t a = 0; a < fe.getCount(); ++a)
        {
          for (size_t j = 0; j < physicalDimension; ++j)
            nodes(j, a) = Real(a + j + 1) / Real(fe.getCount()) * (a % 2 ? -1 : 1);
        }
        ParametricTransformation transformation(nodes, fe);
        std::vector<Math::SpatialPoint> references{traits.getCentroid()};
        // A Point has a zero-dimensional reference domain. Its centroid is
        // the sole reference point; getVertex(0) has one stored coordinate.
        if (traits.getDimension() > 0)
        {
          for (size_t v = 0; v < traits.getVertexCount(); ++v)
            references.push_back(traits.getVertex(v));
        }
        for (const auto& reference : references)
        {
          ASSERT_EQ(static_cast<size_t>(reference.size()), traits.getDimension());
          Math::SpatialPoint actualPoint;
          Math::SpatialMatrix<Real> actualJacobian;
          transformation.transform(actualPoint, reference);
          transformation.jacobian(actualJacobian, reference);
          Math::SpatialPoint expectedPoint(physicalDimension);
          Math::SpatialMatrix<Real> expectedJacobian(
            physicalDimension, traits.getDimension());
          expectedPoint.setZero();
          expectedJacobian.setZero();
          // Preserve the original evaluation and accumulation order as oracle.
          for (size_t a = 0; a < fe.getCount(); ++a)
          {
            expectedPoint += nodes[a] * fe.getBasis(a)(reference);
            for (size_t i = 0; i < traits.getDimension(); ++i)
            {
              for (size_t j = 0; j < physicalDimension; ++j)
              {
                expectedJacobian(j, i) +=
                  nodes(j, a) * fe.getBasis(a).template getDerivative<1>(i)(reference);
              }
            }
          }
          for (size_t j = 0; j < physicalDimension; ++j)
          {
            EXPECT_EQ(actualPoint[j], expectedPoint[j]);
            for (size_t i = 0; i < traits.getDimension(); ++i)
              EXPECT_EQ(actualJacobian(j, i), expectedJacobian(j, i));
          }
        }
      }
    }
  } // namespace

  /**
   * @brief Checks exact agreement with scalar evaluation on all geometries and
   * embeddings.
   */
  TEST(
    Rodin_Geometry_ParametricTransformation, AllGeometryEvaluationMatchesScalarBaseline)
  {
    for (const auto geometry :
      {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      checkEvaluation<1>(geometry);
      checkEvaluation<2>(geometry);
      checkEvaluation<3>(geometry);
      checkEvaluation<4>(geometry);
    }
  }

  /**
   * @brief Evaluates a basis derivative once, independently of physical
   * dimension.
   */
  TEST(Rodin_Geometry_ParametricTransformation, JacobianEvaluatesDerivativeOnce)
  {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      auto check = [&]<size_t K>() {
        SCOPED_TRACE(K);
        size_t evaluations = 0;
        CountedElement<K> fe{Variational::RealH1Element<K>(geometry), &evaluations};
        PointCloud nodes(3, fe.getCount());
        nodes.setZero();
        ParametricTransformation transformation(std::move(nodes), fe);
        Math::SpatialMatrix<Real> J;
        transformation.jacobian(J, Polytope::Traits(geometry).getCentroid());
        EXPECT_EQ(evaluations, fe.getCount() * Polytope::Traits(geometry).getDimension());
        evaluations = 0;
        const auto& formula = QF::PolytopeQuadratureFormula::get(4, geometry);
        std::vector<Math::SpatialMatrix<Real>> matrices;
        transformation.jacobian(matrices, formula);
        EXPECT_EQ(matrices.size(), formula.getSize());
        EXPECT_EQ(evaluations,
          formula.getSize() * fe.getCount() * Polytope::Traits(geometry).getDimension());
      };
      SCOPED_TRACE(static_cast<int>(geometry));
      check.template operator()<1>();
      check.template operator()<2>();
      check.template operator()<3>();
      check.template operator()<4>();
    }
  }

  /**
   * @brief Exercises quadrature Jacobians independently of field and backend.
   *
   * ## Architecture
   *
   * Bulk evaluation is compared with direct basis evaluation using distinct
   * control clouds. Formula-read and basis-read counters independently certify
   * that reference tables, rather than physical Jacobians, are reused across
   * transformations. Logical lifetime changes invalidate reference reuse.
   */
  class QuadratureTransformationTest : public ::testing::Test
  {
    protected:
      /**
       * @brief Compares bulk evaluation with the pointwise transformation.
       * @tparam K Geometry element degree.
       * @param geometry Reference polytope type.
       */
      template <size_t K, class Element = Variational::RealH1Element<K>>
      void checkEvaluation(Polytope::Type geometry)
      {
        SCOPED_TRACE(K);
        Element element(geometry);
        const size_t rdim = Polytope::Traits(geometry).getDimension();
        for (const size_t order : {size_t(8), size_t(16)})
        {
          SCOPED_TRACE(order);
          const auto& formula = QF::PolytopeQuadratureFormula::get(order, geometry);
          for (size_t pdim = std::max(size_t(1), rdim); pdim <= 3; ++pdim)
          {
            for (size_t cell = 0; cell < 2; ++cell)
            {
              PointCloud controls(pdim, element.getCount());
              for (size_t a = 0; a < element.getCount(); ++a)
              {
                for (size_t j = 0; j < pdim; ++j)
                {
                  controls(j, a) = Real((a + 1) * (j + 1) + cell) /
                    Real(element.getCount() + cell + 1) * (a % 2 ? -1 : 1);
                }
              }
              ParametricTransformation transformation(controls, element);
              const PolytopeTransformation& polymorphic = transformation;
              std::vector<Math::SpatialMatrix<Real>> matrices;
              polymorphic.jacobian(matrices, formula);
              ASSERT_EQ(matrices.size(), formula.getSize());
              for (size_t qp = 0; qp < matrices.size(); ++qp)
              {
                Math::SpatialMatrix<Real> expected;
                transformation.jacobian(expected, formula.getPoint(qp));
                ASSERT_EQ(matrices[qp].rows(), expected.rows());
                ASSERT_EQ(matrices[qp].cols(), expected.cols());
                for (size_t j = 0; j < pdim; ++j)
                {
                  for (size_t i = 0; i < rdim; ++i)
                    EXPECT_EQ(matrices[qp](j, i), expected(j, i));
                }
              }
            }
          }
        }
      }

      class ObservedFormula final : public QF::QuadratureFormulaBase
      {
        public:
          /**
           * @brief Constructs an observed one-point formula.
           * @param point Reference coordinate returned by the formula.
           */
          ObservedFormula(Polytope::Type geometry, size_t sample)
            : ObservedFormula(sample ? Polytope::Traits(geometry).getVertex(0)
                                     : Polytope::Traits(geometry).getCentroid())
          {}

          explicit ObservedFormula(const Math::SpatialPoint& point)
            : m_point(point),
              m_reads(0)
          {}

          /**
           * @brief Copies the point with a new logical formula identity.
           * @param other Formula whose point is copied.
           */
          ObservedFormula(const ObservedFormula& other)
            : QuadratureFormulaBase(other),
              m_point(other.m_point),
              m_reads(0)
          {}

          /**
           * @brief Assigns a point and invalidates the old logical identity.
           * @param other Formula whose point is copied.
           * @returns This formula after assignment.
           */
          ObservedFormula& operator=(const ObservedFormula& other)
          {
            QuadratureFormulaBase::operator=(other);
            m_point = other.m_point;
            m_reads = 0;
            return *this;
          }

          /**
           * @brief Returns the single reference sample.
           * @returns The formula size, one.
           */
          size_t getSize() const override
          {
            return 1;
          }

          /**
           * @brief Returns the dummy unit weight used by the cache probe.
           * @param i Sample index; unused because this probe does not integrate.
           * @returns Unit weight.
           */
          Real getWeight([[maybe_unused]] size_t i) const override
          {
            return 1;
          }

          /**
           * @brief Reads and counts the single reference sample.
           * @param i Sample index, which must be zero.
           * @returns The stored reference coordinates.
           */
          const Math::SpatialPoint& getPoint(size_t i) const override
          {
            assert(i == 0);
            ++m_reads;
            return m_point;
          }

          /**
           * @brief Returns the number of reference-point reads.
           * @returns Number of calls to getPoint().
           */
          size_t getReads() const
          {
            return m_reads;
          }

          /**
           * @brief Copies the observed formula polymorphically.
           * @returns A newly owned formula copy.
           */
          ObservedFormula* copy() const noexcept
          {
            return new ObservedFormula(*this);
          }

        private:
          Math::SpatialPoint m_point;
          mutable size_t m_reads;
      };

      template <class Element = Variational::RealH1Element<2>>
      struct ObservedElement
      {
        using RangeType = Real;
        Element element;
        size_t* basisReads;
        size_t* tableReads;
        size_t* countReads;

        /**
         * @brief Observes an H1 geometry element's evaluation entry points.
         * @param geometry Reference geometry type.
         * @param basis Counter for direct basis access.
         * @param table Counter for reference-table access.
         * @param count Counter for basis-count queries.
         */
        ObservedElement(
          Polytope::Type geometry, size_t& basis, size_t& table, size_t& count)
          : element(geometry),
            basisReads(&basis),
            tableReads(&table),
            countReads(&count)
        {}

        /**
         * @brief Gets the underlying reference geometry.
         * @returns The element geometry type.
         */
        auto getGeometry() const
        {
          return element.getGeometry();
        }

        /**
         * @brief Gets the underlying basis count.
         * @returns Number of scalar basis functions.
         */
        size_t getCount() const
        {
          ++*countReads;
          return element.getCount();
        }

        /**
         * @brief Gets the geometry element order.
         * @returns The element order.
         */
        size_t getOrder() const
        {
          return element.getOrder();
        }

        /**
         * @brief Counts direct scalar basis access.
         * @param i Local scalar basis index.
         * @returns The underlying basis function.
         */
        decltype(auto) getBasis(size_t i) const
        {
          ++*basisReads;
          return element.getBasis(i);
        }

        /**
         * @brief Counts access to the existing shared reference table.
         * @param formula Formula defining the reference sample set.
         * @returns The underlying element tabulation.
         */
        const auto& getTabulation(const QF::QuadratureFormula& formula) const
        {
          ++*tableReads;
          return element.getTabulation(formula);
        }

        /**
         * @brief Observes lookup of an existing reference table.
         * @param identity Logical formula identity, including zero.
         * @returns Borrowed table when present on this thread.
         */
        auto findTabulation(const QF::QuadratureFormula& formula) const
        {
          ++*tableReads;
          return element.findTabulation(formula);
        }
      };

      /**
       * @brief Certifies bounded pointwise metadata access and exact evaluation.
       *
       * The basis count is invariant during one map evaluation. Each operation
       * must query it once, independently of the number of basis functions.
       * Exact comparisons retain the reference accumulation order for every
       * admissible physical embedding, including zero-dimensional points.
       *
       * @param geometry Reference polytope type.
       */
      template <class Element>
      void checkPointwiseCount(Polytope::Type geometry)
      {
        const Polytope::Traits traits(geometry);
        size_t basisReads = 0, tableReads = 0, countReads = 0;
        ObservedElement<Element> element(geometry, basisReads, tableReads, countReads);
        const size_t count = element.element.getCount();
        for (size_t dimension = traits.getDimension(); dimension <= 3; ++dimension)
        {
          SCOPED_TRACE(dimension);
          PointCloud controls(dimension, count);
          for (size_t a = 0; a < count; ++a)
          {
            for (size_t j = 0; j < dimension; ++j)
              controls(j, a) = Real((a + 1) * (j + 1)) / count;
          }
          ParametricTransformation observed(controls, element);
          for (const auto& reference : {traits.getCentroid(), traits.getVertex(0)})
          {
            Math::SpatialPoint actual, expected(dimension);
            expected.setZero();
            for (size_t a = 0; a < count; ++a)
              expected += controls[a] * element.element.getBasis(a)(reference);
            const size_t beforeTransform = countReads;
            observed.transform(actual, reference);
            EXPECT_EQ(countReads - beforeTransform, 1);
            ASSERT_EQ(actual.size(), expected.size());
            for (size_t j = 0; j < dimension; ++j)
              EXPECT_EQ(actual(j), expected(j));

            Math::SpatialMatrix<Real> jacobian,
              expectedJacobian(dimension, traits.getDimension());
            expectedJacobian.setZero();
            for (size_t a = 0; a < count; ++a)
            {
              for (size_t i = 0; i < traits.getDimension(); ++i)
              {
                const Real derivative =
                  element.element.getBasis(a).template getDerivative<1>(i)(reference);
                for (size_t j = 0; j < dimension; ++j)
                  expectedJacobian(j, i) += controls(j, a) * derivative;
              }
            }
            const size_t beforeJacobian = countReads;
            observed.jacobian(jacobian, reference);
            EXPECT_EQ(countReads - beforeJacobian, 1);
            ASSERT_EQ(jacobian.rows(), expectedJacobian.rows());
            ASSERT_EQ(jacobian.cols(), expectedJacobian.cols());
            for (size_t j = 0; j < dimension; ++j)
            {
              for (size_t i = 0; i < traits.getDimension(); ++i)
                EXPECT_EQ(jacobian(j, i), expectedJacobian(j, i));
            }
          }
        }
        EXPECT_EQ(tableReads, 0);
      }

      /**
       * @brief Certifies lazy coordinates and reference-sample ownership.
       *
       * ## Invariants
       *
       * Logical reference reuse must perform no direct basis reads. A missing
       * table must evaluate the owned reference sample, never a surviving
       * formula pointer or a table slot reused for another formula. Both paths
       * preserve the scalar accumulation order exactly.
       *
       * @param geometry Positive-dimensional reference geometry.
       */
      template <class Element = Variational::RealH1Element<2>>
      void checkLazyCoordinates(Polytope::Type geometry)
      {
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        Array<size_t> shape(dimension);
        shape.setConstant(3);
        auto mesh = Mesh<Context::Local>::UniformGrid(geometry, shape);
        const auto cell = *mesh.getPolytope(dimension, 0);
        const auto otherCell = *mesh.getPolytope(dimension, 1);
        size_t basisReads = 0, tableReads = 0, countReads = 0;
        ObservedElement<Element> element(geometry, basisReads, tableReads, countReads);
        PointCloud controls(dimension, element.getCount());
        for (size_t a = 0; a < element.getCount(); ++a)
        {
          Point direct(cell, element.element.getNode(a));
          auto physical = direct.getPhysicalCoordinates();
          physical(dimension - 1) += physical(0) * physical(0) / 8;
          for (size_t j = 0; j < dimension; ++j)
            controls(j, a) = physical(j);
        }
        ParametricTransformation map(controls, element);
        ParametricTransformation baseline(controls, element.element);
        mesh.setPolytopeTransformation({dimension, cell.getIndex()}, map.copy());
        const PointCloud otherControls = 2 * controls;
        ParametricTransformation otherMap(otherControls, element);
        mesh.setPolytopeTransformation(
          {dimension, otherCell.getIndex()}, otherMap.copy());

        const auto reference = Polytope::Traits(geometry).getCentroid();
        Math::SpatialPoint expected;
        baseline.transform(expected, reference);
        ASSERT_GT(expected.squaredNorm(), 0);
        const QF::QuadratureFormula* formula =
          &QF::QuadratureFormula::get<ObservedFormula>(geometry, size_t{0});
        std::vector<Math::SpatialMatrix<Real>> matrices;
        map.jacobian(matrices, *formula);
        const auto* original = formula;
        ASSERT_TRUE(element.element.findTabulation(*formula));
        Point pending(matrices[0], cell, *formula, 0);
        Point retained(pending);
        Point moved(std::move(retained));
        EXPECT_EQ(basisReads, 0);
        EXPECT_EQ(tableReads, 1);
        const auto sameCoordinates = [&](const Math::SpatialPoint& actual) {
          ASSERT_EQ(actual.size(), expected.size());
          for (size_t j = 0; j < dimension; ++j)
            EXPECT_EQ(actual(j), expected(j));
        };
        const size_t countBeforeCoordinates = countReads;
        sameCoordinates(pending.getPhysicalCoordinates());
#ifdef NDEBUG
        EXPECT_EQ(countReads, countBeforeCoordinates);
#else
        EXPECT_EQ(countReads, countBeforeCoordinates + 1);
#endif
        EXPECT_EQ(basisReads, 0);
        EXPECT_EQ(tableReads, 2);
        static_cast<void>(pending.getPhysicalCoordinates());
        EXPECT_EQ(tableReads, 2);
        // Assignment at the same address establishes a new logical sample.
        formula = &QF::QuadratureFormula::get<ObservedFormula>(geometry, size_t{1});
        ASSERT_NE(formula, original);
        map.jacobian(matrices, *formula);
        formula = nullptr;
        sameCoordinates(moved.getPhysicalCoordinates());
        EXPECT_EQ(basisReads, 0);
        // A computed copy owns its physical coordinates and performs no lookup.
        Point computedCopy(moved);
        const auto reads = tableReads;
        sameCoordinates(computedCopy.getPhysicalCoordinates());
        EXPECT_EQ(tableReads, reads);

        formula = &QF::QuadratureFormula::get<ObservedFormula>(geometry, size_t{0});
        map.jacobian(matrices, *formula);
        Point evicted(matrices[0], cell, *formula, 0);
        Point transferred(evicted);
        const auto* evictedFormula = formula;
        formula = nullptr;
        // Nine distinct insertions exceed the existing eight-entry cache.
        for (size_t i = 0; i < 9; ++i)
        {
          const auto& replacement =
            QF::QuadratureFormula::get<ObservedFormula>(geometry, i + 2);
          static_cast<void>(element.element.getTabulation(replacement));
        }
        ASSERT_FALSE(element.element.findTabulation(*evictedFormula));
        const auto basisBefore = basisReads;
        sameCoordinates(evicted.getPhysicalCoordinates());
        EXPECT_EQ(basisReads, basisBefore + element.getCount());
        Math::SpatialPoint threadCoordinates;
        std::thread worker([&, point = std::move(transferred)]() mutable {
          threadCoordinates = point.getPhysicalCoordinates();
        });
        worker.join();
        sameCoordinates(threadCoordinates);
        EXPECT_EQ(basisReads, basisBefore + 2 * element.getCount());

        // Rebinding discards provenance even when the reference sample is equal.
        const auto lookupBefore = tableReads;
        pending.setPolytope(otherCell);
        const auto& rebound = pending.getPhysicalCoordinates();
        for (size_t j = 0; j < dimension; ++j)
          EXPECT_EQ(rebound(j), 2 * expected(j));
        EXPECT_EQ(tableReads, lookupBefore);
      }

      /**
       * @brief Certifies reference-only reuse and formula lifetime invalidation.
       * @param geometry Reference polytope type.
       */
      template <class Element = Variational::RealH1Element<2>>
      void checkReuse(Polytope::Type geometry)
      {
        const Polytope::Traits traits(geometry);
        const auto centroid = traits.getCentroid();
        const auto vertex = traits.getVertex(0);
        size_t basisReads = 0, tableReads = 0, countReads = 0;
        ObservedElement<Element> element(geometry, basisReads, tableReads, countReads);
        PointCloud first(3, element.getCount()), second(3, element.getCount());
        for (size_t a = 0; a < element.getCount(); ++a)
        {
          for (size_t j = 0; j < 3; ++j)
          {
            first(j, a) = Real((a + 1) * (j + 1)) / element.getCount();
            second(j, a) = 2 * first(j, a);
          }
        }
        ParametricTransformation firstMap(first, element), secondMap(second, element);
        const QF::QuadratureFormula* formula =
          &QF::QuadratureFormula::get<ObservedFormula>(geometry, size_t{0});

        const auto* initial = formula;
        std::vector<Math::SpatialMatrix<Real>> firstMatrices, secondMatrices;
        const size_t countBefore = countReads;
        firstMap.jacobian(firstMatrices, *formula);
        // Only the debug compatibility assertion may query the basis count.
#ifdef NDEBUG
        EXPECT_EQ(countReads, countBefore);
#else
        EXPECT_EQ(countReads, countBefore + 1);
#endif

        secondMap.jacobian(secondMatrices, *formula);

        EXPECT_EQ(tableReads, 2);
        EXPECT_EQ(basisReads, 0);
        for (size_t j = 0; j < 3; ++j)
        {
          for (size_t i = 0; i < traits.getDimension(); ++i)
            EXPECT_EQ(secondMatrices[0](j, i), 2 * firstMatrices[0](j, i));
        }
        formula = &QF::QuadratureFormula::get<ObservedFormula>(geometry, size_t{1});
        ASSERT_NE(formula, nullptr);
        ASSERT_NE(formula, initial);
        firstMap.jacobian(firstMatrices, *formula);

        const auto assignedMatrices = firstMatrices;
        const auto* assigned = formula;
        formula = &QF::QuadratureFormula::get<ObservedFormula>(geometry, size_t{0});
        ASSERT_NE(formula, nullptr);
        ASSERT_NE(formula, assigned);
        firstMap.jacobian(firstMatrices, *formula);

        EXPECT_EQ(tableReads, 4);
        EXPECT_EQ(basisReads, 0);
        ParametricTransformation baseline(first, element.element);
        Math::SpatialMatrix<Real> expectedCentroid, expectedVertex;
        baseline.jacobian(expectedCentroid, centroid);
        baseline.jacobian(expectedVertex, vertex);
        for (size_t j = 0; j < 3; ++j)
        {
          for (size_t i = 0; i < traits.getDimension(); ++i)
          {
            EXPECT_EQ(firstMatrices[0](j, i), expectedCentroid(j, i));
            EXPECT_EQ(assignedMatrices[0](j, i), expectedVertex(j, i));
          }
        }
      }

      /**
       * @brief Observes geometric evaluation without changing its arithmetic.
       *
       * Pointwise and bulk entry counters distinguish owned point-cache access
       * from repeated transformation evaluation. The wrapped Q2 map owns all
       * geometric coefficients; counters are borrowed for the test lifetime.
       */
      class ObservedTransformation final : public PolytopeTransformation
      {
        public:
          /**
           * @brief Constructs an observed parametric geometry map.
           * @param controls Physical geometry-node coordinates.
           * @param element Scalar Q2 reference element.
           * @param point Counter for pointwise Jacobian evaluation.
           * @param bulk Counter for quadrature Jacobian evaluation.
           */
          ObservedTransformation(const PointCloud& controls,
            const Variational::RealH1Element<2>& element, size_t& point, size_t& bulk)
            : PolytopeTransformation(
                Polytope::Traits(element.getGeometry()).getDimension(), controls.rows()),
              m_map(controls, element),
              m_point(&point),
              m_bulk(&bulk)
          {}

          /**
           * @brief Copies the geometry map and retains its observed counters.
           * @param other Transformation being copied.
           */
          ObservedTransformation(const ObservedTransformation& other)
            : PolytopeTransformation(other),
              m_map(other.m_map),
              m_point(other.m_point),
              m_bulk(other.m_bulk)
          {}

          /**
           * @brief Gets the underlying map order.
           * @returns The geometry order.
           */
          size_t getOrder() const override
          {
            return m_map.getOrder();
          }

          /**
           * @brief Evaluates the underlying physical map.
           * @param[out] physical Mapped physical coordinates.
           * @param[in] reference Reference coordinates.
           */
          void transform(Math::SpatialPoint& physical,
            const Math::SpatialPoint& reference) const override
          {
            m_map.transform(physical, reference);
          }

          /**
           * @brief Counts direct geometric Jacobian evaluation.
           * @param[out] matrix Geometric Jacobian.
           * @param[in] reference Reference coordinates.
           */
          void jacobian(Math::SpatialMatrix<Real>& matrix,
            const Math::SpatialPoint& reference) const override
          {
            ++*m_point;
            m_map.jacobian(matrix, reference);
          }

          /**
           * @brief Counts bulk geometric Jacobian evaluation.
           * @param[out] matrices Geometric Jacobians in formula order.
           * @param[in] formula Reference quadrature formula.
           */
          void jacobian(std::vector<Math::SpatialMatrix<Real>>& matrices,
            const QF::QuadratureFormula& formula) const override
          {
            ++*m_bulk;
            m_map.jacobian(matrices, formula);
          }

          /**
           * @brief Copies the observed transformation polymorphically.
           * @returns A newly owned transformation copy.
           */
          ObservedTransformation* copy() const noexcept override
          {
            return new ObservedTransformation(*this);
          }

        private:
          ParametricTransformation<Variational::RealH1Element<2>> m_map;
          size_t* m_point;
          size_t* m_bulk;
      };

      /**
       * @brief Certifies point ownership and invalidation after bulk evaluation.
       * @param geometry Positive-dimensional reference geometry type.
       */
      void checkPointOwnership(Polytope::Type geometry)
      {
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        Array<size_t> shape(dimension);
        shape.setConstant(3);
        auto mesh = Mesh<Context::Local>::UniformGrid(geometry, shape);
        const std::array<Polytope, 2> cells{
          *mesh.getPolytope(dimension, 0), *mesh.getPolytope(dimension, 1)};
        Variational::RealH1Element<2> element(geometry);
        std::array<PointCloud, 2> controls;
        size_t pointReads = 0, bulkReads = 0;
        // This shear has positive derivative in 1D and unit determinant in 2D/3D.
        constexpr Real ShearAmplitude = Real(1) / 8;
        for (size_t cell = 0; cell < cells.size(); ++cell)
        {
          controls[cell].resize(dimension, element.getCount());
          for (size_t a = 0; a < element.getCount(); ++a)
          {
            Point point(cells[cell], element.getNode(a));
            auto physical = point.getPhysicalCoordinates();
            physical(dimension - 1) += ShearAmplitude * physical(0) * physical(0);
            for (size_t j = 0; j < dimension; ++j)
              controls[cell](j, a) = physical(j);
          }
          mesh.setPolytopeTransformation({dimension, cells[cell].getIndex()},
            new ObservedTransformation(controls[cell], element, pointReads, bulkReads));
        }
        Optional<Point> retained;
        Math::SpatialPoint reference;
        {
          const auto& formula = QF::PolytopeQuadratureFormula::get(8, geometry);
          reference = formula.getPoint(0);
          PolytopeQuadrature quadrature(cells[0], formula);
          EXPECT_EQ(bulkReads, 1);
          for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
            static_cast<void>(quadrature.getPoint(qp).getJacobian());
          EXPECT_EQ(pointReads, 0);
          PolytopeQuadrature moved(std::move(quadrature));
          PolytopeQuadrature copied(moved);
          Point copy(copied.getPoint(0));
          retained.emplace(std::move(copy));
        }
        // The formula and both quadratures have been destroyed. The point owns
        // its coordinates and Jacobian, with no borrowed table or formula.
        const auto& jacobian = retained->getJacobian();
        EXPECT_EQ(pointReads, 0);
        EXPECT_EQ(bulkReads, 1);
        Point direct(cells[0], reference);
        const auto& expected = direct.getJacobian();
        for (size_t j = 0; j < dimension; ++j)
        {
          EXPECT_EQ(
            retained->getPhysicalCoordinates()(j), direct.getPhysicalCoordinates()(j));
          for (size_t i = 0; i < dimension; ++i)
          {
            EXPECT_EQ(jacobian(j, i), expected(j, i));
            EXPECT_EQ(
              retained->getJacobianInverse()(j, i), direct.getJacobianInverse()(j, i));
          }
        }
        EXPECT_EQ(retained->getJacobianDeterminant(), direct.getJacobianDeterminant());
        EXPECT_EQ(retained->getDistortion(), direct.getDistortion());
        EXPECT_EQ(pointReads, 1);
        retained->setPolytope(cells[1]);
        static_cast<void>(retained->getJacobian());
        EXPECT_EQ(pointReads, 2);
        Point rebound(cells[1], reference);
        for (size_t j = 0; j < dimension; ++j)
        {
          EXPECT_EQ(
            retained->getPhysicalCoordinates()(j), rebound.getPhysicalCoordinates()(j));
          for (size_t i = 0; i < dimension; ++i)
          {
            EXPECT_EQ(retained->getJacobian()(j, i), rebound.getJacobian()(j, i));
            EXPECT_EQ(
              retained->getJacobianInverse()(j, i), rebound.getJacobianInverse()(j, i));
          }
        }
        EXPECT_EQ(retained->getJacobianDeterminant(), rebound.getJacobianDeterminant());
        EXPECT_EQ(retained->getDistortion(), rebound.getDistortion());
        EXPECT_EQ(pointReads, 3);
        EXPECT_EQ(bulkReads, 1);
      }
  };

  /**
   * @brief Compares bulk and pointwise Jacobians exactly on every reference
   * geometry.
   */
  TEST_F(QuadratureTransformationTest, AllGeometriesOrdersAndEmbeddings)
  {
    for (const auto geometry :
      {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      checkEvaluation<1, Variational::P1Element<Real>>(geometry);
      checkEvaluation<1>(geometry);
      checkEvaluation<2>(geometry);
      checkEvaluation<3>(geometry);
      checkEvaluation<4>(geometry);
      checkEvaluation<5>(geometry);
      checkEvaluation<6>(geometry);
    }
  }

  /**
   * @brief Certifies retained P1 and H1 outputs across shape and formula changes.
   *
   * Packed input padding and all output storage are deliberately poisoned.
   * Only active physical/reference components may enter the contraction;
   * inactive output components remain zero after shrinking or expanding.
   * Distinct control clouds reuse the same reference table without reusing
   * physical Jacobians. Active entries agree exactly with pointwise evaluation.
   * Native P1 and H1 degrees one, two and six exercise both first-order and
   * high-order reference tables without assuming sample-invariant derivatives.
   */
  TEST_F(QuadratureTransformationTest, P1AndH1RetainedOutputsAndInactivePadding)
  {
    constexpr Real poison = 123;
    for (const auto geometry :
      {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const auto certify = [&]<class Element>() {
        const Element element(geometry);
        SCOPED_TRACE(element.getOrder());
        const size_t d = Polytope::Traits(geometry).getDimension();
        std::vector<Math::SpatialMatrix<Real>> matrices;
        for (const size_t order : {size_t(16), size_t(2), size_t(8)})
        {
          SCOPED_TRACE(order);
          const auto& formula = QF::PolytopeQuadratureFormula::get(order, geometry);
          for (const size_t s : {size_t(3), std::max(size_t(1), d), size_t(3)})
          {
            SCOPED_TRACE(s);
            PointCloud controls(s, element.getCount());
            for (size_t a = 0; a < element.getCount(); ++a)
            {
              for (size_t j = 0; j < 3; ++j)
              {
                controls.point3(a)[j] =
                  j < s ? Real((a + 1) * (j + 2) + order) / 17 : poison;
              }
            }
            const ParametricTransformation transformation(controls, element);
            for (auto& matrix : matrices)
              matrix.getData().setConstant(poison);
            transformation.jacobian(matrices, formula);
            ASSERT_EQ(matrices.size(), formula.getSize());
            for (size_t qp = 0; qp < matrices.size(); ++qp)
            {
              Math::SpatialMatrix<Real> expected;
              transformation.jacobian(expected, formula.getPoint(qp));
              ASSERT_EQ(matrices[qp].rows(), s);
              ASSERT_EQ(matrices[qp].cols(), d);
              for (size_t j = 0; j < 3; ++j)
              {
                for (size_t i = 0; i < 3; ++i)
                {
                  if (j < s && i < d)
                    EXPECT_EQ(matrices[qp](j, i), expected(j, i));
                  else
                    EXPECT_EQ(matrices[qp].getData()(j, i), 0);
                }
              }
            }
          }
        }
      };
      certify.template operator()<Variational::RealP1Element>();
      certify.template operator()<Variational::RealH1Element<1>>();
      certify.template operator()<Variational::RealH1Element<2>>();
      certify.template operator()<Variational::RealH1Element<6>>();
    }
  }

  /**
   * @brief Preserves planar, surface and volume basis-reduction order.
   *
   * For @f$p@f$ binary significand digits, controls cycle through
   * @f$2^{p+1}@f$, small non-dyadic values and @f$-2^{p+1}@f$.
   * Retained bulk outputs must match the original pointwise reduction exactly.
   * A reversed-basis oracle must differ in at least one entry, establishing
   * that the fixture can detect reassociation. Order-four quadrature supplies
   * deterministic interior samples; no integration-exactness claim is made.
   */
  TEST_F(QuadratureTransformationTest, SpatialCancellationPreservesBasisOrder)
  {
    constexpr size_t QuadratureOrder = 4;
    constexpr Real Poison = 42;
    const Real large = std::ldexp(Real(1), std::numeric_limits<Real>::digits + 1);
    size_t reversedDifferences = 0;
    for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
           Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
           Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t referenceDimension = Polytope::Traits(geometry).getDimension();
      for (size_t physicalDimension = referenceDimension; physicalDimension <= 3;
           ++physicalDimension)
      {
        SCOPED_TRACE(physicalDimension);
        const auto certify = [&]<class Element>() {
          const Element element(geometry);
          SCOPED_TRACE(element.getOrder());
          PointCloud controls(physicalDimension, element.getCount());
          for (size_t a = 0; a < element.getCount(); ++a)
          {
            for (size_t j = 0; j < physicalDimension; ++j)
            {
              const size_t group = (a + j) % 3;
              controls(j, a) = group == 0 ? large
                : group == 1              ? Real(a + j + 1) / 7
                                          : -large;
            }
          }
          const auto& formula =
            QF::PolytopeQuadratureFormula::get(QuadratureOrder, geometry);
          const ParametricTransformation transformation(controls, element);
          std::vector<Math::SpatialMatrix<Real>> matrices(formula.getSize());
          for (auto& matrix : matrices)
            matrix.getData().setConstant(Poison);
          transformation.jacobian(matrices, formula);
          ASSERT_EQ(matrices.size(), formula.getSize());
          const auto& table = element.getTabulation(formula);
          for (size_t qp = 0; qp < formula.getSize(); ++qp)
          {
            SCOPED_TRACE(qp);
            Math::SpatialMatrix<Real> expected;
            transformation.jacobian(expected, formula.getPoint(qp));
            Math::SpatialMatrix<Real> reversed(physicalDimension, referenceDimension);
            for (size_t a = element.getCount(); a > 0; --a)
            {
              const auto gradient = table.getGradient(qp, a - 1);
              for (size_t i = 0; i < referenceDimension; ++i)
              {
                for (size_t j = 0; j < physicalDimension; ++j)
                  reversed(j, i) += controls(j, a - 1) * gradient[i];
              }
            }
            ASSERT_EQ(matrices[qp].rows(), physicalDimension);
            ASSERT_EQ(matrices[qp].cols(), referenceDimension);
            for (size_t j = 0; j < physicalDimension; ++j)
            {
              for (size_t i = 0; i < referenceDimension; ++i)
              {
                EXPECT_EQ(matrices[qp](j, i), expected(j, i));
                reversedDifferences += reversed(j, i) != expected(j, i);
              }
            }
          }
        };
        certify.template operator()<Variational::RealP1Element>();
        certify.template operator()<Variational::RealH1Element<2>>();
        certify.template operator()<Variational::RealH1Element<6>>();
      }
    }
    EXPECT_GT(reversedDifferences, 0);
  }

  /**
   * @brief Checks pointwise metadata access for P1 and H1 orders one to six.
   *
   * All eight reference geometries and physical dimensions up to three are
   * evaluated at the centroid and a vertex with exact arithmetic oracles.
   */
  TEST_F(QuadratureTransformationTest, PointwiseCountIsIndependentOfBasisSize)
  {
    for (const auto geometry :
      {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      checkPointwiseCount<Variational::P1Element<Real>>(geometry);
      checkPointwiseCount<Variational::RealH1Element<1>>(geometry);
      checkPointwiseCount<Variational::RealH1Element<2>>(geometry);
      checkPointwiseCount<Variational::RealH1Element<3>>(geometry);
      checkPointwiseCount<Variational::RealH1Element<4>>(geometry);
      checkPointwiseCount<Variational::RealH1Element<5>>(geometry);
      checkPointwiseCount<Variational::RealH1Element<6>>(geometry);
    }
  }

  /// @brief Certifies logical formula identity and cross-cell reference reuse.
  TEST_F(QuadratureTransformationTest, ReferenceReuseAndFormulaLifetime)
  {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      checkReuse(geometry);
    }
  }

  /**
   * @brief Verifies owned mapped points survive source lifetimes and invalidate
   * on rebinding.
   */
  TEST_F(QuadratureTransformationTest, PointOwnershipCopyMoveAndRebinding)
  {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      checkPointOwnership(geometry);
    }
  }

  /**
   * @brief Verifies lazy reference reuse, eviction, lifetimes and thread
   * transfer.
   */
  TEST_F(QuadratureTransformationTest, LazyCoordinatesOwnReferenceSamples)
  {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      checkLazyCoordinates(geometry);
    }
  }

  /**
   * @brief Certifies P1 reference reuse, invalidation and safe lazy fallbacks.
   *
   * Integer access counts distinguish shared reference tables from pointwise
   * evaluation. Exact Jacobian and coordinate comparisons cover all supported
   * positive-dimensional geometries, including rational pyramids, as well as
   * formula destruction, eviction, rebinding and cross-thread point transfer.
   */
  TEST_F(QuadratureTransformationTest, P1ReferenceReuseAndOwnedSamples)
  {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      checkReuse<Variational::P1Element<Real>>(geometry);
      checkLazyCoordinates<Variational::P1Element<Real>>(geometry);
    }
  }
} // namespace Rodin::Tests::Unit
