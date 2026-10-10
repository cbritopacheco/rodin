#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/Variational/H1.h>
#include <Rodin/Variational/P1.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Test/Random/RandomPointOnTriangle.h>

using namespace Rodin;
using namespace Rodin::Test;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Verifies sanity test reference triangle for geometry parametric transformation by checking tolerance-based numerical results.
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

  /// @brief Verifies sanity test triangle 1 for geometry parametric transformation by checking tolerance-based numerical results.
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

  /// @brief Verifies P 2 curves interface edge for geometry parametric transformation by checking tolerance-based numerical results.
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

  /// @brief Verifies P 3 curved edge improves fit for geometry parametric transformation.
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

  /// @brief Verifies curved triangle jacobian stays positive for geometry parametric transformation.
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
          for (size_t v = 0; v < traits.getVertexCount(); ++v)
            references.push_back(traits.getVertex(v));
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
  }

  /// @brief Checks exact agreement with scalar evaluation on all geometries and embeddings.
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

  /// @brief Evaluates a basis derivative once, independently of physical dimension.
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
      template <size_t K>
      void checkEvaluation(Polytope::Type geometry)
      {
        SCOPED_TRACE(K);
        Variational::RealH1Element<K> element(geometry);
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
          ObservedFormula* copy() const noexcept override
          {
            return new ObservedFormula(*this);
          }

        private:
          Math::SpatialPoint m_point;
          mutable size_t m_reads;
      };

      struct ObservedElement
      {
          using RangeType = Real;
          Variational::RealH1Element<2> element;
          size_t* basisReads;
          size_t* tableReads;
          /**
         * @brief Observes an H1 geometry element's evaluation entry points.
         * @param geometry Reference geometry type.
         * @param basis Counter for direct basis access.
         * @param table Counter for reference-table access.
         */
          ObservedElement(Polytope::Type geometry, size_t& basis, size_t& table)
            : element(geometry),
              basisReads(&basis),
              tableReads(&table)
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
          const auto& getBasis(size_t i) const
          {
            ++*basisReads;
            return element.getBasis(i);
          }
          /**
         * @brief Counts access to the existing shared reference table.
         * @param formula Formula defining the reference sample set.
         * @returns The underlying element tabulation.
         */
          const auto& getTabulation(const QF::QuadratureFormulaBase& formula) const
          {
            ++*tableReads;
            return element.getTabulation(formula);
          }
      };

      /**
       * @brief Certifies reference-only reuse and formula lifetime invalidation.
       * @param geometry Reference polytope type.
       */
      void checkReuse(Polytope::Type geometry)
      {
        const Polytope::Traits traits(geometry);
        const auto centroid = traits.getCentroid();
        const auto vertex = traits.getVertex(0);
        size_t basisReads = 0, tableReads = 0;
        ObservedElement element(geometry, basisReads, tableReads);
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
        Optional<ObservedFormula> formula;
        formula.emplace(centroid);
        const auto* address = &*formula;
        const auto initialIdentity = formula->getCacheIdentity();
        std::vector<Math::SpatialMatrix<Real>> firstMatrices, secondMatrices;
        firstMap.jacobian(firstMatrices, *formula);
        ASSERT_EQ(formula->getReads(), 1);
        secondMap.jacobian(secondMatrices, *formula);
        EXPECT_EQ(formula->getReads(), 1);
        EXPECT_EQ(tableReads, 2);
        EXPECT_EQ(basisReads, 0);
        for (size_t j = 0; j < 3; ++j)
        {
          for (size_t i = 0; i < traits.getDimension(); ++i)
            EXPECT_EQ(secondMatrices[0](j, i), 2 * firstMatrices[0](j, i));
        }
        *formula = ObservedFormula(vertex);
        ASSERT_EQ(&*formula, address);
        ASSERT_NE(formula->getCacheIdentity(), initialIdentity);
        firstMap.jacobian(firstMatrices, *formula);
        EXPECT_EQ(formula->getReads(), 1);
        const auto assignedMatrices = firstMatrices;
        const auto assignedIdentity = formula->getCacheIdentity();
        formula.emplace(centroid);
        ASSERT_EQ(&*formula, address);
        ASSERT_NE(formula->getCacheIdentity(), assignedIdentity);
        firstMap.jacobian(firstMatrices, *formula);
        EXPECT_EQ(formula->getReads(), 1);
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
            const QF::QuadratureFormulaBase& formula) const override
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
          QF::PolytopeQuadratureFormula formula(8, geometry);
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

  /// @brief Compares bulk and pointwise Jacobians exactly on every reference geometry.
  TEST_F(QuadratureTransformationTest, AllGeometriesOrdersAndEmbeddings)
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
      checkEvaluation<5>(geometry);
      checkEvaluation<6>(geometry);
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

  /// @brief Verifies owned mapped points survive source lifetimes and invalidate on rebinding.
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
}
