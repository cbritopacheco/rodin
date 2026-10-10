#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include <Rodin/Geometry.h>
#include <Rodin/Variational/H1.h>
#include <Rodin/Variational/P1.h>
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
      };
      SCOPED_TRACE(static_cast<int>(geometry));
      check.template operator()<1>();
      check.template operator()<2>();
      check.template operator()<3>();
      check.template operator()<4>();
    }
  }
}
