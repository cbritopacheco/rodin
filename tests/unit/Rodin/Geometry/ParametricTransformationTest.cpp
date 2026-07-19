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
}
