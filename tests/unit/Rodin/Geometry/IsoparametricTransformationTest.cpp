#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational/P1.h>
#include <Rodin/Test/Random/RandomPointOnTriangle.h>

using namespace Rodin;
using namespace Rodin::Test;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Geometry_IsoparametricTransformation, SanityTest_ReferenceTriangle)
  {
    constexpr const size_t sdim = 2;
    constexpr const size_t n = 3;

    Math::Matrix<Real> pm(sdim, n);
    pm << 0, 1, 0,
          0, 0, 1;

    Variational::RealP1Element fe(Polytope::Type::Triangle);
    IsoparametricTransformation trans(pm, fe);

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

  TEST(Rodin_Geometry_IsoparametricTransformation, SanityTest_Triangle_1)
  {
    constexpr const size_t rdim = 2;
    constexpr const size_t sdim = 2;
    constexpr const size_t n = 3;

    Math::Matrix<Real> pm(sdim, n);
    pm << -1, 1, 0,
          -1, 1, 1;

    Variational::RealP1Element fe(Polytope::Type::Triangle);
    IsoparametricTransformation trans(pm, fe);

    Math::SpatialPoint rc(rdim);
    Math::SpatialPoint res;
    Math::SpatialPoint inv;

    {
      rc << 0, 0;
      trans.transform(res, rc);
      trans.inverse(inv, res);
      EXPECT_NEAR((res - pm.col(0)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      rc << 1, 0;
      trans.transform(res, rc);
      trans.inverse(inv, res);
      EXPECT_NEAR((res - pm.col(1)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      rc << 0, 1;
      trans.transform(res, rc);
      trans.inverse(inv, res);
      EXPECT_NEAR((res - pm.col(2)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::SpatialPoint rc(rdim);
      rc << (1.0 / 3.0), (1.0 / 3.0);

      Math::SpatialPoint pc(sdim);
      pc << 0, (1.0 / 3.0);

      trans.transform(res, rc);
      EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);

      trans.inverse(inv, res);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::SpatialPoint rc(rdim);
      rc << 0.5, 0;

      Math::SpatialPoint pc(sdim);
      pc << 0, 0;

      trans.transform(res, rc);
      EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);

      trans.inverse(inv, res);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::SpatialPoint rc(rdim);
      rc << 0.5, 0.5;

      Math::SpatialPoint pc(sdim);
      pc << 0.5, 1;

      trans.transform(res, rc);
      EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);

      trans.inverse(inv, res);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::SpatialPoint rc(rdim);
      rc << 0.5, 0.5;

      Math::SpatialPoint pc(sdim);
      pc << 0.5, 1;

      trans.transform(res, rc);
      EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);

      trans.inverse(inv, res);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }
}
