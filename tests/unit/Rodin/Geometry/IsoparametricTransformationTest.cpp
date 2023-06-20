#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational/P1.h>
#include <Rodin/Test/Random/RandomPointOnTriangle.h>

using namespace Rodin;
using namespace Rodin::Test;
using namespace Rodin::Geometry;

TEST(Rodin_Geometry_IsoparametricTransformation, SanityTest_ReferenceTriangle)
{
  constexpr const size_t sdim = 2;
  constexpr const size_t n = 3;

  Math::Matrix pm(sdim, n);
  pm << 0, 1, 0,
        0, 0, 1;

  Variational::ScalarP1Element fe(Polytope::Geometry::Triangle);
  IsoparametricTransformation trans(pm, fe);

  {
    const auto res = trans.transform(pm.col(0));
    EXPECT_NEAR((res - pm.col(0)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  {
    const auto res = trans.transform(pm.col(1));
    EXPECT_NEAR((res - pm.col(1)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  {
    const auto res = trans.transform(pm.col(2));
    EXPECT_NEAR((res - pm.col(2)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }
}

TEST(Rodin_Geometry_IsoparametricTransformation, SanityTest_Triangle_1)
{
  constexpr const size_t rdim = 2;
  constexpr const size_t sdim = 2;
  constexpr const size_t n = 3;

  Math::Matrix pm(sdim, n);
  pm << -1, 1, 0,
        -1, 1, 1;

  Variational::ScalarP1Element fe(Polytope::Geometry::Triangle);
  IsoparametricTransformation trans(pm, fe);

  {
    Math::Vector r(rdim);
    r << 0, 0;
    const auto res = trans.transform(r);
    EXPECT_NEAR((res - pm.col(0)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  {
    Math::Vector r(rdim);
    r << 1, 0;
    const auto res = trans.transform(r);
    EXPECT_NEAR((res - pm.col(1)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  {
    Math::Vector r(rdim);
    r << 0, 1;
    const auto res = trans.transform(r);
    EXPECT_NEAR((res - pm.col(2)).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  {
    Math::Vector r(rdim);
    r << (1.0 / 3.0), (1.0 / 3.0);

    Math::Vector p(sdim);
    p << 0, (1.0 / 3.0);

    const auto res = trans.transform(r);
    EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  {
    Math::Vector r(rdim);
    r << 0.5, 0;

    Math::Vector p(sdim);
    p << 0, 0;
    const auto res = trans.transform(r);
    EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  {
    Math::Vector r(rdim);
    r << 0.5, 0.5;

    Math::Vector p(sdim);
    p << 0.5, 1;

    const auto res = trans.transform(r);
    EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  {
    Math::Vector r(rdim);
    r << 0.5, 0.5;

    Math::Vector p(sdim);
    p << 0.5, 1;

    const auto res = trans.transform(r);
    EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }
}

TEST(Rodin_Geometry_IsoparametricTransformation, PermutationTest_Triangle_1)
{
  constexpr const size_t sdim = 2;
  constexpr const size_t n = 3;

  Math::Matrix pm1(sdim, n);
  pm1 << -1, 1, 0,
         -1, 1, 1;

  Math::Matrix pm2(sdim, n);
  pm2 << 1, -1, 0,
         1, -1, 1;

  Variational::ScalarP1Element fe(Polytope::Geometry::Triangle);
  IsoparametricTransformation trans1(pm1, fe);
  IsoparametricTransformation trans2(pm2, fe);

  Random::PointOnTriangle gen;
  const auto rc = gen();
  const auto p1 = trans1.transform(rc);
  const auto p2 = trans2.transform(rc);

  std::cout << "-----------\n";
  std::cout << rc << std::endl;
  std::cout << "--- p1: \n";
  std::cout << p1 << std::endl;
  std::cout << "--- p2: \n";
  std::cout << p2 << std::endl;

  EXPECT_NEAR((p1 - p2).norm(), 0.0, RODIN_FUZZY_CONSTANT);
}

