#include <gtest/gtest.h>

#include <limits>
#include <type_traits>

#include <Eigen/Dense>

#include "Rodin/Math/Common.h"
#include "Rodin/Math/Constants.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

class CommonAdvancedTest : public ::testing::Test
{};

TEST_F(CommonAdvancedTest, PowerAndSqrtUtilities)
{
  EXPECT_DOUBLE_EQ(pow2(3.0), 9.0);
  EXPECT_DOUBLE_EQ(pow<0>(2.0), 1.0);
  EXPECT_DOUBLE_EQ(pow<1>(2.0), 2.0);
  EXPECT_DOUBLE_EQ(pow<5>(2.0), 32.0);
  EXPECT_DOUBLE_EQ(pow(2.0, 3.0), 8.0);
  EXPECT_DOUBLE_EQ(sqrt(16.0), 4.0);
}

TEST_F(CommonAdvancedTest, NaNAndInfUtilities)
{
  const auto rnan = nan<Real>();
  EXPECT_TRUE(isNaN(rnan));
  EXPECT_FALSE(isInf(rnan));

  const auto cnan = nan<Complex>();
  EXPECT_TRUE(isNaN(cnan));

  const Real inf = std::numeric_limits<Real>::infinity();
  EXPECT_TRUE(isInf(inf));
  EXPECT_FALSE(isNaN(inf));
}

TEST_F(CommonAdvancedTest, TrigonometricAndHyperbolicUtilities)
{
  EXPECT_NEAR(cos(0.0), 1.0, 1e-12);
  EXPECT_NEAR(sin(Constants::pi() / 2.0), 1.0, 1e-12);
  EXPECT_NEAR(tan(Constants::pi() / 4.0), 1.0, 1e-12);
  EXPECT_NEAR(cosh(0.0), 1.0, 1e-12);
  EXPECT_NEAR(sinh(0.0), 0.0, 1e-12);
}

TEST_F(CommonAdvancedTest, SignAndCombinatoricsUtilities)
{
  EXPECT_EQ(sgn(-3), -1);
  EXPECT_EQ(sgn(0), 0);
  EXPECT_EQ(sgn(5), 1);

  EXPECT_EQ(binom(5, 2), 10);
  EXPECT_EQ(factorial(0), 1);
  EXPECT_EQ(factorial(6), 720);
  EXPECT_EQ(permutation(5, 3), 60);
}

TEST_F(CommonAdvancedTest, ArithmeticHelperUtilities)
{
  EXPECT_DOUBLE_EQ(sum(1.5, 2.0), 3.5);
  EXPECT_DOUBLE_EQ(minus(3.5, 1.0), 2.5);
  EXPECT_DOUBLE_EQ(minus(2.0), -2.0);
  EXPECT_DOUBLE_EQ(mult(2.0, 4.0), 8.0);
  EXPECT_DOUBLE_EQ(division(10.0, 2.0), 5.0);
}

TEST_F(CommonAdvancedTest, DotOverloadFamilies)
{
  EXPECT_DOUBLE_EQ(dot(2.0, 3.0), 6.0);
  const Complex z1(1.0, 2.0);
  const Complex z2(3.0, -1.0);
  EXPECT_EQ(dot(z1, z2), z1 * std::conj(z2));

  Eigen::Vector3d a;
  a << 1.0, 2.0, 3.0;
  Eigen::Vector3d b;
  b << 4.0, 5.0, 6.0;
  EXPECT_DOUBLE_EQ(dot(a, b), 32.0);

  SpatialVector<Real> sv{1.0, 2.0, 3.0};
  EXPECT_DOUBLE_EQ(dot(sv, b), 32.0);
  EXPECT_DOUBLE_EQ(dot(a, sv), 32.0);

  SpatialMatrix<Real> sm(2, 2);
  sm(0, 0) = 1.0;
  sm(0, 1) = 2.0;
  sm(1, 0) = 3.0;
  sm(1, 1) = 4.0;

  Eigen::Matrix2d em;
  em << 2.0, 0.0,
        1.0, 2.0;
  EXPECT_DOUBLE_EQ(dot(sm, em), 13.0);
  EXPECT_DOUBLE_EQ(dot(em, sm), 13.0);
}

TEST_F(CommonAdvancedTest, CompileTimePowIntegralConstantOverload)
{
  constexpr Real x = 2.0;
  constexpr std::integral_constant<size_t, 4> p{};
  constexpr Real v = pow(x, p);
  static_assert(v == 16.0);
  EXPECT_DOUBLE_EQ(v, 16.0);
}
