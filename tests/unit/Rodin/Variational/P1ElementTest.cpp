#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational/P1.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

TEST(Rodin_Variational_ScalarP1Element_SanityTest, 1D_Reference_Segment)
{
  ScalarP1Element k(Polytope::Geometry::Segment);

  {
    EXPECT_NEAR(k.getBasis(0)(Math::Vector{{0}}), 1, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(1)(Math::Vector{{0}}), 0, RODIN_FUZZY_CONSTANT);
  }

  {
    EXPECT_NEAR(k.getBasis(0)(Math::Vector{{0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(1)(Math::Vector{{0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
  }

  {
    EXPECT_NEAR(k.getBasis(0)(Math::Vector{{1}}), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(1)(Math::Vector{{1}}), 1, RODIN_FUZZY_CONSTANT);
  }
}

TEST(Rodin_Variational_ScalarP1Element_FuzzyTest, 1D_Reference_Segment)
{
  constexpr size_t n = 25;

  RandomFloat gen(0.0, 1.0);
  ScalarP1Element k(Polytope::Geometry::Segment);

  {
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      Math::Vector p = (1 - s) * Math::Vector{{1}};
      EXPECT_NEAR(k.getBasis(0)(p), 1 - p.x(), RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(p), p.x(), RODIN_FUZZY_CONSTANT);
    }
  }
}

TEST(Rodin_Variational_ScalarP1Element_SanityTest, 2D_Reference_Triangle)
{
  ScalarP1Element k(Polytope::Geometry::Triangle);

  {
    EXPECT_NEAR(k.getBasis(0)(Math::Vector{{0, 0}}), 1, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(1)(Math::Vector{{0, 0}}), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(2)(Math::Vector{{0, 0}}), 0, RODIN_FUZZY_CONSTANT);
  }

  {
    EXPECT_NEAR(k.getBasis(0)(Math::Vector{{1, 0}}), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(1)(Math::Vector{{1, 0}}), 1, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(2)(Math::Vector{{1, 0}}), 0, RODIN_FUZZY_CONSTANT);
  }

  {
    EXPECT_NEAR(k.getBasis(0)(Math::Vector{{0, 1}}), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(1)(Math::Vector{{0, 1}}), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(2)(Math::Vector{{0, 1}}), 1, RODIN_FUZZY_CONSTANT);
  }

  {
    EXPECT_NEAR(k.getBasis(0)(Math::Vector{{0.5, 0}}), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(1)(Math::Vector{{0.5, 0}}), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(2)(Math::Vector{{0.5, 0}}), 0, RODIN_FUZZY_CONSTANT);
  }

  {
    EXPECT_NEAR(k.getBasis(0)(Math::Vector{{0.5, 0.5}}), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(1)(Math::Vector{{0.5, 0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(2)(Math::Vector{{0.5, 0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
  }

  {
    EXPECT_NEAR(k.getBasis(0)(Math::Vector{{0, 0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(1)(Math::Vector{{0, 0.5}}), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(k.getBasis(2)(Math::Vector{{0, 0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
  }
}

TEST(Rodin_Variational_ScalarP1Element_FuzzyTest, 2D_Reference_Triangle)
{
  constexpr size_t n = 25;

  RandomFloat gen(0.0, 1.0);
  ScalarP1Element k(Polytope::Geometry::Triangle);

  {
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      Math::Vector p = s * Math::Vector{{0, 0}} + (1 - s) * Math::Vector{{1, 0}};
      EXPECT_NEAR(k.getBasis(0)(p), -p.x() - p.y() + 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(p), p.x(), RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(p), p.y(), RODIN_FUZZY_CONSTANT);
    }
  }
}

