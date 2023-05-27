#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

TEST(Rodin_Variational_P1Element_SanityTest, 1D_Reference_Segment)
{
  P1Element k(Polytope::Geometry::Segment);

  {
    const auto phi = k.getBasis(Math::Vector{{0}});
    EXPECT_EQ(phi.size(), 2);
    EXPECT_NEAR(phi(0), 1, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(1), 0, RODIN_FUZZY_CONSTANT);
  }

  {
    const auto phi = k.getBasis(Math::Vector{{0.5}});
    EXPECT_EQ(phi.size(), 2);
    EXPECT_NEAR(phi(0), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(1), 0.5, RODIN_FUZZY_CONSTANT);
  }

  {
    const auto phi = k.getBasis(Math::Vector{{1}});
    EXPECT_EQ(phi.size(), 2);
    EXPECT_NEAR(phi(0), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(1), 1, RODIN_FUZZY_CONSTANT);
  }
}

TEST(Rodin_Variational_P1Element_FuzzyTest, 1D_Reference_Segment)
{
  constexpr size_t n = 25;

  RandomFloat gen(0.0, 1.0);
  P1Element k(Polytope::Geometry::Segment);

  {
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      Math::Vector p = (1 - s) * Math::Vector{{1}};
      const auto phi = k.getBasis(p);
      EXPECT_EQ(phi.size(), 2);
      EXPECT_NEAR(phi(0), 1 - p.x(), RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(phi(1), p.x(), RODIN_FUZZY_CONSTANT);
    }
  }
}

TEST(Rodin_Variational_P1Element_SanityTest, 2D_Reference_Triangle)
{
  P1Element k(Polytope::Geometry::Triangle);

  {
    const auto phi = k.getBasis(Math::Vector{{0, 0}});
    EXPECT_EQ(phi.size(), 3);
    EXPECT_NEAR(phi(0), 1, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(1), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(2), 0, RODIN_FUZZY_CONSTANT);
  }

  {
    const auto phi = k.getBasis(Math::Vector{{1, 0}});
    EXPECT_EQ(phi.size(), 3);
    EXPECT_NEAR(phi(0), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(1), 1, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(2), 0, RODIN_FUZZY_CONSTANT);
  }

  {
    const auto phi = k.getBasis(Math::Vector{{0, 1}});
    EXPECT_EQ(phi.size(), 3);
    EXPECT_NEAR(phi(0), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(1), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(2), 1, RODIN_FUZZY_CONSTANT);
  }

  {
    const auto phi = k.getBasis(Math::Vector{{0.5, 0}});
    EXPECT_EQ(phi.size(), 3);
    EXPECT_NEAR(phi(0), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(1), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(2), 0, RODIN_FUZZY_CONSTANT);
  }

  {
    const auto phi = k.getBasis(Math::Vector{{0.5, 0.5}});
    EXPECT_EQ(phi.size(), 3);
    EXPECT_NEAR(phi(0), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(1), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(2), 0.5, RODIN_FUZZY_CONSTANT);
  }

  {
    const auto phi = k.getBasis(Math::Vector{{0, 0.5}});
    EXPECT_EQ(phi.size(), 3);
    EXPECT_NEAR(phi(0), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(1), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(phi(2), 0.5, RODIN_FUZZY_CONSTANT);
  }
}

TEST(Rodin_Variational_P1Element_FuzzyTest, 2D_Reference_Triangle)
{
  constexpr size_t n = 25;

  RandomFloat gen(0.0, 1.0);
  P1Element k(Polytope::Geometry::Triangle);

  {
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      Math::Vector p = s * Math::Vector{{0, 0}} + (1 - s) * Math::Vector{{1, 0}};
      const auto phi = k.getBasis(p);
      EXPECT_EQ(phi.size(), 3);
      EXPECT_NEAR(phi(0), -p.x() - p.y() + 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(phi(1), p.x(), RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(phi(2), p.y(), RODIN_FUZZY_CONSTANT);
    }
  }

  {
    RandomFloat gen(0.0, 1.0);
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      Math::Vector p = s * Math::Vector{{1, 0}} + (1 - s) * Math::Vector{{0, 1}};
      const auto phi = k.getBasis(p);
      EXPECT_EQ(phi.size(), 3);
      EXPECT_NEAR(phi(0), -p.x() - p.y() + 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(phi(1), p.x(), RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(phi(2), p.y(), RODIN_FUZZY_CONSTANT);
    }
  }

  {
    RandomFloat gen(0.0, 1.0);
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      Math::Vector p = s * Math::Vector{{0, 1}} + (1 - s) * Math::Vector{{0, 0}};
      const auto phi = k.getBasis(p);
      EXPECT_EQ(phi.size(), 3);
      EXPECT_NEAR(phi(0), -p.x() - p.y() + 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(phi(1), p.x(), RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(phi(2), p.y(), RODIN_FUZZY_CONSTANT);
    }
  }
}

TEST(Rodin_Variational_P1_SanityTest, 2D_Square_Build)
{
  constexpr size_t vdim = 1;
  constexpr size_t mdim = 2;
  constexpr size_t sdim = 2;

  Mesh mesh =
    Mesh<Rodin::Context::Serial>::Builder()
    .initialize(sdim)
    .nodes(4)
    .vertex({0, 0})
    .vertex({1, 0})
    .vertex({0, 1})
    .vertex({1, 1})
    .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
    .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
    .finalize();

  EXPECT_EQ(mesh.getDimension(), mdim);
  EXPECT_EQ(mesh.getSpaceDimension(), sdim);

  P1 fes(mesh);

  EXPECT_EQ(fes.getVectorDimension(), vdim);
  EXPECT_EQ(fes.getSize(), mesh.getVertexCount());
}

TEST(Rodin_Variational_P1_GridFunction_SanityTest, 2D_Square_Project_Constant)
{
  constexpr size_t sdim = 2;

  Mesh mesh =
    Mesh<Rodin::Context::Serial>::Builder()
    .initialize(sdim)
    .nodes(4)
    .vertex({0, 0})
    .vertex({1, 0})
    .vertex({0, 1})
    .vertex({1, 1})
    .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
    .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
    .finalize();

  P1 fes(mesh);
  GridFunction gf(fes);

  ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
  gf.project(c);

  EXPECT_EQ(gf.getRangeShape(), RangeShape(1, 1));
}
