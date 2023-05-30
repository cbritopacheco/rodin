#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

TEST(Rodin_Variational_P1_SanityTest, 2D_Square_Build)
{
  constexpr size_t vdim = 1;
  constexpr size_t sdim = 2;
  constexpr size_t mdim = 2;

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

  EXPECT_EQ(fes.getSize(), 4);
  EXPECT_EQ(fes.getFiniteElement(mdim, 0).getGeometry(), Polytope::Geometry::Triangle);
  EXPECT_EQ(fes.getFiniteElement(mdim, 1).getGeometry(), Polytope::Geometry::Triangle);
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

  EXPECT_NEAR(gf.getValue(0), 0, RODIN_FUZZY_CONSTANT);
  EXPECT_NEAR(gf.getValue(1), 1, RODIN_FUZZY_CONSTANT);
  EXPECT_NEAR(gf.getValue(2), 1, RODIN_FUZZY_CONSTANT);
  EXPECT_NEAR(gf.getValue(3), 2, RODIN_FUZZY_CONSTANT);
}
