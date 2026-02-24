#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/MarchingTriangles.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Geometry_MarchingTriangles, SplitSingleTriangleAndTagInterface)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0.0, 0.0})
      .vertex({1.0, 0.0})
      .vertex({0.0, 1.0})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .attribute({2, 0}, 5)
      .finalize();

    P1 fes(mesh);
    GridFunction ls(fes);
    ls[0] = -1.0;
    ls[1] = 1.0;
    ls[2] = 1.0;

    MarchingTriangles mt(ls);
    mt.setInterface(2, 42)
      .split(2, 5, {7, 8});
    const auto split = mt.discretize();

    EXPECT_GT(split.getCellCount(), 1);

    size_t negativeCount = 0;
    size_t positiveCount = 0;
    for (auto it = split.getCell(); !it.end(); ++it)
    {
      if (it->getAttribute() == 7)
        negativeCount++;
      if (it->getAttribute() == 8)
        positiveCount++;
      EXPECT_EQ(it->getGeometry(), Polytope::Type::Triangle);
    }
    EXPECT_GT(negativeCount, 0);
    EXPECT_GT(positiveCount, 0);

    size_t interfaceCount = 0;
    for (auto it = split.getPolytope(1); !it.end(); ++it)
      if (it->getAttribute() == 42)
        interfaceCount++;
    EXPECT_GT(interfaceCount, 0);
  }

  TEST(Rodin_Geometry_MarchingTriangles, NoSplitCellAttributeKeepsElement)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0.0, 0.0})
      .vertex({1.0, 0.0})
      .vertex({0.0, 1.0})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .attribute({2, 0}, 5)
      .finalize();

    P1 fes(mesh);
    GridFunction ls(fes);
    ls[0] = -1.0;
    ls[1] = 1.0;
    ls[2] = 1.0;

    MarchingTriangles mt(ls);
    mt.noSplit(2, 5);
    const auto split = mt.discretize();

    ASSERT_EQ(split.getCellCount(), 1);
    EXPECT_EQ(split.getCell(0)->getAttribute(), 5);
  }

  TEST(Rodin_Geometry_MarchingTriangles, SplitEdgeAttributesBySign)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0.0, 0.0})
      .vertex({1.0, 0.0})
      .vertex({0.0, 1.0})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .attribute({2, 0}, 5)
      .attribute({1, 0}, 5)
      .attribute({1, 1}, 5)
      .attribute({1, 2}, 5)
      .finalize();

    P1 fes(mesh);
    GridFunction ls(fes);
    ls[0] = -1.0;
    ls[1] = -1.0;
    ls[2] = -1.0;

    MarchingTriangles mt(ls);
    mt.split(2, 5, {7, 8})
      .split(1, 5, {50, 60});
    const auto split = mt.discretize();

    for (auto it = split.getPolytope(1); !it.end(); ++it)
      EXPECT_EQ(it->getAttribute(), 50);
  }
}
