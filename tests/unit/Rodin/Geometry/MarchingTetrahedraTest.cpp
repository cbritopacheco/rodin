#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/MarchingTetrahedra.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Geometry_MarchingTetrahedra, SplitSingleTetrahedronAndTagInterface)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(3)
      .nodes(4)
      .vertex({0.0, 0.0, 0.0})
      .vertex({1.0, 0.0, 0.0})
      .vertex({0.0, 1.0, 0.0})
      .vertex({0.0, 0.0, 1.0})
      .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
      .attribute({3, 0}, 5)
      .finalize();

    P1 fes(mesh);
    GridFunction ls(fes);
    ls[0] = -1.0;
    ls[1] = 1.0;
    ls[2] = 1.0;
    ls[3] = 1.0;

    MarchingTetrahedra mt(ls);
    mt.split(3, 5, {7, 8})
      .setInterface(42);
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
      EXPECT_EQ(it->getGeometry(), Polytope::Type::Tetrahedron);
    }
    EXPECT_GT(negativeCount, 0);
    EXPECT_GT(positiveCount, 0);

    size_t interfaceCount = 0;
    for (auto it = split.getFace(); !it.end(); ++it)
    {
      if (it->getAttribute() == 42)
        interfaceCount++;
    }
    EXPECT_GT(interfaceCount, 0);
  }

  TEST(Rodin_Geometry_MarchingTetrahedra, NoSplitCellAttributeKeepsElement)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(3)
      .nodes(4)
      .vertex({0.0, 0.0, 0.0})
      .vertex({1.0, 0.0, 0.0})
      .vertex({0.0, 1.0, 0.0})
      .vertex({0.0, 0.0, 1.0})
      .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
      .attribute({3, 0}, 5)
      .finalize();

    P1 fes(mesh);
    GridFunction ls(fes);
    ls[0] = -1.0;
    ls[1] = 1.0;
    ls[2] = 1.0;
    ls[3] = 1.0;

    MarchingTetrahedra mt(ls);
    mt.noSplit(3, 5);
    const auto split = mt.discretize();

    ASSERT_EQ(split.getCellCount(), 1);
    EXPECT_EQ(split.getCell(0)->getAttribute(), 5);
  }

  TEST(Rodin_Geometry_MarchingTetrahedra, SplitFaceAndEdgeAttributesBySign)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(3)
      .nodes(4)
      .vertex({0.0, 0.0, 0.0})
      .vertex({1.0, 0.0, 0.0})
      .vertex({0.0, 1.0, 0.0})
      .vertex({0.0, 0.0, 1.0})
      .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
      .attribute({3, 0}, 5)
      .finalize();

    P1 fes(mesh);
    GridFunction ls(fes);
    ls[0] = -1.0;
    ls[1] = -1.0;
    ls[2] = -1.0;
    ls[3] = -1.0;

    MarchingTetrahedra mt(ls);
    mt.split(3, 5, {7, 8})
      .split(2, 5, {50, 60})
      .split(1, 5, {70, 80});
    const auto split = mt.discretize();

    for (auto it = split.getFace(); !it.end(); ++it)
      EXPECT_EQ(it->getAttribute(), 50);
    for (auto it = split.getPolytope(1); !it.end(); ++it)
      EXPECT_EQ(it->getAttribute(), 70);
  }
}
