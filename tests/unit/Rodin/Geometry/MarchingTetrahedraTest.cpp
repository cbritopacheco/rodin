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
      .finalize();

    P1 fes(mesh);
    GridFunction ls(fes);
    ls[0] = -1.0;
    ls[1] = 1.0;
    ls[2] = 1.0;
    ls[3] = 1.0;

    MarchingTetrahedra mt;
    mt.setNegativeAttribute(7)
      .setPositiveAttribute(8)
      .setInterfaceAttribute(42);
    const auto split = mt.discretize(ls);

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
}
