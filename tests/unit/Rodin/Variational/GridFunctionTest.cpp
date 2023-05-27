#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

TEST(Rodin_Variational_GridFunction_SanityTest, 1D_Segment_Cos)
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

  GridFunction gf(fes);
}
