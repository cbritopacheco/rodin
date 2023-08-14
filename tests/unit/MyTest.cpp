#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin::Geometry;

TEST(Cat, Miaow) {
  Mesh mesh;
  mesh.initialize(2, 2).finalize();
  EXPECT_EQ(mesh.getDimension(), 2);
}

