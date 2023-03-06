#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin::Geometry;

TEST(Cat, Miaow) {
  Mesh mesh;
  mesh.load("/Users/carlos/Projects/rodin/resources/mfem/StarSquare.mfem.mesh");
  std::cout << mesh.getSpaceDimension() << std::endl;
  EXPECT_EQ(mesh.getDimension(), 2);
}

