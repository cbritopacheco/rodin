#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin::Geometry;

TEST(Rodin_Geometry_Connectivity, SanityTest_2D_3Nodes_Triangles)
{
  constexpr const size_t meshDim = 2;
  constexpr const size_t nodes = 3;

  MeshConnectivity connectivity;
  connectivity.initialize(meshDim)
              .nodes(nodes)
              .polytope(Polytope::Geometry::Triangle, {0, 1, 2});

  EXPECT_EQ(connectivity.getMeshDimension(), 2);

  size_t d;

  d = 2;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(d), 1);
  EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));
  EXPECT_EQ(connectivity.getIncidence({d, 0}, 0), IndexSet({0, 1, 2}));
  EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Geometry::Triangle);

  EXPECT_EQ(connectivity.getIncidence(d, d).size(), 1);
  EXPECT_EQ(connectivity.getIncidence({d, d}, 0).size(), 0);

  connectivity.compute(0, d);
  EXPECT_EQ(connectivity.getIncidence({0, d}, 0), IndexSet({ 0 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 1), IndexSet({ 0 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 2), IndexSet({ 0 }));

  connectivity.compute(d, d - 1);
  EXPECT_EQ(connectivity.getIncidence({d, d - 1}, 0), IndexSet({ {0, 1, 2} }));

  d = 1;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(d), 3);
  EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));

  EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Geometry::Segment);
  EXPECT_EQ(connectivity.getGeometry(d, 1), Polytope::Geometry::Segment);
  EXPECT_EQ(connectivity.getGeometry(d, 2), Polytope::Geometry::Segment);

  d = 0;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(d), 3);
  EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));
  EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Geometry::Point);
  EXPECT_EQ(connectivity.getGeometry(d, 1), Polytope::Geometry::Point);
  EXPECT_EQ(connectivity.getGeometry(d, 2), Polytope::Geometry::Point);
}

TEST(Rodin_Geometry_Connectivity, SanityTest2D_4Nodes_Triangles)
{
  constexpr const size_t meshDim = 2;
  constexpr const size_t nodes = 4;

  MeshConnectivity connectivity;
  connectivity.initialize(meshDim)
              .nodes(nodes)
              .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
              .polytope(Polytope::Geometry::Triangle, {1, 2, 3});

  size_t d;

  d = 2;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(d), 2);
  EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));
  EXPECT_EQ(connectivity.getIncidence({d, 0}, 0), IndexSet({0, 1, 2}));
  EXPECT_EQ(connectivity.getIncidence({d, 0}, 1), IndexSet({1, 2, 3}));

  EXPECT_EQ(connectivity.getIncidence(d, d).size(), 2);
  EXPECT_EQ(connectivity.getIncidence({d, d}, 0).size(), 1);
  EXPECT_EQ(connectivity.getIncidence({d, d}, 1).size(), 1);

  connectivity.compute(0, d);
  EXPECT_EQ(connectivity.getIncidence({0, d}, 0), IndexSet({ 0 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 1), IndexSet({ 0, 1 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 2), IndexSet({ 0, 1 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 3), IndexSet({ 1 }));

  connectivity.compute(d, d - 1);
  EXPECT_EQ(connectivity.getIncidence({d, d - 1}, 0).size(), 3);
  EXPECT_EQ(connectivity.getIncidence({d, d - 1}, 1).size(), 3);

  d = 1;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(d), 5);
  EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));

  d = 0;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(d), 4);
  EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));
}

TEST(Rodin_Geometry_Connectivity, SanityTest2D_5Nodes_Triangles)
{
  constexpr const size_t meshDim = 2;
  constexpr const size_t nodes = 5;

  MeshConnectivity connectivity;
  connectivity.initialize(meshDim)
              .nodes(nodes)
              .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
              .polytope(Polytope::Geometry::Triangle, {1, 2, 3})
              .polytope(Polytope::Geometry::Triangle, {2, 3, 4});

  size_t d;

  d = 2;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(d), 3);
  EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));
  EXPECT_EQ(connectivity.getIncidence({d, 0}, 0), IndexSet({0, 1, 2}));
  EXPECT_EQ(connectivity.getIncidence({d, 0}, 1), IndexSet({1, 2, 3}));
  EXPECT_EQ(connectivity.getIncidence({d, 0}, 2), IndexSet({2, 3, 4}));
  EXPECT_EQ(connectivity.getIncidence(d, d).size(), 3);
  EXPECT_EQ(connectivity.getIncidence({d, d}, 0), IndexSet({1, 2}));
  EXPECT_EQ(connectivity.getIncidence({d, d}, 1), IndexSet({0, 2}));
  EXPECT_EQ(connectivity.getIncidence({d, d}, 2), IndexSet({0, 1}));

  connectivity.compute(0, d);
  EXPECT_EQ(connectivity.getIncidence({0, d}, 0), IndexSet({ 0 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 1), IndexSet({ 0, 1 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 2), IndexSet({ 0, 1, 2 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 3), IndexSet({ 1, 2 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 4), IndexSet({ 2 }));

  connectivity.compute(d, d - 1);
  EXPECT_EQ(connectivity.getIncidence({d, d - 1}, 0).size(), 3);
  EXPECT_EQ(connectivity.getIncidence({d, d - 1}, 1).size(), 3);

  d = 1;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(d), 7);
  EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));

  d = 0;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(d), 5);
  EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));
}

TEST(Rodin_Geometry_Connectivity, SanityTest2D_5Nodes_Mixed)
{
  const size_t meshDim = 2;
  constexpr const size_t nodes = 5;

  MeshConnectivity connectivity;
  connectivity.initialize(meshDim)
              .nodes(nodes)
              .polytope(Polytope::Geometry::Triangle, {0, 1, 4})
              .polytope(Polytope::Geometry::Quadrilateral, {1, 2, 3, 4});

  size_t d;

  d = 2;

  EXPECT_EQ(connectivity.getCount(d), 2);

  connectivity.compute(0, d);
  EXPECT_EQ(connectivity.getIncidence({0, d}, 0), IndexSet({ 0 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 1), IndexSet({ 0, 1 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 2), IndexSet({ 1 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 3), IndexSet({ 1 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 4), IndexSet({ 0, 1 }));
}

TEST(Rodin_Geometry_Connectivity, SanityTest2D_9Nodes_Mixed)
{
  constexpr const size_t meshDim = 2;
  constexpr const size_t nodes = 9;

  MeshConnectivity connectivity;
  connectivity.initialize(meshDim)
              .nodes(nodes)
              .polytope(Polytope::Geometry::Triangle, {0, 1, 4})
              .polytope(Polytope::Geometry::Quadrilateral, {1, 2, 3, 4})
              .polytope(Polytope::Geometry::Quadrilateral, {0, 6, 7, 8})
              .polytope(Polytope::Geometry::Quadrilateral, {0, 4, 5, 6});

  size_t d;

  d = 2;

  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(d), 4);

  connectivity.compute(0, d);
  EXPECT_EQ(connectivity.getIncidence({0, d}, 0), IndexSet({ 0, 2, 3 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 1), IndexSet({ 0, 1 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 2), IndexSet({ 1 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 3), IndexSet({ 1 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 4), IndexSet({ 0, 1, 3 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 5), IndexSet({ 3 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 6), IndexSet({ 2, 3 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 7), IndexSet({ 2 }));
  EXPECT_EQ(connectivity.getIncidence({0, d}, 8), IndexSet({ 2 }));

  d = 1;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(1), 12);

  d = 0;
  connectivity.compute(d, 0);
  EXPECT_EQ(connectivity.getCount(0), 9);
}
