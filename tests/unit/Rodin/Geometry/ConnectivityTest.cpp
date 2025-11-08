#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Geometry_Connectivity, SanityTest_2D_3Nodes_Triangles)
  {
    constexpr const size_t meshDim = 2;
    constexpr const size_t nodes = 3;

    Connectivity<Context::Local> connectivity;
    connectivity.initialize(meshDim)
                .nodes(nodes)
                .polytope(Polytope::Type::Triangle, {0, 1, 2});

    EXPECT_EQ(connectivity.getDimension(), 2);

    size_t d;

    d = 2;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 1);
    EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));
    EXPECT_EQ(connectivity.getIncidence({d, 0}, 0), IndexVector({0, 1, 2}));
    EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Type::Triangle);

    connectivity.compute(d, d);
    EXPECT_EQ(connectivity.getIncidence(d, d).size(), 1);
    EXPECT_EQ(connectivity.getIncidence({d, d}, 0).size(), 0);

    connectivity.compute(0, d);
    EXPECT_EQ(connectivity.getIncidence({0, d}, 0), IndexVector({ 0 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 1), IndexVector({ 0 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 2), IndexVector({ 0 }));

    connectivity.compute(d, d - 1);
    EXPECT_EQ(connectivity.getIncidence({d, d - 1}, 0), IndexVector({ {0, 1, 2} }));

    d = 1;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 3);
    EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));

    EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Type::Segment);
    EXPECT_EQ(connectivity.getGeometry(d, 1), Polytope::Type::Segment);
    EXPECT_EQ(connectivity.getGeometry(d, 2), Polytope::Type::Segment);

    d = 0;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 3);
    EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));
    EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Type::Point);
    EXPECT_EQ(connectivity.getGeometry(d, 1), Polytope::Type::Point);
    EXPECT_EQ(connectivity.getGeometry(d, 2), Polytope::Type::Point);
  }

  TEST(Rodin_Geometry_Connectivity, SanityTest2D_4Nodes_Triangles)
  {
    constexpr const size_t meshDim = 2;
    constexpr const size_t nodes = 4;

    Connectivity<Context::Local> connectivity;
    connectivity.initialize(meshDim)
                .nodes(nodes)
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Triangle, {1, 2, 3});

    size_t d;

    d = 2;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 2);
    EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));
    EXPECT_EQ(connectivity.getIncidence({d, 0}, 0), IndexVector({0, 1, 2}));
    EXPECT_EQ(connectivity.getIncidence({d, 0}, 1), IndexVector({1, 2, 3}));

    connectivity.compute(d, d);
    EXPECT_EQ(connectivity.getIncidence(d, d).size(), 2);
    EXPECT_EQ(connectivity.getIncidence({d, d}, 0).size(), 1);
    EXPECT_EQ(connectivity.getIncidence({d, d}, 1).size(), 1);

    connectivity.compute(0, d);
    EXPECT_EQ(connectivity.getIncidence({0, d}, 0), IndexVector({ 0 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 1), IndexVector({ 0, 1 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 2), IndexVector({ 0, 1 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 3), IndexVector({ 1 }));

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

    Connectivity<Context::Local> connectivity;
    connectivity.initialize(meshDim)
                .nodes(nodes)
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Triangle, {1, 2, 3})
                .polytope(Polytope::Type::Triangle, {2, 3, 4});

    size_t d;

    d = 2;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 3);
    EXPECT_EQ(connectivity.getIncidence(d, 0).size(), connectivity.getCount(d));
    EXPECT_EQ(connectivity.getIncidence({d, 0}, 0), IndexVector({0, 1, 2}));
    EXPECT_EQ(connectivity.getIncidence({d, 0}, 1), IndexVector({1, 2, 3}));
    EXPECT_EQ(connectivity.getIncidence({d, 0}, 2), IndexVector({2, 3, 4}));

    connectivity.compute(d, d);
    EXPECT_EQ(connectivity.getIncidence(d, d).size(), 3);
    EXPECT_EQ(connectivity.getIncidence({d, d}, 0), IndexVector({1, 2}));
    EXPECT_EQ(connectivity.getIncidence({d, d}, 1), IndexVector({0, 2}));
    EXPECT_EQ(connectivity.getIncidence({d, d}, 2), IndexVector({0, 1}));

    connectivity.compute(0, d);
    EXPECT_EQ(connectivity.getIncidence({0, d}, 0), IndexVector({ 0 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 1), IndexVector({ 0, 1 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 2), IndexVector({ 0, 1, 2 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 3), IndexVector({ 1, 2 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 4), IndexVector({ 2 }));

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

    Connectivity<Context::Local> connectivity;
    connectivity.initialize(meshDim)
                .nodes(nodes)
                .polytope(Polytope::Type::Triangle, {0, 1, 4})
                .polytope(Polytope::Type::Quadrilateral, {1, 2, 4, 3});

    size_t d;

    d = 2;

    EXPECT_EQ(connectivity.getCount(d), 2);

    connectivity.compute(0, d);
    EXPECT_EQ(connectivity.getIncidence({0, d}, 0), IndexVector({ 0 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 1), IndexVector({ 0, 1 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 2), IndexVector({ 1 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 3), IndexVector({ 1 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 4), IndexVector({ 0, 1 }));
  }

  TEST(Rodin_Geometry_Connectivity, SanityTest2D_9Nodes_Mixed)
  {
    constexpr const size_t meshDim = 2;
    constexpr const size_t nodes = 9;

    Connectivity<Context::Local> connectivity;
    connectivity.initialize(meshDim)
                .nodes(nodes)
                .polytope(Polytope::Type::Triangle, {0, 1, 4})
                .polytope(Polytope::Type::Quadrilateral, {1, 2, 4, 3})
                .polytope(Polytope::Type::Quadrilateral, {8, 0, 7, 6})
                .polytope(Polytope::Type::Quadrilateral, {0, 4, 6, 5});

    size_t d;

    d = 2;

    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 4);

    connectivity.compute(0, d);
    EXPECT_EQ(connectivity.getIncidence({0, d}, 0), IndexVector({ 0, 2, 3 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 1), IndexVector({ 0, 1 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 2), IndexVector({ 1 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 3), IndexVector({ 1 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 4), IndexVector({ 0, 1, 3 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 5), IndexVector({ 3 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 6), IndexVector({ 2, 3 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 7), IndexVector({ 2 }));
    EXPECT_EQ(connectivity.getIncidence({0, d}, 8), IndexVector({ 2 }));

    d = 1;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(1), 14);

    d = 0;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(0), 9);
  }

  // Additional comprehensive tests

  TEST(Rodin_Geometry_Connectivity, Quadrilateral_2D)
  {
    constexpr const size_t meshDim = 2;
    constexpr const size_t nodes = 4;

    Connectivity<Context::Local> connectivity;
    connectivity.initialize(meshDim)
                .nodes(nodes)
                .polytope(Polytope::Type::Quadrilateral, {0, 1, 3, 2});

    EXPECT_EQ(connectivity.getDimension(), 2);

    size_t d = 2;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 1);
    EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Type::Quadrilateral);

    d = 1;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 4);

    d = 0;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 4);
  }

  TEST(Rodin_Geometry_Connectivity, Tetrahedron_3D)
  {
    constexpr const size_t meshDim = 3;
    constexpr const size_t nodes = 4;

    Connectivity<Context::Local> connectivity;
    connectivity.initialize(meshDim)
                .nodes(nodes)
                .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3});

    EXPECT_EQ(connectivity.getDimension(), 3);

    size_t d = 3;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 1);
    EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Type::Tetrahedron);

    d = 2;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 4);

    d = 1;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 6);

    d = 0;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 4);
  }

  TEST(Rodin_Geometry_Connectivity, EmptyMesh)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2).nodes(0);

    size_t d = 2;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 0);
  }

  TEST(Rodin_Geometry_Connectivity, MultipleComputeCalls)
  {
    constexpr const size_t meshDim = 2;
    constexpr const size_t nodes = 3;

    Connectivity<Context::Local> connectivity;
    connectivity.initialize(meshDim)
                .nodes(nodes)
                .polytope(Polytope::Type::Triangle, {0, 1, 2});

    size_t d = 2;
    connectivity.compute(d, 0);
    size_t count1 = connectivity.getCount(d);

    // Compute again, should be idempotent
    connectivity.compute(d, 0);
    size_t count2 = connectivity.getCount(d);

    EXPECT_EQ(count1, count2);
  }

  // EDGE CASE TESTS - Starting with the most basic edge cases

  TEST(Rodin_Geometry_Connectivity, EdgeCase_ZeroNodes)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2).nodes(0);

    EXPECT_EQ(connectivity.getDimension(), 0);
    EXPECT_EQ(connectivity.getCount(0), 0);
    EXPECT_EQ(connectivity.getCount(1), 0);
    EXPECT_EQ(connectivity.getCount(2), 0);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_SingleNode_NoPolytopes)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2).nodes(1);

    connectivity.compute(0, 0);
    EXPECT_EQ(connectivity.getCount(0), 1);
    EXPECT_EQ(connectivity.getCount(2), 0);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_MultipleNodes_NoPolytopes)
  {
    constexpr const size_t nodes = 10;
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2).nodes(nodes);

    connectivity.compute(0, 0);
    EXPECT_EQ(connectivity.getCount(0), nodes);
    EXPECT_EQ(connectivity.getCount(1), 0);
    EXPECT_EQ(connectivity.getCount(2), 0);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_1D_SingleSegment)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(1)
                .nodes(2)
                .polytope(Polytope::Type::Segment, {0, 1});

    EXPECT_EQ(connectivity.getDimension(), 1);

    size_t d = 1;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 1);
    EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Type::Segment);

    d = 0;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 2);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_1D_MultipleSegments)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(1)
                .nodes(4)
                .polytope(Polytope::Type::Segment, {0, 1})
                .polytope(Polytope::Type::Segment, {1, 2})
                .polytope(Polytope::Type::Segment, {2, 3});

    size_t d = 1;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 3);

    connectivity.compute(d, d);
    EXPECT_EQ(connectivity.getIncidence(d, d).size(), 3);
    // Middle segments should have neighbors
    EXPECT_EQ(connectivity.getIncidence({d, d}, 1).size(), 2);

    d = 0;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 4);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_3D_SingleTetrahedron)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(3)
                .nodes(4)
                .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3});

    EXPECT_EQ(connectivity.getDimension(), 3);

    size_t d = 3;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 1);
    EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Type::Tetrahedron);

    // 4 triangular faces
    d = 2;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 4);
    EXPECT_EQ(connectivity.getGeometry(d, 0), Polytope::Type::Triangle);

    // 6 edges
    d = 1;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 6);

    // 4 vertices
    d = 0;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 4);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_CopyConstruction)
  {
    Connectivity<Context::Local> connectivity1;
    connectivity1.initialize(2)
                 .nodes(3)
                 .polytope(Polytope::Type::Triangle, {0, 1, 2});

    connectivity1.compute(2, 0);
    connectivity1.compute(1, 0);
    connectivity1.compute(0, 0);

    // Copy construct
    Connectivity<Context::Local> connectivity2(connectivity1);

    EXPECT_EQ(connectivity2.getDimension(), connectivity1.getDimension());
    EXPECT_EQ(connectivity2.getCount(0), connectivity1.getCount(0));
    EXPECT_EQ(connectivity2.getCount(1), connectivity1.getCount(1));
    EXPECT_EQ(connectivity2.getCount(2), connectivity1.getCount(2));
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_MoveConstruction)
  {
    Connectivity<Context::Local> connectivity1;
    connectivity1.initialize(2)
                 .nodes(3)
                 .polytope(Polytope::Type::Triangle, {0, 1, 2});

    connectivity1.compute(2, 0);
    size_t originalCount = connectivity1.getCount(2);

    // Move construct
    Connectivity<Context::Local> connectivity2(std::move(connectivity1));

    EXPECT_EQ(connectivity2.getDimension(), 2);
    EXPECT_EQ(connectivity2.getCount(2), originalCount);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_CopyAssignment)
  {
    Connectivity<Context::Local> connectivity1;
    connectivity1.initialize(2)
                 .nodes(3)
                 .polytope(Polytope::Type::Triangle, {0, 1, 2});

    connectivity1.compute(2, 0);
    connectivity1.compute(1, 0);

    Connectivity<Context::Local> connectivity2;
    connectivity2 = connectivity1;

    EXPECT_EQ(connectivity2.getDimension(), connectivity1.getDimension());
    EXPECT_EQ(connectivity2.getCount(1), connectivity1.getCount(1));
    EXPECT_EQ(connectivity2.getCount(2), connectivity1.getCount(2));
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_MoveAssignment)
  {
    Connectivity<Context::Local> connectivity1;
    connectivity1.initialize(2)
                 .nodes(4)
                 .polytope(Polytope::Type::Quadrilateral, {0, 1, 3, 2});

    connectivity1.compute(2, 0);
    size_t originalCount = connectivity1.getCount(2);

    Connectivity<Context::Local> connectivity2;
    connectivity2 = std::move(connectivity1);

    EXPECT_EQ(connectivity2.getDimension(), 2);
    EXPECT_EQ(connectivity2.getCount(2), originalCount);
    EXPECT_EQ(connectivity2.getGeometry(2, 0), Polytope::Type::Quadrilateral);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_DegenerateTriangle_DuplicateVertex)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(2)
                .polytope(Polytope::Type::Triangle, {0, 0, 1});

    connectivity.compute(2, 0);
    EXPECT_EQ(connectivity.getCount(2), 1);

    connectivity.compute(1, 0);
    // Degenerate triangle should still produce edges (even if degenerate)
    EXPECT_GT(connectivity.getCount(1), 0);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_ComputeOrder_Dependency)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(3)
                .polytope(Polytope::Type::Triangle, {0, 1, 2});

    // Compute higher dimension first
    connectivity.compute(2, 0);
    EXPECT_EQ(connectivity.getCount(2), 1);

    // Then compute lower dimensions
    connectivity.compute(1, 0);
    EXPECT_EQ(connectivity.getCount(1), 3);

    connectivity.compute(0, 0);
    EXPECT_EQ(connectivity.getCount(0), 3);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_GetCountByGeometry_Triangle)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(5)
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Triangle, {1, 2, 3})
                .polytope(Polytope::Type::Quadrilateral, {2, 3, 4, 1});

    connectivity.compute(2, 0);

    EXPECT_EQ(connectivity.getCount(Polytope::Type::Triangle), 2);
    EXPECT_EQ(connectivity.getCount(Polytope::Type::Quadrilateral), 1);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_AllQuadrilaterals)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(6)
                .polytope(Polytope::Type::Quadrilateral, {0, 1, 4, 3})
                .polytope(Polytope::Type::Quadrilateral, {1, 2, 5, 4});

    size_t d = 2;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 2);
    EXPECT_EQ(connectivity.getCount(Polytope::Type::Quadrilateral), 2);
    EXPECT_EQ(connectivity.getCount(Polytope::Type::Triangle), 0);

    // Check shared edge (interface)
    connectivity.compute(d, d);
    EXPECT_EQ(connectivity.getIncidence({d, d}, 0).size(), 1);
    EXPECT_EQ(connectivity.getIncidence({d, d}, 1).size(), 1);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_3D_TwoTets_SharedFace)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(3)
                .nodes(5)
                .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
                .polytope(Polytope::Type::Tetrahedron, {1, 2, 3, 4});

    size_t d = 3;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 2);

    // Check for shared face
    connectivity.compute(d, d);
    EXPECT_EQ(connectivity.getIncidence({d, d}, 0).size(), 1);
    EXPECT_EQ(connectivity.getIncidence({d, d}, 1).size(), 1);

    // 7 total faces (4+4-1 shared)
    d = 2;
    connectivity.compute(d, 0);
    EXPECT_EQ(connectivity.getCount(d), 7);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_ComputeTranspose)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(4)
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Triangle, {1, 2, 3});

    // Compute cell to vertex
    connectivity.compute(2, 0);
    EXPECT_EQ(connectivity.getIncidence({2, 0}, 0).size(), 3);

    // Compute transpose: vertex to cell
    connectivity.compute(0, 2);
    EXPECT_EQ(connectivity.getIncidence({0, 2}, 0).size(), 1);
    EXPECT_EQ(connectivity.getIncidence({0, 2}, 1).size(), 2);
    EXPECT_EQ(connectivity.getIncidence({0, 2}, 2).size(), 2);
    EXPECT_EQ(connectivity.getIncidence({0, 2}, 3).size(), 1);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_Clear_Specific_Connectivity)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(3)
                .polytope(Polytope::Type::Triangle, {0, 1, 2});

    connectivity.compute(2, 0);
    connectivity.compute(1, 0);

    EXPECT_EQ(connectivity.getCount(2), 1);
    EXPECT_EQ(connectivity.getCount(1), 3);

    // Clear edge connectivity
    connectivity.clear(1, 0);

    // Cell connectivity should remain
    EXPECT_EQ(connectivity.getCount(2), 1);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_Reserve_Performance)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2).nodes(100);

    // Reserve space for many cells
    connectivity.reserve(2, 50);

    // Add multiple triangles
    for (size_t i = 0; i < 10; ++i)
    {
      connectivity.polytope(Polytope::Type::Triangle, {i, i+1, i+2});
    }

    connectivity.compute(2, 0);
    EXPECT_EQ(connectivity.getCount(2), 10);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_GetGeometry_AllDimensions)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(3)
                .polytope(Polytope::Type::Triangle, {0, 1, 2});

    connectivity.compute(2, 0);
    connectivity.compute(1, 0);
    connectivity.compute(0, 0);

    // Check geometry types at all dimensions
    EXPECT_EQ(connectivity.getGeometry(2, 0), Polytope::Type::Triangle);
    EXPECT_EQ(connectivity.getGeometry(1, 0), Polytope::Type::Segment);
    EXPECT_EQ(connectivity.getGeometry(0, 0), Polytope::Type::Point);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_GetPolytope_VertexIndices)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(4)
                .polytope(Polytope::Type::Quadrilateral, {0, 1, 3, 2});

    connectivity.compute(2, 0);

    const auto& polytope = connectivity.getPolytope(2, 0);
    EXPECT_EQ(polytope.size(), 4);
    EXPECT_EQ(polytope[0], 0);
    EXPECT_EQ(polytope[1], 1);
    EXPECT_EQ(polytope[2], 3);
    EXPECT_EQ(polytope[3], 2);
  }

  TEST(Rodin_Geometry_Connectivity, EdgeCase_InterfaceBoundary_SingleCell)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(3)
                .polytope(Polytope::Type::Triangle, {0, 1, 2});

    size_t d = 2;
    connectivity.compute(d, d);

    // Single cell has no neighbors (interfaces)
    EXPECT_EQ(connectivity.getIncidence({d, d}, 0).size(), 0);
  }
}
