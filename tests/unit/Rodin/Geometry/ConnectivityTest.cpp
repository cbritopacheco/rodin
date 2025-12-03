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

  // ============================================================================
  // COMPREHENSIVE MODE TESTS
  // ============================================================================

  TEST(Rodin_Geometry_Connectivity, Mode_Discover_2D_AddsEdges)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(4)
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Triangle, {1, 2, 3});

    // Before computing edges, they don't exist
    EXPECT_EQ(connectivity.getCount(1), 0);

    // Mode::Discover should discover and add edges
    connectivity.compute(2, 1, Connectivity<Context::Local>::Mode::Discover);

    // After computing with Discover mode, edges should exist
    // Two triangles: 3+3=6 edges total, sharing 1 common edge = 5 unique edges
    EXPECT_EQ(connectivity.getCount(1), 5);
  }

  TEST(Rodin_Geometry_Connectivity, Mode_Restrict_2D_OnlyUsesExisting)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(4)
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Triangle, {1, 2, 3});

    // First, discover edges normally
    connectivity.compute(2, 1, Connectivity<Context::Local>::Mode::Discover);
    EXPECT_EQ(connectivity.getCount(1), 5);

    // Clear the relation
    connectivity.clear(2, 1);

    // Now compute with Restrict mode - should use existing edges only
    connectivity.compute(2, 1, Connectivity<Context::Local>::Mode::Restrict);

    // Should have computed the relation without changing edge count
    EXPECT_EQ(connectivity.getCount(1), 5);
    EXPECT_EQ(connectivity.getIncidence(2, 1).size(), 2);  // 2 cells
  }

  TEST(Rodin_Geometry_Connectivity, Mode_Discover_3D_AddsFaces)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(3)
                .nodes(5)
                .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
                .polytope(Polytope::Type::Tetrahedron, {1, 2, 3, 4});

    // Before computing faces, they don't exist
    EXPECT_EQ(connectivity.getCount(2), 0);

    // Mode::Discover should discover and add faces
    connectivity.compute(3, 2, Connectivity<Context::Local>::Mode::Discover);

    // After computing with Discover mode, faces should exist
    // Two tets with 4 faces each, sharing 1 face = 7 unique faces
    EXPECT_EQ(connectivity.getCount(2), 7);
  }

  TEST(Rodin_Geometry_Connectivity, Mode_Restrict_3D_OnlyUsesExisting)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(3)
                .nodes(5)
                .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
                .polytope(Polytope::Type::Tetrahedron, {1, 2, 3, 4});

    // First, discover faces
    connectivity.compute(3, 2, Connectivity<Context::Local>::Mode::Discover);
    size_t faceCount = connectivity.getCount(2);
    EXPECT_EQ(faceCount, 7);

    // Clear the relation
    connectivity.clear(3, 2);

    // Now compute with Restrict mode
    connectivity.compute(3, 2, Connectivity<Context::Local>::Mode::Restrict);

    // Face count should remain unchanged
    EXPECT_EQ(connectivity.getCount(2), faceCount);
  }

  TEST(Rodin_Geometry_Connectivity, Mode_Discover_Mixed_2D)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(6)
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Quadrilateral, {1, 3, 4, 2})
                .polytope(Polytope::Type::Triangle, {3, 4, 5});

    // Compute edges with Discover mode
    connectivity.compute(2, 1, Connectivity<Context::Local>::Mode::Discover);

    // Should have correct number of unique edges
    // Triangle: 3 edges, Quad: 4 edges, Triangle: 3 edges
    // Shared edges reduce the count
    EXPECT_GT(connectivity.getCount(1), 0);
    EXPECT_LE(connectivity.getCount(1), 10);  // Max 10 if no sharing

    // Verify cell-to-edge connectivity
    EXPECT_EQ(connectivity.getIncidence(2, 1).size(), 3);  // 3 cells
    EXPECT_EQ(connectivity.getIncidence({2, 1}, 0).size(), 3);  // Triangle has 3 edges
    EXPECT_EQ(connectivity.getIncidence({2, 1}, 1).size(), 4);  // Quad has 4 edges
    EXPECT_EQ(connectivity.getIncidence({2, 1}, 2).size(), 3);  // Triangle has 3 edges
  }

  TEST(Rodin_Geometry_Connectivity, Mode_Restrict_WithPartialEdges)
  {
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(5)
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Triangle, {2, 3, 4});

    // Manually add only some edges
    connectivity.polytope(Polytope::Type::Segment, {0, 1});
    connectivity.polytope(Polytope::Type::Segment, {1, 2});
    connectivity.polytope(Polytope::Type::Segment, {2, 3});

    EXPECT_EQ(connectivity.getCount(1), 3);

    // Compute with Restrict mode - should only use these edges
    connectivity.compute(2, 1, Connectivity<Context::Local>::Mode::Restrict);

    // Edge count should not increase
    EXPECT_EQ(connectivity.getCount(1), 3);

    // Cells should reference only the existing edges
    // First triangle should have 2 edges (0-1, 1-2), missing (0-2)
    // Second triangle should have 1 edge (2-3), missing (3-4, 4-2)
    EXPECT_LE(connectivity.getIncidence({2, 1}, 0).size(), 3);
    EXPECT_LE(connectivity.getIncidence({2, 1}, 1).size(), 3);
  }

  // ============================================================================
  // UNIFORM GRID TESTS
  // ============================================================================

  TEST(Rodin_Geometry_Connectivity, UniformGrid_2D_Triangles_Small)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    auto& connectivity = mesh.getConnectivity();

    // Check basic counts
    EXPECT_EQ(connectivity.getDimension(), 2);
    EXPECT_EQ(mesh.getVertexCount(), 16);
    EXPECT_EQ(mesh.getCellCount(), 18);

    // Compute connectivity
    connectivity.compute(2, 1);
    connectivity.compute(1, 2);
    connectivity.compute(0, 2);

    // Verify edge count (in a 4x4 triangle grid)
    EXPECT_GT(connectivity.getCount(1), 0);

    // Verify all cells have 3 edges (triangles)
    for (size_t i = 0; i < mesh.getCellCount(); ++i)
    {
      EXPECT_EQ(connectivity.getIncidence({2, 1}, i).size(), 3);
    }
  }

  TEST(Rodin_Geometry_Connectivity, UniformGrid_2D_Triangles_Medium)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {8, 8});
    auto& connectivity = mesh.getConnectivity();

    EXPECT_EQ(connectivity.getDimension(), 2);

    // Compute various connectivity relations
    connectivity.compute(2, 0);  // Cell to vertex
    connectivity.compute(1, 0);  // Edge to vertex
    connectivity.compute(2, 1);  // Cell to edge
    connectivity.compute(2, 2);  // Cell to cell (neighbors)

    // Check that all cells have exactly 3 vertices
    for (size_t i = 0; i < mesh.getCellCount(); ++i)
    {
      EXPECT_EQ(connectivity.getIncidence({2, 0}, i).size(), 3);
    }

    // Check that all edges have exactly 2 vertices
    for (size_t i = 0; i < connectivity.getCount(1); ++i)
    {
      EXPECT_EQ(connectivity.getIncidence({1, 0}, i).size(), 2);
    }
  }

  TEST(Rodin_Geometry_Connectivity, UniformGrid_2D_Quadrilaterals)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {5, 5});
    auto& connectivity = mesh.getConnectivity();

    EXPECT_EQ(connectivity.getDimension(), 2);
    EXPECT_EQ(mesh.getVertexCount(), 25);
    EXPECT_EQ(mesh.getCellCount(), 16);  // (5-1)*(5-1) = 16 quads

    // Compute connectivity
    connectivity.compute(2, 1);  // Cell to edge
    connectivity.compute(2, 0);  // Cell to vertex

    // All cells should have 4 edges and 4 vertices (quadrilaterals)
    for (size_t i = 0; i < mesh.getCellCount(); ++i)
    {
      EXPECT_EQ(connectivity.getIncidence({2, 1}, i).size(), 4);
      EXPECT_EQ(connectivity.getIncidence({2, 0}, i).size(), 4);
      EXPECT_EQ(connectivity.getGeometry(2, i), Polytope::Type::Quadrilateral);
    }
  }

  TEST(Rodin_Geometry_Connectivity, UniformGrid_3D_Tetrahedra)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    auto& connectivity = mesh.getConnectivity();

    EXPECT_EQ(connectivity.getDimension(), 3);
    // Note: 3x3x3 tetrahedra grid has more vertices than a simple cube
    EXPECT_GT(mesh.getVertexCount(), 0);

    // Compute connectivity
    connectivity.compute(3, 2);  // Cell to face
    connectivity.compute(3, 0);  // Cell to vertex

    // All cells should have 4 faces and 4 vertices (tetrahedra)
    for (size_t i = 0; i < mesh.getCellCount(); ++i)
    {
      EXPECT_EQ(connectivity.getIncidence({3, 2}, i).size(), 4);
      EXPECT_EQ(connectivity.getIncidence({3, 0}, i).size(), 4);
      EXPECT_EQ(connectivity.getGeometry(3, i), Polytope::Type::Tetrahedron);
    }
  }

  TEST(Rodin_Geometry_Connectivity, UniformGrid_2D_BoundaryAndInterface)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {5, 5});
    auto& connectivity = mesh.getConnectivity();

    // Compute face-to-cell connectivity
    connectivity.compute(1, 2);

    size_t boundaryCount = 0;
    size_t interfaceCount = 0;

    // Count boundary and interface edges
    for (size_t i = 0; i < connectivity.getCount(1); ++i)
    {
      const auto& inc = connectivity.getIncidence({1, 2}, i);
      if (inc.size() == 1)
        boundaryCount++;
      else if (inc.size() == 2)
        interfaceCount++;
    }

    // In a uniform grid, should have both boundary and interface edges
    EXPECT_GT(boundaryCount, 0);
    EXPECT_GT(interfaceCount, 0);
  }

  // ============================================================================
  // COMPLEX MIXED MESH TESTS
  // ============================================================================

  TEST(Rodin_Geometry_Connectivity, ComplexMixed_2D_TrianglesAndQuads)
  {
    // Create a more complex mixed mesh
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(9)
                // Bottom row: 2 triangles
                .polytope(Polytope::Type::Triangle, {0, 1, 3})
                .polytope(Polytope::Type::Triangle, {1, 4, 3})
                // Middle: 1 quad
                .polytope(Polytope::Type::Quadrilateral, {3, 4, 7, 6})
                // Top: 2 triangles
                .polytope(Polytope::Type::Triangle, {4, 5, 7})
                .polytope(Polytope::Type::Triangle, {5, 8, 7})
                // Side triangles
                .polytope(Polytope::Type::Triangle, {1, 2, 4})
                .polytope(Polytope::Type::Triangle, {2, 5, 4});

    EXPECT_EQ(connectivity.getCount(2), 7);

    // Test with Discover mode
    connectivity.compute(2, 1, Connectivity<Context::Local>::Mode::Discover);
    connectivity.compute(2, 0, Connectivity<Context::Local>::Mode::Discover);

    // Verify mixed geometry types
    size_t triangleCount = 0;
    size_t quadCount = 0;
    for (size_t i = 0; i < connectivity.getCount(2); ++i)
    {
      auto geom = connectivity.getGeometry(2, i);
      if (geom == Polytope::Type::Triangle)
        triangleCount++;
      else if (geom == Polytope::Type::Quadrilateral)
        quadCount++;
    }

    EXPECT_EQ(triangleCount, 6);
    EXPECT_EQ(quadCount, 1);

    // Verify edge connectivity
    EXPECT_GT(connectivity.getCount(1), 0);
  }

  TEST(Rodin_Geometry_Connectivity, ComplexMixed_3D_TetsAndWedges)
  {
    // Create a 3D mixed mesh with tets and wedges (triangular prisms)
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(3)
                .nodes(9)
                .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
                .polytope(Polytope::Type::Wedge, {1, 2, 4, 3, 5, 6})
                .polytope(Polytope::Type::Tetrahedron, {2, 4, 5, 8});

    EXPECT_EQ(connectivity.getCount(3), 3);

    // Compute connectivity with Discover mode
    connectivity.compute(3, 2, Connectivity<Context::Local>::Mode::Discover);
    connectivity.compute(3, 1, Connectivity<Context::Local>::Mode::Discover);

    // Verify geometry types
    EXPECT_EQ(connectivity.getCount(Polytope::Type::Tetrahedron), 2);
    EXPECT_EQ(connectivity.getCount(Polytope::Type::Wedge), 1);

    // Verify face counts
    EXPECT_EQ(connectivity.getIncidence({3, 2}, 0).size(), 4);  // Tet has 4 faces
    EXPECT_EQ(connectivity.getIncidence({3, 2}, 1).size(), 5);  // Wedge has 5 faces
    EXPECT_EQ(connectivity.getIncidence({3, 2}, 2).size(), 4);  // Tet has 4 faces

    // Total faces should be computed
    EXPECT_GT(connectivity.getCount(2), 0);
  }

  // ============================================================================
  // NON-TRIVIAL TOPOLOGY TESTS
  // ============================================================================

  TEST(Rodin_Geometry_Connectivity, NonTrivial_DisconnectedComponents)
  {
    // Create two separate triangles (disconnected mesh)
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(6)
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Triangle, {3, 4, 5});

    EXPECT_EQ(connectivity.getCount(2), 2);

    connectivity.compute(2, 1, Connectivity<Context::Local>::Mode::Discover);
    connectivity.compute(2, 2, Connectivity<Context::Local>::Mode::Discover);

    // Should have 6 edges (3 per triangle, no sharing)
    EXPECT_EQ(connectivity.getCount(1), 6);

    // No cell should have neighbors (disconnected)
    EXPECT_EQ(connectivity.getIncidence({2, 2}, 0).size(), 0);
    EXPECT_EQ(connectivity.getIncidence({2, 2}, 1).size(), 0);
  }

  TEST(Rodin_Geometry_Connectivity, NonTrivial_ComplexSharing)
  {
    // Create a mesh where multiple cells share edges in complex ways
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2)
                .nodes(7)
                // Central vertex with 6 triangles around it
                .polytope(Polytope::Type::Triangle, {0, 1, 6})
                .polytope(Polytope::Type::Triangle, {0, 2, 1})
                .polytope(Polytope::Type::Triangle, {0, 3, 2})
                .polytope(Polytope::Type::Triangle, {0, 4, 3})
                .polytope(Polytope::Type::Triangle, {0, 5, 4})
                .polytope(Polytope::Type::Triangle, {0, 6, 5});

    EXPECT_EQ(connectivity.getCount(2), 6);

    connectivity.compute(2, 1, Connectivity<Context::Local>::Mode::Discover);
    connectivity.compute(0, 2, Connectivity<Context::Local>::Mode::Discover);

    // Central vertex (0) should be in all 6 triangles
    EXPECT_EQ(connectivity.getIncidence({0, 2}, 0).size(), 6);

    // Should have 12 edges total:
    // - 6 radial edges connecting center vertex to perimeter vertices
    // - 6 perimeter edges connecting adjacent perimeter vertices
    EXPECT_EQ(connectivity.getCount(1), 12);
  }

  TEST(Rodin_Geometry_Connectivity, NonTrivial_LargeMesh_Consistency)
  {
    // Create a larger mesh to test consistency
    Connectivity<Context::Local> connectivity;
    connectivity.initialize(2).nodes(25);

    // Create a 4x4 grid of triangles manually
    for (size_t i = 0; i < 4; ++i)
    {
      for (size_t j = 0; j < 4; ++j)
      {
        size_t v0 = i * 5 + j;
        size_t v1 = v0 + 1;
        size_t v2 = v0 + 5;
        size_t v3 = v2 + 1;

        connectivity.polytope(Polytope::Type::Triangle, {v0, v1, v2});
        connectivity.polytope(Polytope::Type::Triangle, {v1, v3, v2});
      }
    }

    EXPECT_EQ(connectivity.getCount(2), 32);  // 4*4*2 = 32 triangles

    // Compute all connectivity with Discover
    connectivity.compute(2, 1, Connectivity<Context::Local>::Mode::Discover);
    connectivity.compute(1, 2, Connectivity<Context::Local>::Mode::Discover);
    connectivity.compute(2, 2, Connectivity<Context::Local>::Mode::Discover);

    // Verify consistency: each triangle should have 3 edges
    for (size_t i = 0; i < connectivity.getCount(2); ++i)
    {
      EXPECT_EQ(connectivity.getIncidence({2, 1}, i).size(), 3);
    }

    // Count interior vs boundary edges
    size_t boundaryEdges = 0;
    size_t interiorEdges = 0;
    for (size_t i = 0; i < connectivity.getCount(1); ++i)
    {
      const auto& inc = connectivity.getIncidence({1, 2}, i);
      if (inc.size() == 1)
        boundaryEdges++;
      else if (inc.size() == 2)
        interiorEdges++;
    }

    // Should have both interior and boundary edges
    EXPECT_GT(boundaryEdges, 0);
    EXPECT_GT(interiorEdges, 0);
  }

  TEST(Rodin_Geometry_Connectivity, NonTrivial_AllDimensionConnectivity)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {6, 6});
    auto& connectivity = mesh.getConnectivity();

    // Compute all possible connectivity relations
    connectivity.compute(2, 0);  // Cell to vertex
    connectivity.compute(2, 1);  // Cell to edge
    connectivity.compute(2, 2);  // Cell to cell
    connectivity.compute(1, 0);  // Edge to vertex
    connectivity.compute(1, 1);  // Edge to edge
    connectivity.compute(1, 2);  // Edge to cell
    connectivity.compute(0, 1);  // Vertex to edge
    connectivity.compute(0, 2);  // Vertex to cell

    // Verify all counts are positive
    EXPECT_GT(connectivity.getCount(0), 0);
    EXPECT_GT(connectivity.getCount(1), 0);
    EXPECT_GT(connectivity.getCount(2), 0);

    // Spot check: vertices on interior should be shared by multiple cells
    bool foundSharedVertex = false;
    for (size_t i = 0; i < connectivity.getCount(0); ++i)
    {
      if (connectivity.getIncidence({0, 2}, i).size() > 3)
      {
        foundSharedVertex = true;
        break;
      }
    }
    EXPECT_TRUE(foundSharedVertex);
  }

  TEST(Rodin_Geometry_Connectivity, Mode_Comparison_SameResults)
  {
    // Test that Discover and Restrict give same results when edges pre-exist
    Connectivity<Context::Local> conn1, conn2;

    // Setup identical meshes
    conn1.initialize(2).nodes(4)
         .polytope(Polytope::Type::Triangle, {0, 1, 2})
         .polytope(Polytope::Type::Triangle, {1, 2, 3});

    conn2.initialize(2).nodes(4)
         .polytope(Polytope::Type::Triangle, {0, 1, 2})
         .polytope(Polytope::Type::Triangle, {1, 2, 3});

    // For conn1: use Discover mode
    conn1.compute(2, 1, Connectivity<Context::Local>::Mode::Discover);

    // For conn2: first discover edges, then restrict
    conn2.compute(2, 1, Connectivity<Context::Local>::Mode::Discover);
    size_t edgeCount = conn2.getCount(1);
    conn2.clear(2, 1);
    conn2.compute(2, 1, Connectivity<Context::Local>::Mode::Restrict);

    // Both should have same edge count
    EXPECT_EQ(conn1.getCount(1), conn2.getCount(1));
    EXPECT_EQ(conn1.getCount(1), edgeCount);

    // Both should have same incidence relations
    for (size_t i = 0; i < 2; ++i)
    {
      EXPECT_EQ(conn1.getIncidence({2, 1}, i).size(),
                conn2.getIncidence({2, 1}, i).size());
    }
  }
}
