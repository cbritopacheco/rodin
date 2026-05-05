/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SubMeshGridFunctionTest.cpp
 * @brief Unit tests for GridFunction and FES on SubMesh.
 *
 * Tests cover:
 *  - P1/P0 FES construction on a SubMesh (child mesh as support)
 *  - GridFunction projection onto a FES whose underlying mesh is a SubMesh
 *  - GridFunction evaluation via the restriction path:
 *      FES mesh is a SubMesh, point comes from an ancestor (parent) mesh
 *  - GridFunction evaluation via the inclusion path:
 *      FES mesh is the parent mesh, point comes from a child SubMesh
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  // ---- Helper: build a 4x4 triangle mesh with first-half cells labelled 2 ----

  static Mesh<Context::Local> buildSplitMesh()
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t halfCells = mesh.getCellCount() / 2;
    for (size_t i = 0; i < halfCells; ++i)
      mesh.setAttribute({mesh.getDimension(), static_cast<Index>(i)}, 2);
    return mesh;
  }

  // ==========================================================================
  // Group 1 — P1 FES on SubMesh (child mesh as support)
  // ==========================================================================

  TEST(Rodin_Variational_P1_SubMesh, FES_Construction)
  {
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    P1 fes(sub);

    // P1 DOFs = number of vertices in the SubMesh
    EXPECT_EQ(fes.getSize(), sub.getVertexCount());
    EXPECT_EQ(fes.getVectorDimension(), 1u);
    EXPECT_EQ(&fes.getMesh(), &static_cast<const Mesh<Context::Local>&>(sub));
  }

  TEST(Rodin_Variational_P1_SubMesh, GridFunction_ProjectConstant)
  {
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    P1 fes(sub);
    GridFunction gf(fes);

    gf = RealFunction(7.0);

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); ++i)
      EXPECT_NEAR(gf[i], 7.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_P1_SubMesh, GridFunction_ProjectLinear)
  {
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    P1 fes(sub);
    GridFunction gf(fes);

    // f(x, y) = x + y  is exactly representable by P1
    RealFunction f([](const Point& p) { return p.x() + p.y(); });
    gf = f;

    // Verify each DOF matches the vertex coordinate
    const size_t vdim = sub.getSpaceDimension();
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); ++i)
    {
      const auto& coords = sub.getVertices();
      const Real expected = coords(0, i) + (vdim >= 2 ? coords(1, i) : 0.0);
      EXPECT_NEAR(gf[i], expected, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_P1_SubMesh, GridFunction_EvaluateAtLocalPoint)
  {
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    P1 fes(sub);
    GridFunction gf(fes);
    gf = RealFunction(3.5);

    // Evaluate at the centroid of the first cell of the SubMesh (local point)
    auto it = sub.getPolytope(sub.getDimension(), 0);
    const auto& polytope = *it;
    const Math::SpatialPoint rc{{1.0 / 3.0, 1.0 / 3.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(gf.getValue(p), 3.5, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_P1_SubMesh, GridFunction_EvaluateAtParentPoint_Restriction)
  {
    // The FES mesh is a SubMesh; we evaluate using a Point whose polytope
    // lives in the parent mesh — this exercises the restriction path in
    // GridFunctionBase::getValue().
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    P1 fes(sub);
    GridFunction gf(fes);
    RealFunction f([](const Point& p) { return p.x() + p.y(); });
    gf = f;

    // Find the first parent-mesh cell that was kept (attribute 2).
    const size_t mdim = mesh.getDimension();
    const auto& pmap = sub.getPolytopeMap(mdim);

    // pmap.left[0] is the parent-mesh index of sub cell 0
    const Index parentCellIdx = pmap.left[0];

    auto it = mesh.getPolytope(mdim, parentCellIdx);
    const auto& parentPolytope = *it;
    const Math::SpatialPoint rc{{1.0 / 3.0, 1.0 / 3.0}};
    // Build a point whose polytope belongs to the *parent* mesh
    Point parentPoint(parentPolytope, rc);

    // Transform reference coordinates to physical coordinates to get expected
    const auto& trans = mesh.getPolytopeTransformation(mdim, parentCellIdx);
    Math::SpatialPoint pc;
    trans.transform(pc, rc);
    const Real expected = pc(0) + pc(1);

    EXPECT_NEAR(gf.getValue(parentPoint), expected, RODIN_FUZZY_CONSTANT);
  }

  // ==========================================================================
  // Group 2 — P1 FES on parent mesh, evaluated at SubMesh points (inclusion)
  // ==========================================================================

  TEST(Rodin_Variational_P1_ParentMesh, GridFunction_EvaluateAtSubMeshPoint_Inclusion)
  {
    // The FES mesh is the full parent mesh; we evaluate using a Point whose
    // polytope lives in a child SubMesh — this exercises the inclusion path
    // in GridFunctionBase::getValue().
    Mesh mesh = buildSplitMesh();

    P1 fes(mesh);
    GridFunction gf(fes);
    RealFunction f([](const Point& p) { return p.x() + p.y(); });
    gf = f;

    // Build a SubMesh over the "kept" half
    SubMesh<Context::Local> sub = mesh.keep(2);

    const size_t mdim = sub.getDimension();
    // Take the first cell from the SubMesh
    auto it = sub.getPolytope(mdim, 0);
    const auto& subPolytope = *it;
    const Math::SpatialPoint rc{{1.0 / 3.0, 1.0 / 3.0}};
    // Build a point whose polytope belongs to the *SubMesh*
    Point subPoint(subPolytope, rc);

    // The expected value is f evaluated at the physical coordinates of subPoint.
    // Since f(x,y)=x+y is linear and P1 represents it exactly, the interpolated
    // value must equal f at the physical location of the point.
    const Real expected = subPoint.x() + subPoint.y();

    EXPECT_NEAR(gf.getValue(subPoint), expected, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_P1_ParentMesh, GridFunction_EvaluateAtBoundarySubMeshPoint_Inclusion)
  {
    // Skin SubMesh: (d-1)-dimensional boundary of the full 2D mesh.
    // The FES is on the parent 2D mesh; evaluate at a boundary edge point.
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    GridFunction gf(fes);
    RealFunction f([](const Point& p) { return p.x() + p.y(); });
    gf = f;

    SubMesh<Context::Local> skin = mesh.skin();

    const size_t edgeDim = skin.getDimension();   // = 1 for 2D mesh
    auto it = skin.getPolytope(edgeDim, 0);
    const auto& edgePolytope = *it;
    // Midpoint of the segment in reference coordinates
    const Math::SpatialPoint rc{{0.5}};
    Point skinPoint(edgePolytope, rc);

    // The expected value is f at the physical location of skinPoint.
    // Since f(x,y)=x+y is linear and P1 represents it exactly, the interpolated
    // value must equal f at the physical midpoint of the boundary edge.
    const Real expected = skinPoint.x() + skinPoint.y();

    EXPECT_NEAR(gf.getValue(skinPoint), expected, RODIN_FUZZY_CONSTANT);
  }

  // ==========================================================================
  // Group 3 — P0 FES on SubMesh (child mesh as support)
  // ==========================================================================

  TEST(Rodin_Variational_P0_SubMesh, FES_Construction)
  {
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    P0 fes(sub);

    // P0 DOFs = one per cell in the SubMesh
    EXPECT_EQ(fes.getSize(), sub.getCellCount());
    EXPECT_EQ(fes.getVectorDimension(), 1u);
    EXPECT_EQ(&fes.getMesh(), &static_cast<const Mesh<Context::Local>&>(sub));
  }

  TEST(Rodin_Variational_P0_SubMesh, GridFunction_ProjectConstant)
  {
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    P0 fes(sub);
    GridFunction gf(fes);

    gf = RealFunction(4.2);

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); ++i)
      EXPECT_NEAR(gf[i], 4.2, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_P0_SubMesh, GridFunction_ProjectSpatialFunction)
  {
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    P0 fes(sub);
    GridFunction gf(fes);

    // Project f = 1 (constant function) — piecewise-constant projection
    RealFunction f([](const Point&) { return 1.0; });
    gf = f;

    // All DOFs must equal 1.0 since f is constant
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); ++i)
      EXPECT_NEAR(gf[i], 1.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_P0_SubMesh, GridFunction_EvaluateAtParentPoint_Restriction)
  {
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    P0 fes(sub);
    GridFunction gf(fes);
    gf = RealFunction(9.0);

    // Evaluate at a parent-mesh cell that belongs to the SubMesh
    const size_t mdim = mesh.getDimension();
    const Index parentCellIdx = sub.getPolytopeMap(mdim).left[0];

    auto it = mesh.getPolytope(mdim, parentCellIdx);
    const auto& parentPolytope = *it;
    const Math::SpatialPoint rc{{1.0 / 3.0, 1.0 / 3.0}};
    Point parentPoint(parentPolytope, rc);

    EXPECT_NEAR(gf.getValue(parentPoint), 9.0, RODIN_FUZZY_CONSTANT);
  }

  // ==========================================================================
  // Group 4 — Vector-valued P1 FES on SubMesh
  // ==========================================================================

  TEST(Rodin_Variational_Vector_P1_SubMesh, FES_Construction)
  {
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    constexpr size_t vdim = 2;
    P1 fes(sub, vdim);

    EXPECT_EQ(fes.getSize(), vdim * sub.getVertexCount());
    EXPECT_EQ(fes.getVectorDimension(), vdim);
  }

  TEST(Rodin_Variational_Vector_P1_SubMesh, GridFunction_ProjectVectorConstant)
  {
    Mesh mesh = buildSplitMesh();
    SubMesh<Context::Local> sub = mesh.keep(2);

    constexpr size_t vdim = 2;
    P1 fes(sub, vdim);
    GridFunction gf(fes);

    VectorFunction v{1.0, 2.0};
    gf = v;

    // DOF layout is component-major:
    //   DOFs [0 .. nv)    → component 0 (value = 1.0)
    //   DOFs [nv .. 2*nv) → component 1 (value = 2.0)
    const size_t nv = sub.getVertexCount();
    for (Index i = 0; i < static_cast<Index>(nv); ++i)
      EXPECT_NEAR(gf[i], 1.0, RODIN_FUZZY_CONSTANT);
    for (Index i = static_cast<Index>(nv); i < static_cast<Index>(2 * nv); ++i)
      EXPECT_NEAR(gf[i], 2.0, RODIN_FUZZY_CONSTANT);
  }
}
