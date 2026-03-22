/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/T8Code/Mesh.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ---- Construction -------------------------------------------------------

  TEST(Rodin_T8Code_Mesh, ConstructFromTriangleMesh)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Triangle, {3, 3});

    T8Code::Mesh amr(std::move(base));

    // Should have the same number of cells and vertices as the base mesh
    EXPECT_GT(amr.getCellCount(), 0u);
    EXPECT_GT(amr.getVertexCount(), 0u);
    EXPECT_EQ(amr.getMaxRefinementLevel(), 0u);
  }

  TEST(Rodin_T8Code_Mesh, ConstructFromQuadMesh)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Quadrilateral, {3, 3});

    T8Code::Mesh amr(std::move(base));

    EXPECT_GT(amr.getCellCount(), 0u);
    EXPECT_GT(amr.getVertexCount(), 0u);
    EXPECT_EQ(amr.getMaxRefinementLevel(), 0u);
  }

  // ---- Uniform refinement -------------------------------------------------

  TEST(Rodin_T8Code_Mesh, UniformRefineTriangles)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Triangle, {2, 2});

    const size_t baseCellCount = base.getCellCount();

    T8Code::Mesh amr(std::move(base));
    amr.refine();

    // After one level of uniform refinement, each triangle splits into 4
    EXPECT_EQ(amr.getCellCount(), baseCellCount * 4);
    EXPECT_EQ(amr.getMaxRefinementLevel(), 1u);

    // All cells should be at level 1
    for (Index i = 0; i < static_cast<Index>(amr.getCellCount()); ++i)
      EXPECT_EQ(amr.getRefinementLevel(i), 1u);
  }

  TEST(Rodin_T8Code_Mesh, UniformRefineQuads)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Quadrilateral, {2, 2});

    const size_t baseCellCount = base.getCellCount();

    T8Code::Mesh amr(std::move(base));
    amr.refine();

    // After one level of uniform refinement, each quad splits into 4
    EXPECT_EQ(amr.getCellCount(), baseCellCount * 4);
    EXPECT_EQ(amr.getMaxRefinementLevel(), 1u);
  }

  TEST(Rodin_T8Code_Mesh, TwoLevelsUniformRefinement)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Quadrilateral, {2, 2});

    const size_t baseCellCount = base.getCellCount();

    T8Code::Mesh amr(std::move(base));
    amr.refine();
    amr.refine();

    // After two levels: each quad -> 4 -> 16
    EXPECT_EQ(amr.getCellCount(), baseCellCount * 16);
    EXPECT_EQ(amr.getMaxRefinementLevel(), 2u);
  }

  // ---- Local (predicate-based) refinement ----------------------------------

  TEST(Rodin_T8Code_Mesh, LocalRefinement)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Quadrilateral, {3, 3});

    const size_t baseCellCount = base.getCellCount();

    T8Code::Mesh amr(std::move(base));

    // Refine only cells in the lower-left quadrant
    amr.refine([](const Cell& cell) {
      auto vIt = cell.getVertex();
      auto coords = vIt->getCoordinates();
      return coords[0] < 0.5 && coords[1] < 0.5;
    });

    // Some cells should be refined, so total count > base count
    EXPECT_GT(amr.getCellCount(), baseCellCount);
    // But not all cells are refined, so count < 4 * base count
    EXPECT_LT(amr.getCellCount(), baseCellCount * 4);
    EXPECT_EQ(amr.getMaxRefinementLevel(), 1u);
  }

  // ---- Hanging nodes -------------------------------------------------------

  TEST(Rodin_T8Code_Mesh, HangingNodesAfterLocalRefinement)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Quadrilateral, {3, 3});

    T8Code::Mesh amr(std::move(base));

    // Refine only cells in the lower-left quadrant
    amr.refine([](const Cell& cell) {
      auto vIt = cell.getVertex();
      auto coords = vIt->getCoordinates();
      return coords[0] < 0.5 && coords[1] < 0.5;
    });

    // After local refinement, there should be hanging nodes
    auto hangingNodes = amr.getHangingNodes();
    EXPECT_GT(hangingNodes.size(), 0u);

    // Each hanging node should have valid constraining vertices
    for (Index hn : hangingNodes)
    {
      EXPECT_TRUE(amr.isHangingNode(hn));
      auto [v0, v1] = amr.getConstrainingVertices(hn);
      EXPECT_NE(v0, v1);
      EXPECT_LT(static_cast<size_t>(v0), amr.getVertexCount());
      EXPECT_LT(static_cast<size_t>(v1), amr.getVertexCount());
    }
  }

  TEST(Rodin_T8Code_Mesh, NoHangingNodesAfterUniformRefinement)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Quadrilateral, {3, 3});

    T8Code::Mesh amr(std::move(base));
    amr.refine();

    // Uniform refinement should produce no hanging nodes
    auto hangingNodes = amr.getHangingNodes();
    EXPECT_EQ(hangingNodes.size(), 0u);
  }

  // ---- Balance -------------------------------------------------------------

  TEST(Rodin_T8Code_Mesh, Balance)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Quadrilateral, {3, 3});

    T8Code::Mesh amr(std::move(base));

    // Create a 2-level difference by refining one cell twice
    amr.refine([](const Cell& cell) {
      auto vIt = cell.getVertex();
      auto coords = vIt->getCoordinates();
      return coords[0] < 0.3 && coords[1] < 0.3;
    });
    amr.refine([](const Cell& cell) {
      auto vIt = cell.getVertex();
      auto coords = vIt->getCoordinates();
      return coords[0] < 0.15 && coords[1] < 0.15;
    });

    const size_t preBalanceCellCount = amr.getCellCount();

    // Balance should add cells to enforce 2:1 constraint
    amr.balance();

    // After balance, cell count should be >= pre-balance count
    EXPECT_GE(amr.getCellCount(), preBalanceCellCount);
  }

  // ---- Coarsening ----------------------------------------------------------

  TEST(Rodin_T8Code_Mesh, CoarsenAfterUniformRefine)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Quadrilateral, {2, 2});

    const size_t baseCellCount = base.getCellCount();

    T8Code::Mesh amr(std::move(base));
    amr.refine();

    EXPECT_EQ(amr.getCellCount(), baseCellCount * 4);

    // Coarsen all cells
    amr.coarsen([](const Cell&) { return true; });

    // Should return to original cell count
    EXPECT_EQ(amr.getCellCount(), baseCellCount);
    EXPECT_EQ(amr.getMaxRefinementLevel(), 0u);
  }

  // ---- Copy / Move ---------------------------------------------------------

  TEST(Rodin_T8Code_Mesh, CopyConstruction)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Triangle, {3, 3});

    T8Code::Mesh amr(std::move(base));
    amr.refine();

    T8Code::Mesh copy(amr);

    EXPECT_EQ(copy.getCellCount(), amr.getCellCount());
    EXPECT_EQ(copy.getVertexCount(), amr.getVertexCount());
    EXPECT_EQ(copy.getMaxRefinementLevel(), amr.getMaxRefinementLevel());
  }

  TEST(Rodin_T8Code_Mesh, MoveConstruction)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Triangle, {3, 3});

    T8Code::Mesh amr(std::move(base));
    amr.refine();

    const size_t cellCount = amr.getCellCount();
    const size_t vertexCount = amr.getVertexCount();

    T8Code::Mesh moved(std::move(amr));

    EXPECT_EQ(moved.getCellCount(), cellCount);
    EXPECT_EQ(moved.getVertexCount(), vertexCount);
  }

  // ---- Accessors -----------------------------------------------------------

  TEST(Rodin_T8Code_Mesh, ForestAndCmeshAccessors)
  {
    Mesh base = Mesh<Context::Local>::UniformGrid(
      Polytope::Type::Triangle, {3, 3});

    T8Code::Mesh amr(std::move(base));

    EXPECT_NE(amr.getForest(), nullptr);
    EXPECT_NE(amr.getCoarseMesh(), nullptr);
  }
}
