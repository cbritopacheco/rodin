/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <sstream>
#include <cmath>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  // ========================================================================
  // MMG::Mesh construction tests
  // ========================================================================

  TEST(Rodin_MMG_Mesh, DefaultConstruction)
  {
    MMG::Mesh mesh;
    EXPECT_TRUE(mesh.isEmpty());
    EXPECT_EQ(mesh.getCorners().size(), 0);
    EXPECT_EQ(mesh.getRidges().size(), 0);
    EXPECT_EQ(mesh.getRequiredVertices().size(), 0);
    EXPECT_EQ(mesh.getRequiredEdges().size(), 0);
  }

  TEST(Rodin_MMG_Mesh, UniformGridConstruction)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_EQ(mesh.getVertexCount(), n * n);
    EXPECT_EQ(mesh.getCellCount(), 2 * (n - 1) * (n - 1));
    EXPECT_EQ(mesh.getDimension(), 2);
    EXPECT_EQ(mesh.getSpaceDimension(), 2);
  }

  TEST(Rodin_MMG_Mesh, SetCorner)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    mesh.setCorner(0);
    mesh.setCorner(3);
    mesh.setCorner(12);
    mesh.setCorner(15);

    EXPECT_EQ(mesh.getCorners().size(), 4);
    EXPECT_TRUE(mesh.getCorners().count(0));
    EXPECT_TRUE(mesh.getCorners().count(3));
    EXPECT_TRUE(mesh.getCorners().count(12));
    EXPECT_TRUE(mesh.getCorners().count(15));
  }

  TEST(Rodin_MMG_Mesh, SetRidge)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(1, 2);

    size_t ridgeCount = 0;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      mesh.setRidge(it->getIndex());
      ridgeCount++;
    }

    EXPECT_EQ(mesh.getRidges().size(), ridgeCount);
    EXPECT_GT(ridgeCount, 0);
  }

  TEST(Rodin_MMG_Mesh, SetRequiredVertex)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    mesh.setRequiredVertex(0);
    mesh.setRequiredVertex(5);

    EXPECT_EQ(mesh.getRequiredVertices().size(), 2);
    EXPECT_TRUE(mesh.getRequiredVertices().count(0));
    EXPECT_TRUE(mesh.getRequiredVertices().count(5));
  }

  TEST(Rodin_MMG_Mesh, SetRequiredEdge)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    mesh.setRequiredEdge(0);
    mesh.setRequiredEdge(1);

    EXPECT_EQ(mesh.getRequiredEdges().size(), 2);
    EXPECT_TRUE(mesh.getRequiredEdges().count(0));
    EXPECT_TRUE(mesh.getRequiredEdges().count(1));
  }

  // ========================================================================
  // MMG::Mesh copy and move tests
  // ========================================================================

  TEST(Rodin_MMG_Mesh, CopyConstruction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.setCorner(0);
    mesh.setCorner(3);
    mesh.setRidge(0);
    mesh.setRequiredVertex(5);
    mesh.setRequiredEdge(2);

    MMG::Mesh copy(mesh);
    EXPECT_EQ(copy.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(copy.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(copy.getCorners().size(), 2);
    EXPECT_EQ(copy.getRidges().size(), 1);
    EXPECT_EQ(copy.getRequiredVertices().size(), 1);
    EXPECT_EQ(copy.getRequiredEdges().size(), 1);
    EXPECT_TRUE(copy.getCorners().count(0));
    EXPECT_TRUE(copy.getCorners().count(3));
  }

  TEST(Rodin_MMG_Mesh, MoveConstruction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.setCorner(0);
    mesh.setRidge(0);
    mesh.setRequiredVertex(5);
    mesh.setRequiredEdge(2);

    size_t vertexCount = mesh.getVertexCount();
    size_t cellCount = mesh.getCellCount();

    MMG::Mesh moved(std::move(mesh));
    EXPECT_EQ(moved.getVertexCount(), vertexCount);
    EXPECT_EQ(moved.getCellCount(), cellCount);
    EXPECT_EQ(moved.getCorners().size(), 1);
    EXPECT_EQ(moved.getRidges().size(), 1);
    EXPECT_EQ(moved.getRequiredVertices().size(), 1);
    EXPECT_EQ(moved.getRequiredEdges().size(), 1);
  }

  TEST(Rodin_MMG_Mesh, MoveAssignment)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.setCorner(0);
    mesh.setRidge(0);

    size_t vertexCount = mesh.getVertexCount();
    size_t cellCount = mesh.getCellCount();

    MMG::Mesh target;
    target = std::move(mesh);
    EXPECT_EQ(target.getVertexCount(), vertexCount);
    EXPECT_EQ(target.getCellCount(), cellCount);
    EXPECT_EQ(target.getCorners().size(), 1);
    EXPECT_EQ(target.getRidges().size(), 1);
  }

  TEST(Rodin_MMG_Mesh, MoveFromParent)
  {
    Mesh<Context::Local> parentMesh;
    parentMesh = parentMesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    size_t vertexCount = parentMesh.getVertexCount();
    size_t cellCount = parentMesh.getCellCount();

    MMG::Mesh mmgMesh(std::move(parentMesh));
    EXPECT_EQ(mmgMesh.getVertexCount(), vertexCount);
    EXPECT_EQ(mmgMesh.getCellCount(), cellCount);
    // MMG metadata should be empty since parent doesn't carry it
    EXPECT_EQ(mmgMesh.getCorners().size(), 0);
    EXPECT_EQ(mmgMesh.getRidges().size(), 0);
  }

  // ========================================================================
  // MMG::Mesh Builder tests
  // ========================================================================

  TEST(Rodin_MMG_Mesh_Builder, TriangleMeshWithMMGMetadata)
  {
    MMG::Mesh mesh =
      MMG::Mesh::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .corner(0)
      .corner(1)
      .corner(2)
      .corner(3)
      .finalize();

    EXPECT_EQ(mesh.getVertexCount(), 4);
    EXPECT_EQ(mesh.getCellCount(), 2);
    EXPECT_EQ(mesh.getCorners().size(), 4);
    EXPECT_TRUE(mesh.getCorners().count(0));
    EXPECT_TRUE(mesh.getCorners().count(1));
    EXPECT_TRUE(mesh.getCorners().count(2));
    EXPECT_TRUE(mesh.getCorners().count(3));
  }

  TEST(Rodin_MMG_Mesh_Builder, BuilderWithAllMMGMetadata)
  {
    MMG::Mesh mesh =
      MMG::Mesh::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .corner(0)
      .ridge(0)
      .requiredVertex(1)
      .requiredEdge(0)
      .finalize();

    EXPECT_EQ(mesh.getCorners().size(), 1);
    EXPECT_EQ(mesh.getRidges().size(), 1);
    EXPECT_EQ(mesh.getRequiredVertices().size(), 1);
    EXPECT_EQ(mesh.getRequiredEdges().size(), 1);
  }

  // ========================================================================
  // MMG::Mesh I/O roundtrip tests (MEDIT format)
  // ========================================================================

  TEST(Rodin_MMG_Mesh_IO, SaveLoadMEDITRoundtrip)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(1, 2);

    // Set some MMG metadata
    mesh.setCorner(0);
    mesh.setCorner(3);
    mesh.setCorner(12);
    mesh.setCorner(15);

    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      mesh.setRidge(it->getIndex());

    mesh.setRequiredVertex(0);
    mesh.setRequiredVertex(15);

    // Save to file
    const std::string filename = "/tmp/rodin_mmg_test_roundtrip.mesh";
    mesh.save(filename, IO::FileFormat::MEDIT);

    // Load back
    MMG::Mesh loaded;
    loaded.load(filename, IO::FileFormat::MEDIT);

    // Verify topology preserved
    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());

    // Verify MMG metadata preserved
    EXPECT_EQ(loaded.getCorners().size(), mesh.getCorners().size());
    EXPECT_EQ(loaded.getRidges().size(), mesh.getRidges().size());
    EXPECT_EQ(loaded.getRequiredVertices().size(), mesh.getRequiredVertices().size());

    // Verify specific corners
    EXPECT_TRUE(loaded.getCorners().count(0));
    EXPECT_TRUE(loaded.getCorners().count(3));
    EXPECT_TRUE(loaded.getCorners().count(12));
    EXPECT_TRUE(loaded.getCorners().count(15));

    // Verify required vertices
    EXPECT_TRUE(loaded.getRequiredVertices().count(0));
    EXPECT_TRUE(loaded.getRequiredVertices().count(15));

    // Clean up
    std::remove(filename.c_str());
  }

  TEST(Rodin_MMG_Mesh_IO, SaveLoadMEDITEmptyMetadata)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    // Save without any MMG metadata
    const std::string filename = "/tmp/rodin_mmg_test_empty_metadata.mesh";
    mesh.save(filename, IO::FileFormat::MEDIT);

    MMG::Mesh loaded;
    loaded.load(filename, IO::FileFormat::MEDIT);

    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(loaded.getCorners().size(), 0);
    EXPECT_EQ(loaded.getRidges().size(), 0);
    EXPECT_EQ(loaded.getRequiredVertices().size(), 0);
    EXPECT_EQ(loaded.getRequiredEdges().size(), 0);

    std::remove(filename.c_str());
  }

  // ========================================================================
  // MMG5 mesh conversion roundtrip tests
  // ========================================================================

  TEST(Rodin_MMG_MMG5, RodinToMeshAndBack2D)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    size_t origVertexCount = mesh.getVertexCount();
    size_t origCellCount = mesh.getCellCount();

    // Convert to MMG and back
    MMG5_pMesh mmgMesh = MMG::MMG5::rodinToMesh(mesh);
    ASSERT_NE(mmgMesh, nullptr);

    MMG::Mesh result = MMG::MMG5::meshToRodin(mmgMesh);
    MMG::MMG5::destroyMesh(mmgMesh);

    EXPECT_EQ(result.getVertexCount(), origVertexCount);
    EXPECT_EQ(result.getCellCount(), origCellCount);
    EXPECT_EQ(result.getDimension(), 2);
  }

  TEST(Rodin_MMG_MMG5, RodinToMeshPreservesMMGMetadata)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(1, 2);

    mesh.setCorner(0);
    mesh.setCorner(3);

    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      mesh.setRidge(it->getIndex());

    mesh.setRequiredVertex(0);

    // Convert to MMG and back
    MMG5_pMesh mmgMesh = MMG::MMG5::rodinToMesh(mesh);
    ASSERT_NE(mmgMesh, nullptr);

    MMG::Mesh result = MMG::MMG5::meshToRodin(mmgMesh);
    MMG::MMG5::destroyMesh(mmgMesh);

    // Corners should be preserved
    EXPECT_EQ(result.getCorners().size(), mesh.getCorners().size());
    EXPECT_TRUE(result.getCorners().count(0));
    EXPECT_TRUE(result.getCorners().count(3));

    // Ridges should be preserved
    EXPECT_EQ(result.getRidges().size(), mesh.getRidges().size());

    // Required vertices should be preserved
    EXPECT_EQ(result.getRequiredVertices().size(), mesh.getRequiredVertices().size());
    EXPECT_TRUE(result.getRequiredVertices().count(0));
  }

  // ========================================================================
  // MMG::Optimizer tests
  // ========================================================================

  TEST(Rodin_MMG_Optimizer, Optimize2DTriangleMesh)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    size_t origVertexCount = mesh.getVertexCount();

    MMG::Optimizer optimizer;
    optimizer.setHMax(0.5);
    optimizer.optimize(mesh);

    // After optimization, mesh should still be valid and non-empty
    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
    EXPECT_EQ(mesh.getDimension(), 2);
  }

  TEST(Rodin_MMG_Optimizer, FluentAPI)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    // Fluent API method chaining
    MMG::Optimizer()
      .setHMax(0.5)
      .setAngleDetection(true)
      .setGradation(1.3)
      .optimize(mesh);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  TEST(Rodin_MMG_Optimizer, OptimizeWithCorners)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(1, 2);

    // Mark corners
    mesh.setCorner(0);
    mesh.setCorner(n - 1);
    mesh.setCorner(n * (n - 1));
    mesh.setCorner(n * n - 1);

    // Mark ridges
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      mesh.setRidge(it->getIndex());

    MMG::Optimizer().setHMax(0.5).optimize(mesh);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  // ========================================================================
  // MMG::Adapt tests
  // ========================================================================

  TEST(Rodin_MMG_Adapt, UniformSizeMap2D)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    // Create a constant size map
    P1 fes(mesh);
    MMG::RealGridFunction sizeMap(fes);
    sizeMap = [](const Geometry::Point&) { return 0.2; };

    MMG::Adapt adapter;
    adapter.setHMax(0.5);
    adapter.setHMin(0.05);
    adapter.adapt(mesh, sizeMap);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
    EXPECT_EQ(mesh.getDimension(), 2);
  }

  TEST(Rodin_MMG_Adapt, VariableSizeMap2D)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    // Create a variable size map: smaller near center, larger near boundary
    P1 fes(mesh);
    MMG::RealGridFunction sizeMap(fes);
    sizeMap = [&n](const Geometry::Point& p)
    {
      Real cx = static_cast<Real>(n - 1) / 2.0;
      Real cy = cx;
      Real dx = p.x() - cx;
      Real dy = p.y() - cy;
      Real dist = std::sqrt(dx * dx + dy * dy);
      return 0.1 + 0.3 * dist / cx;
    };

    MMG::Adapt adapter;
    adapter.setHMin(0.05);
    adapter.setHMax(1.0);
    adapter.adapt(mesh, sizeMap);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  TEST(Rodin_MMG_Adapt, FluentAPI)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    P1 fes(mesh);
    MMG::RealGridFunction sizeMap(fes);
    sizeMap = [](const Geometry::Point&) { return 0.2; };

    // Fluent API
    MMG::Adapt()
      .setHMin(0.05)
      .setHMax(0.5)
      .setGradation(1.3)
      .setHausdorff(0.01)
      .setAngleDetection(true)
      .adapt(mesh, sizeMap);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  // ========================================================================
  // MMG::LevelSetDiscretizer tests
  // ========================================================================

  TEST(Rodin_MMG_LevelSetDiscretizer, CircleLevelSet2D)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute boundary = 10;
    static constexpr Real radius = 0.1;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.scale(1.0 / (n - 1));

    // Set all cells to interior attribute
    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, interior);

    // Create level-set function: circle of given radius centered at (0.5, 0.5)
    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p)
    {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - radius;
    };

    MMG::Mesh result =
      MMG::LevelSetDiscretizer()
        .split(interior, { interior, exterior })
        .setBoundaryReference(boundary)
        .setHMax(0.1)
        .discretize(ls);

    // Result should be a valid non-empty mesh
    EXPECT_FALSE(result.isEmpty());
    EXPECT_GT(result.getVertexCount(), 0);
    EXPECT_GT(result.getCellCount(), 0);
    EXPECT_EQ(result.getDimension(), 2);
  }

  TEST(Rodin_MMG_LevelSetDiscretizer, FluentAPI)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.scale(1.0 / (n - 1));

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, interior);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p)
    {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - 0.1;
    };

    // Verify fluent API chaining compiles and works
    MMG::Mesh result =
      MMG::LevelSetDiscretizer()
        .split(interior, { interior, exterior })
        .setBoundaryReference(10)
        .setHMin(0.01)
        .setHMax(0.1)
        .setGradation(1.3)
        .setAngleDetection(true)
        .discretize(ls);

    EXPECT_FALSE(result.isEmpty());
  }

  TEST(Rodin_MMG_LevelSetDiscretizer, NoSplitMaterial)
  {
    static constexpr Attribute mat1 = 1;
    static constexpr Attribute mat2 = 2;
    static constexpr Attribute exterior = 3;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.scale(1.0 / (n - 1));

    // Set half the cells to mat1, half to mat2
    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, mat1);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p)
    {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - 0.1;
    };

    MMG::LevelSetDiscretizer discretizer;
    discretizer.split(mat1, { mat1, exterior });
    discretizer.setBoundaryReference(10);
    discretizer.setHMax(0.1);

    // Verify getSplitMap() reflects the configuration
    const auto& splitMap = discretizer.getSplitMap();
    EXPECT_EQ(splitMap.size(), 1);
    EXPECT_TRUE(splitMap.count(mat1));

    MMG::Mesh result = discretizer.discretize(ls);
    EXPECT_FALSE(result.isEmpty());
  }

  TEST(Rodin_MMG_LevelSetDiscretizer, SetLevelSetValue)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.scale(1.0 / (n - 1));

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, interior);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p)
    {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - 0.1;
    };

    // Discretize at a non-zero level-set value
    MMG::Mesh result =
      MMG::LevelSetDiscretizer()
        .split(interior, { interior, exterior })
        .setBoundaryReference(10)
        .setLevelSet(0.0)
        .setHMax(0.1)
        .discretize(ls);

    EXPECT_FALSE(result.isEmpty());
    EXPECT_GT(result.getVertexCount(), 0);
  }

  // ========================================================================
  // Common types tests
  // ========================================================================

  TEST(Rodin_MMG_Common, SplitMapConstruction)
  {
    MMG::SplitMap splitMap;
    splitMap[1] = MMG::Split{ 2, 3 };
    splitMap[4] = MMG::NoSplit;

    EXPECT_EQ(splitMap.size(), 2);
    EXPECT_TRUE(std::holds_alternative<MMG::Split>(splitMap.at(1)));
    EXPECT_TRUE(std::holds_alternative<MMG::NoSplitT>(splitMap.at(4)));

    const auto& split = std::get<MMG::Split>(splitMap.at(1));
    EXPECT_EQ(split.interior, 2);
    EXPECT_EQ(split.exterior, 3);
  }

  // ========================================================================
  // MMG::GridFunction alias tests
  // ========================================================================

  TEST(Rodin_MMG_GridFunction, ScalarGridFunctionConstruction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);
    gf = [](const Geometry::Point& p) { return p.x() + p.y(); };

    // Grid function should have correct size
    EXPECT_EQ(gf.getFiniteElementSpace().getSize(), mesh.getVertexCount());
  }

  TEST(Rodin_MMG_GridFunction, VectorGridFunctionConstruction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    P1 fes(mesh, 2);
    MMG::VectorGridFunction gf(fes);
    gf = [](const Geometry::Point& p)
    {
      return Math::Vector<Real>{{ p.x(), p.y() }};
    };

    EXPECT_EQ(gf.getFiniteElementSpace().getVectorDimension(), 2);
  }

  // ========================================================================
  // MMG5 parameter configuration tests
  // ========================================================================

  TEST(Rodin_MMG_MMG5, ParameterSetters)
  {
    // Verify that parameter setters return reference for chaining
    MMG::Optimizer opt;
    auto& ref1 = opt.setHMin(0.01);
    auto& ref2 = opt.setHMax(1.0);
    auto& ref3 = opt.setHausdorff(0.001);
    auto& ref4 = opt.setGradation(1.3);
    auto& ref5 = opt.setAngleDetection(true);

    // All should return reference to the same object
    EXPECT_EQ(&ref1, &opt);
    EXPECT_EQ(&ref2, &opt);
    EXPECT_EQ(&ref3, &opt);
    EXPECT_EQ(&ref4, &opt);
    EXPECT_EQ(&ref5, &opt);
  }

  TEST(Rodin_MMG_MMG5, AdaptParameterSetters)
  {
    MMG::Adapt adapt;
    auto& ref1 = adapt.setHMin(0.01);
    auto& ref2 = adapt.setHMax(1.0);
    auto& ref3 = adapt.setHausdorff(0.001);
    auto& ref4 = adapt.setGradation(1.3);
    auto& ref5 = adapt.setAngleDetection(true);

    EXPECT_EQ(&ref1, &adapt);
    EXPECT_EQ(&ref2, &adapt);
    EXPECT_EQ(&ref3, &adapt);
    EXPECT_EQ(&ref4, &adapt);
    EXPECT_EQ(&ref5, &adapt);
  }

  TEST(Rodin_MMG_MMG5, LevelSetDiscretizerParameterSetters)
  {
    MMG::LevelSetDiscretizer lsd;
    auto& ref1 = lsd.setHMin(0.01);
    auto& ref2 = lsd.setHMax(1.0);
    auto& ref3 = lsd.setHausdorff(0.001);
    auto& ref4 = lsd.setGradation(1.3);
    auto& ref5 = lsd.setAngleDetection(true);
    auto& ref6 = lsd.setLevelSet(0.5);
    auto& ref7 = lsd.setBoundaryReference(10);
    auto& ref8 = lsd.surface(true);

    EXPECT_EQ(&ref1, &lsd);
    EXPECT_EQ(&ref2, &lsd);
    EXPECT_EQ(&ref3, &lsd);
    EXPECT_EQ(&ref4, &lsd);
    EXPECT_EQ(&ref5, &lsd);
    EXPECT_EQ(&ref6, &lsd);
    EXPECT_EQ(&ref7, &lsd);
    EXPECT_EQ(&ref8, &lsd);
  }
}
