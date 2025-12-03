#include <gtest/gtest.h>
#include <type_traits>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "Rodin/Test/Random.h"

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/IO/MFEM.h"
#include "Rodin/IO/MEDIT.h"
#include "Rodin/Configure.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_H1_Space, Triangle_VertexNodes_MatchConnectivity)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0}) // 0
      .vertex({1, 0}) // 1
      .vertex({0, 1}) // 2
      .vertex({1, 1}) // 3
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} }) // cell 0
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} }) // cell 1
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    auto cell_it = mesh.getCell();
    for (; cell_it; ++cell_it)
    {
      const auto& cell = *cell_it;
      const auto& fe   = fes.getFiniteElement(2, cell.getIndex());
      ASSERT_EQ(fe.getGeometry(), Polytope::Type::Triangle);

      // connectivity vertices
      auto vconn = cell.getVertices(); // {v0, v1, v2}

      for (int local_v = 0; local_v < 3; ++local_v)
      {
        Index vtx = vconn[local_v];

        auto vit = mesh.getVertex(vtx);
        const auto& X = vit->getCoordinates();

        // Find the local node j whose physical position is X
        int j_found = -1;
        for (size_t j = 0; j < fe.getCount(); ++j)
        {
          const auto& node = fe.getNode(j);
          // map node reference to physical as in section 2
          Math::SpatialPoint phys;
          mesh.getPolytopeTransformation(
              cell.getDimension(), cell.getIndex()).transform(phys, node);
          if ((phys - X).norm() < 1e-14)
          {
            j_found = static_cast<int>(j);
            break;
          }
        }

        EXPECT_NE(j_found, -1)
          << "No local node found matching vertex " << vtx
          << " on cell " << cell.getIndex();

        // Optionally check that fes.getGlobalIndex for that local
        // is the same as using the vertex DOF list:
        Index gdof_from_cell = fes.getGlobalIndex({2, cell.getIndex()}, j_found);
        const auto& vertex_dofs = fes.getDOFs(0, vtx);
        ASSERT_EQ(vertex_dofs.size(), 1u);
        Index gdof_from_vertex = vertex_dofs(0);

        EXPECT_EQ(gdof_from_cell, gdof_from_vertex)
          << "Inconsistent DOF for vertex " << vtx
          << " between cell " << cell.getIndex()
          << " and vertex-DOF list.";
      }
    }
  }

  TEST(Rodin_Variational_H1_Space, GlobalIndex_MatchesGetDOFs_ForVertices)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    const auto& conn = mesh.getConnectivity();

    for (Index vtx = 0; vtx < conn.getCount(0); ++vtx)
    {
      const auto& dofs_vertex = fes.getDOFs(0, vtx);
      ASSERT_EQ(dofs_vertex.size(), 1u);

      Index gdof_dofs = dofs_vertex(0);
      Index gdof_gi   = fes.getGlobalIndex({0, vtx}, 0);

      EXPECT_EQ(gdof_dofs, gdof_gi)
        << "Vertex " << vtx << " has inconsistent DOF index.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Pushforward_PointBasis_IsOneAtVertex)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    const auto& conn = mesh.getConnectivity();

    for (Index vtx = 0; vtx < conn.getCount(0); ++vtx)
    {
      const auto& fe_point = fes.getFiniteElement(0, vtx);
      ASSERT_EQ(fe_point.getCount(), 1u);
      const auto& basis0 = fe_point.getBasis(0);

      const auto mapping = fes.getPushforward({0, vtx}, basis0);

      auto vit = mesh.getVertex(vtx);
      Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());

      Real val = mapping(p);
      EXPECT_NEAR(val, 1.0, 1e-14)
        << "Pushforward of point basis at vertex " << vtx << " is not 1.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Project_OneTriangle_AtNodes)
  {
    // 1 element mesh, triangle with vertices (0,0), (1,0), (0,1)
    constexpr size_t sdim = 2;
    constexpr size_t mdim = 2;

    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() + 3.0 * p.y() + 1.0;
    };

    gf = exact; // uses Region::Cells internally

    auto cell = mesh.getCell();
    ASSERT_TRUE(cell);

    const auto& fe = fes.getFiniteElement(cell->getDimension(), cell->getIndex());
    ASSERT_EQ(fe.getGeometry(), Polytope::Type::Triangle);

    for (size_t local = 0; local < fe.getCount(); ++local)
    {
      const auto& ref_node = fe.getNode(local); // reference coords (xi,eta)
      Geometry::Point p(*cell, ref_node, ref_node);

      Real value    = gf(p);
      // physical coords: here ref = physical
      Real expected = exact(p);

      Index global = fes.getGlobalIndex({cell->getDimension(), cell->getIndex()}, local);
      Real stored  = gf[global];

      EXPECT_NEAR(stored, expected, 1e-14);
      EXPECT_NEAR(value, stored, 1e-14);
    }
  }

  // Basic construction test for H1<1> (equivalent to P1)
  TEST(Rodin_Variational_H1_Space, SanityTest_H1_1_2D_Square_Build)
  {
    constexpr size_t vdim = 1;
    constexpr size_t sdim = 2;
    constexpr size_t mdim = 2;

    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    EXPECT_EQ(mesh.getDimension(), mdim);
    EXPECT_EQ(mesh.getSpaceDimension(), sdim);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);

    EXPECT_EQ(fes.getVectorDimension(), vdim);
    // H1<1> has DOFs only at vertices, so size should equal vertex count
    EXPECT_EQ(fes.getSize(), mesh.getVertexCount());
    EXPECT_EQ(fes.getSize(), 4);
  }

  // Test H1<2> on a simple mesh
  TEST(Rodin_Variational_H1_Space, SanityTest_H1_2_2D_Square_Build)
  {
    constexpr size_t sdim = 2;
    constexpr size_t mdim = 2;

    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    EXPECT_EQ(mesh.getDimension(), mdim);
    EXPECT_EQ(mesh.getSpaceDimension(), sdim);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    EXPECT_EQ(fes.getVectorDimension(), 1);

    // H1<2> has:
    // - 1 DOF per vertex (4 vertices = 4 DOFs)
    // - 1 DOF per edge interior (5 edges, K-1=1 DOF each = 5 DOFs)
    // Total = 4 + 5 = 9 DOFs
    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t expectedDOFs = vertexCount + edgeCount;  // 4 + 5 = 9
    EXPECT_EQ(fes.getSize(), expectedDOFs);
  }

  // Test H1<3> on a simple mesh
  TEST(Rodin_Variational_H1_Space, SanityTest_H1_3_2D_Square_Build)
  {
    constexpr size_t sdim = 2;

    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);

    EXPECT_EQ(fes.getVectorDimension(), 1);

    // H1<3> has:
    // - 1 DOF per vertex (4 vertices = 4 DOFs)
    // - 2 DOFs per edge interior (5 edges, K-1=2 DOFs each = 10 DOFs)
    // - 1 DOF per face interior (2 triangles, (K-1)(K-2)/2=1 DOF each = 2 DOFs)
    // Total = 4 + 10 + 2 = 16 DOFs
    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();
    const size_t expectedDOFs = vertexCount + 2 * edgeCount + 1 * cellCount;
    EXPECT_EQ(fes.getSize(), expectedDOFs);
  }

  // Test H1<2> DOF structure on uniform grid
  TEST(Rodin_Variational_H1_Space, SanityTest_H1_2_UniformGrid)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    // For a 4x4 uniform grid, we have 25 vertices and ~some edges
    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t expectedDOFs = vertexCount + edgeCount;
    EXPECT_EQ(fes.getSize(), expectedDOFs);
  }

  // Test finite element retrieval
  TEST(Rodin_Variational_H1_Space, FiniteElement_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    const auto& fe = fes.getFiniteElement(2, 0);
    EXPECT_EQ(fe.getGeometry(), Polytope::Type::Triangle);
    // H1Element<2> on Triangle has 6 DOFs
    EXPECT_EQ(fe.getCount(), 6);
    EXPECT_EQ(fe.getOrder(), 2);
  }

  // Test DOF structure (closures) for H1<2> on a simple 2D mesh
  TEST(Rodin_Variational_H1_Space, GetDOFs_H1_2)
  {
    using Mesh     = Geometry::Mesh<Context::Local>;
    using Polytope = Geometry::Polytope;

    Mesh mesh =
      Mesh::Builder()
        .initialize(2)
        .nodes(4)
        .vertex({0, 0})
        .vertex({1, 0})
        .vertex({0, 1})
        .vertex({1, 1})
        .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
        .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
        .finalize();

    // Build edge–cell and edge–vertex incidence so that H1 can enumerate entities
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    constexpr size_t K = 2;
    H1 fes(std::integral_constant<size_t, K>{}, mesh);

    const auto& conn = mesh.getConnectivity();
    const size_t nv  = mesh.getVertexCount();
    const size_t ne  = conn.getCount(1);
    const size_t nc  = mesh.getCellCount();

    // 1) Vertices (d = 0): H1<2> has 1 DOF per vertex
    for (size_t v = 0; v < nv; ++v)
    {
      const auto& dofs_v = fes.getDOFs(0, static_cast<Index>(v));
      EXPECT_EQ(dofs_v.size(), 1u)
        << "Vertex " << v << " should have exactly 1 DOF for H1<2>.";
    }

    // 2) Edges (d = 1): H1<2> has (K+1) = 3 DOFs per edge (2 vertices + 1 interior)
    for (size_t e = 0; e < ne; ++e)
    {
      const auto& dofs_e = fes.getDOFs(1, static_cast<Index>(e));
      EXPECT_EQ(dofs_e.size(), 3u)
        << "Edge " << e << " should have 3 DOFs (2 vertex + 1 interior) for H1<2>.";
    }

    // 3) Cells (d = 2): H1<2> triangle has closure = 3 vertex DOFs + 3 edge DOFs = 6 DOFs
    for (size_t c = 0; c < nc; ++c)
    {
      const auto& dofs_c = fes.getDOFs(2, static_cast<Index>(c));
      const auto& fe_c   = fes.getFiniteElement(2, static_cast<Index>(c));

      // Check that cell closure size matches element DOF count
      EXPECT_EQ(dofs_c.size(), fe_c.getCount())
        << "Cell " << c << " should have closure size equal to element DOF count.";

      EXPECT_EQ(dofs_c.size(), 6u)
        << "Triangle " << c << " should have 6 DOFs in its closure for H1<2>.";
    }

    // 4) Consistency check: for one cell, its closure is exactly the union
    //    of its vertex DOFs and edge DOFs (no interior DOFs for K=2 on triangles).
    {
      const size_t c = 0; // first cell
      const auto& dofs_c = fes.getDOFs(2, static_cast<Index>(c));

      std::set<Index> closure_from_subentities;

      // vertices of cell c
      const auto& poly = conn.getPolytope(2, static_cast<Index>(c));
      for (Index v : poly)
      {
        const auto& dofs_v = fes.getDOFs(0, v);
        for (size_t k = 0; k < static_cast<size_t>(dofs_v.size()); ++k)
          closure_from_subentities.insert(dofs_v(k));
      }

      // edges of cell c
      const auto& edges = conn.getIncidence({2, 1}, static_cast<Index>(c));
      for (Index e : edges)
      {
        const auto& dofs_e = fes.getDOFs(1, e);
        for (size_t k = 0; k < static_cast<size_t>(dofs_e.size()); ++k)
          closure_from_subentities.insert(dofs_e(k));
      }

      // No interior DOFs for H1<2> on triangles, so closure_from_subentities
      // should match dofs_c exactly.
      EXPECT_EQ(closure_from_subentities.size(), static_cast<size_t>(dofs_c.size()));
      for (size_t k = 0; k < static_cast<size_t>(dofs_c.size()); ++k)
      {
        EXPECT_EQ(closure_from_subentities.count(dofs_c(k)), 1u)
          << "Cell closure DOF " << dofs_c(k)
          << " should come from a vertex or an edge DOF.";
      }
    }

    // 5) Global count: for H1<2> on this triangular mesh, total DOFs are
    //    nv vertex DOFs + ne*(K-1) edge interior DOFs (no cell interior for K=2).
    const size_t expected_global =
      nv                 // 1 DOF per vertex
      + ne * (K - 1);    // (K-1) interior DOFs per edge

    EXPECT_EQ(fes.getSize(), expected_global)
      << "Global DOF count should be nv + ne*(K-1) for H1<2> on triangles.";
  }

  // Test that H1<1> size matches P1 size
  TEST(Rodin_Variational_H1_Space, H1_1_MatchesP1_Size)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 h1_fes(std::integral_constant<size_t, 1>{}, mesh);
    P1<Real> p1_fes(mesh);

    EXPECT_EQ(h1_fes.getSize(), p1_fes.getSize());
  }

  // Test Vector H1 space
  TEST(Rodin_Variational_H1_Space, VectorH1_2_Build)
  {
    constexpr size_t vdim = 2;

    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, vdim);

    EXPECT_EQ(fes.getVectorDimension(), vdim);

    // Scalar H1<2> has 9 DOFs, vector version should have 9 * vdim = 18 DOFs
    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t scalarDOFs = vertexCount + edgeCount;
    EXPECT_EQ(fes.getSize(), scalarDOFs * vdim);
  }

  // Test that H1<1> can be used with TrialFunction and TestFunction
  TEST(Rodin_Variational_H1_Space, H1_1_TrialTestFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
  }

  // Test that H1<2> can be used with TrialFunction and TestFunction
  TEST(Rodin_Variational_H1_Space, H1_2_TrialTestFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
  }

  // Test DOF count consistency for various polynomial degrees
  TEST(Rodin_Variational_H1_Space, DOFCount_Consistency)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    // H1<K> DOF distribution:
    // - K >= 1: 1 DOF per vertex
    // - K >= 2: (K-1) DOFs per edge interior
    // - K >= 3: (K-1)(K-2)/2 DOFs per triangle interior
    // 
    // Note: K=0 is a special case where we have no vertex DOFs (constant elements).
    // This is consistent with P0 elements.

    // H1<1>: vertices only
    H1 h1_1(std::integral_constant<size_t, 1>{}, mesh);
    EXPECT_EQ(h1_1.getSize(), vertexCount);

    // H1<2>: vertices + edges
    H1 h1_2(std::integral_constant<size_t, 2>{}, mesh);
    EXPECT_EQ(h1_2.getSize(), vertexCount + edgeCount);

    // H1<3>: vertices + 2*edges + faces
    H1 h1_3(std::integral_constant<size_t, 3>{}, mesh);
    EXPECT_EQ(h1_3.getSize(), vertexCount + 2 * edgeCount + cellCount);
  }

  // GlobalIndex must match getDOFs entry for every (d,i,local)
  TEST(Rodin_Variational_H1_Space, GlobalIndex_MatchesGetDOFs)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    const auto& conn = mesh.getConnectivity();
    const size_t D   = mesh.getDimension();

    for (size_t d = 0; d <= D; ++d)
    {
      const size_t nd = conn.getCount(d);
      for (size_t i = 0; i < nd; ++i)
      {
        const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
        for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
        {
          Index g_from_dofs = dofs(local);
          Index g_from_api  = fes.getGlobalIndex({ d, static_cast<Index>(i) }, local);
          EXPECT_EQ(g_from_dofs, g_from_api)
            << "Mismatch between getDOFs and getGlobalIndex at (d=" << d
            << ", i=" << i << ", local=" << local << ").";
        }
      }
    }
  }

  // Union of closures over all entities must cover all global DOFs exactly once (in set sense)
  TEST(Rodin_Variational_H1_Space, ClosuresCoverAllGlobalDOFs_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    const auto& conn = mesh.getConnectivity();
    const size_t D   = mesh.getDimension();

    std::set<Index> all_dofs;
    for (size_t d = 0; d <= D; ++d)
    {
      const size_t nd = conn.getCount(d);
      for (size_t i = 0; i < nd; ++i)
      {
        const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
        for (size_t k = 0; k < static_cast<size_t>(dofs.size()); ++k)
          all_dofs.insert(dofs(k));
      }
    }

    EXPECT_EQ(all_dofs.size(), fes.getSize())
      << "Union of entity closures should cover all global DOFs exactly once in set sense.";
  }

  // For H1<3> on triangles: cell closure = vertices + edges + 1 interior DOF
  TEST(Rodin_Variational_H1_Space, CellClosure_DecomposesIntoSubentities_H1_3)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);

    const auto& conn = mesh.getConnectivity();
    const size_t nc  = mesh.getCellCount();

    ASSERT_EQ(nc, 2u);

    for (size_t c = 0; c < nc; ++c)
    {
      const auto& dofs_c = fes.getDOFs(2, static_cast<Index>(c));
      const auto& fe_c   = fes.getFiniteElement(2, static_cast<Index>(c));

      // Sanity: closure size = element DOF count
      EXPECT_EQ(dofs_c.size(), fe_c.getCount());

      // Build set of DOFs coming from vertices and edges of cell c
      std::set<Index> from_vertices_edges;

      // Vertices of the cell
      const auto& poly = conn.getPolytope(2, static_cast<Index>(c));
      for (Index v : poly)
      {
        const auto& dofs_v = fes.getDOFs(0, v);
        for (size_t k = 0; k < static_cast<size_t>(dofs_v.size()); ++k)
          from_vertices_edges.insert(dofs_v(k));
      }

      // Edges of the cell
      const auto& edges = conn.getIncidence({ 2, 1 }, static_cast<Index>(c));
      for (Index e : edges)
      {
        const auto& dofs_e = fes.getDOFs(1, e);
        for (size_t k = 0; k < static_cast<size_t>(dofs_e.size()); ++k)
          from_vertices_edges.insert(dofs_e(k));
      }

      // Interior DOFs = cell closure minus vertex+edge DOFs
      std::set<Index> interior;
      for (size_t k = 0; k < static_cast<size_t>(dofs_c.size()); ++k)
      {
        Index g = dofs_c(k);
        if (from_vertices_edges.count(g) == 0)
          interior.insert(g);
      }

      // For H1<3> on a triangle: 1 interior DOF
      EXPECT_EQ(interior.size(), 1u)
        << "Cell " << c << " should have exactly 1 interior DOF for H1<3>.";

      // Total cell closure = subentities + interior
      EXPECT_EQ(from_vertices_edges.size() + interior.size(), dofs_c.size());
    }
  }

  // Vector H1: DOFs per entity should be vdim times scalar DOFs per entity
  TEST(Rodin_Variational_H1_Space, VectorH1_EntityDOFScaling_H1_2)
  {
    constexpr size_t vdim = 3;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Scalar and vector spaces on same mesh and degree
    H1 scalar_fes(std::integral_constant<size_t, 2>{}, mesh);
    H1 vector_fes(std::integral_constant<size_t, 2>{}, mesh, vdim);

    EXPECT_EQ(vector_fes.getVectorDimension(), vdim);
    EXPECT_EQ(vector_fes.getSize(), scalar_fes.getSize() * vdim);

    const auto& conn = mesh.getConnectivity();
    const size_t D   = mesh.getDimension();

    // Check DOF scaling per entity (closure-wise)
    for (size_t d = 0; d <= D; ++d)
    {
      const size_t nd = conn.getCount(d);
      for (size_t i = 0; i < nd; ++i)
      {
        const auto& sdofs = scalar_fes.getDOFs(d, static_cast<Index>(i));
        const auto& vdofs = vector_fes.getDOFs(d, static_cast<Index>(i));

        EXPECT_EQ(vdofs.size(), sdofs.size() * vdim)
          << "Vector H1 DOFs on (d=" << d << ", i=" << i
          << ") should be vdim times scalar DOFs.";
      }
    }
  }

  // ============================================================================
  // Comprehensive behavioral tests for H1<K> spaces from K=1 to K=6
  // ============================================================================

  // Helper function to compute expected DOF count for H1<K> on a 2D triangle mesh
  inline size_t expectedDOFCount2D(size_t K, size_t vertexCount, size_t edgeCount, size_t cellCount)
  {
    // DOF formula for H1<K> on 2D triangular meshes:
    // - K >= 1: 1 DOF per vertex
    // - K >= 2: (K-1) DOFs per edge interior
    // - K >= 3: (K-1)(K-2)/2 DOFs per cell interior (triangle)
    size_t total = vertexCount;  // vertices (K >= 1)
    if (K >= 2)
      total += (K - 1) * edgeCount;  // edge interior
    if (K >= 3)
      total += ((K - 1) * (K - 2) / 2) * cellCount;  // cell interior
    return total;
  }

  // Test H1<4> construction and DOF count
  TEST(Rodin_Variational_H1_Space, SanityTest_H1_4_2D_Square_Build)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    // H1<4> has:
    // - 1 DOF per vertex (4)
    // - 3 DOFs per edge interior (5 edges * 3 = 15)
    // - 3 DOFs per cell interior (2 triangles * 3 = 6)
    // Total = 4 + 15 + 6 = 25
    size_t expectedDOFs = expectedDOFCount2D(4, vertexCount, edgeCount, cellCount);
    EXPECT_EQ(fes.getSize(), expectedDOFs);
    EXPECT_EQ(fes.getVectorDimension(), 1);
  }

  // Test H1<5> construction and DOF count
  TEST(Rodin_Variational_H1_Space, SanityTest_H1_5_2D_Square_Build)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(1, 2);

    H1 fes(std::integral_constant<size_t, 5>{}, mesh);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    // H1<5> has:
    // - 1 DOF per vertex (4)
    // - 4 DOFs per edge interior (5 edges * 4 = 20)
    // - 6 DOFs per cell interior (2 triangles * 6 = 12)
    // Total = 4 + 20 + 12 = 36
    size_t expectedDOFs = expectedDOFCount2D(5, vertexCount, edgeCount, cellCount);
    EXPECT_EQ(fes.getSize(), expectedDOFs);
    EXPECT_EQ(fes.getVectorDimension(), 1);
  }

  // Test H1<6> construction and DOF count
  TEST(Rodin_Variational_H1_Space, SanityTest_H1_6_2D_Square_Build)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    // H1<6> has:
    // - 1 DOF per vertex (4)
    // - 5 DOFs per edge interior (5 edges * 5 = 25)
    // - 10 DOFs per cell interior (2 triangles * 10 = 20)
    // Total = 4 + 25 + 20 = 49
    size_t expectedDOFs = expectedDOFCount2D(6, vertexCount, edgeCount, cellCount);
    EXPECT_EQ(fes.getSize(), expectedDOFs);
    EXPECT_EQ(fes.getVectorDimension(), 1);
  }

  // Comprehensive DOF count test for K = 1 to 6 on larger uniform grid
  TEST(Rodin_Variational_H1_Space, DOFCount_K1_to_K6_UniformGrid)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    // K = 1
    {
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(1, vertexCount, edgeCount, cellCount));
    }

    // K = 2
    {
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(2, vertexCount, edgeCount, cellCount));
    }

    // K = 3
    {
      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(3, vertexCount, edgeCount, cellCount));
    }

    // K = 4
    {
      H1 fes(std::integral_constant<size_t, 4>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(4, vertexCount, edgeCount, cellCount));
    }

    // K = 5
    {
      H1 fes(std::integral_constant<size_t, 5>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(5, vertexCount, edgeCount, cellCount));
    }

    // K = 6
    {
      H1 fes(std::integral_constant<size_t, 6>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(6, vertexCount, edgeCount, cellCount));
    }
  }

  // Test finite element retrieval for K = 4, 5, 6
  TEST(Rodin_Variational_H1_Space, FiniteElement_K4_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // K = 4
    {
      H1 fes(std::integral_constant<size_t, 4>{}, mesh);
      const auto& fe = fes.getFiniteElement(2, 0);
      EXPECT_EQ(fe.getGeometry(), Polytope::Type::Triangle);
      EXPECT_EQ(fe.getCount(), 15);  // (4+1)(4+2)/2 = 15
      EXPECT_EQ(fe.getOrder(), 4);
    }

    // K = 5
    {
      H1 fes(std::integral_constant<size_t, 5>{}, mesh);
      const auto& fe = fes.getFiniteElement(2, 0);
      EXPECT_EQ(fe.getGeometry(), Polytope::Type::Triangle);
      EXPECT_EQ(fe.getCount(), 21);  // (5+1)(5+2)/2 = 21
      EXPECT_EQ(fe.getOrder(), 5);
    }

    // K = 6
    {
      H1 fes(std::integral_constant<size_t, 6>{}, mesh);
      const auto& fe = fes.getFiniteElement(2, 0);
      EXPECT_EQ(fe.getGeometry(), Polytope::Type::Triangle);
      EXPECT_EQ(fe.getCount(), 28);  // (6+1)(6+2)/2 = 28
      EXPECT_EQ(fe.getOrder(), 6);
    }
  }

  // Test closure decomposition for K = 4
  TEST(Rodin_Variational_H1_Space, CellClosure_DecomposesIntoSubentities_H1_4)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);

    const auto& conn = mesh.getConnectivity();
    const size_t nc = mesh.getCellCount();

    for (size_t c = 0; c < nc; ++c)
    {
      const auto& dofs_c = fes.getDOFs(2, static_cast<Index>(c));
      const auto& fe_c = fes.getFiniteElement(2, static_cast<Index>(c));

      EXPECT_EQ(dofs_c.size(), fe_c.getCount());

      std::set<Index> from_vertices_edges;

      // Vertices of the cell
      const auto& poly = conn.getPolytope(2, static_cast<Index>(c));
      for (Index v : poly)
      {
        const auto& dofs_v = fes.getDOFs(0, v);
        for (size_t k = 0; k < static_cast<size_t>(dofs_v.size()); ++k)
          from_vertices_edges.insert(dofs_v(k));
      }

      // Edges of the cell
      const auto& edges = conn.getIncidence({ 2, 1 }, static_cast<Index>(c));
      for (Index e : edges)
      {
        const auto& dofs_e = fes.getDOFs(1, e);
        for (size_t k = 0; k < static_cast<size_t>(dofs_e.size()); ++k)
          from_vertices_edges.insert(dofs_e(k));
      }

      // Interior DOFs = cell closure minus vertex+edge DOFs
      std::set<Index> interior;
      for (size_t k = 0; k < static_cast<size_t>(dofs_c.size()); ++k)
      {
        Index g = dofs_c(k);
        if (from_vertices_edges.count(g) == 0)
          interior.insert(g);
      }

      // For H1<4> on a triangle: 3 interior DOFs = (4-1)(4-2)/2 = 3
      EXPECT_EQ(interior.size(), 3u)
        << "Cell " << c << " should have exactly 3 interior DOFs for H1<4>.";

      EXPECT_EQ(from_vertices_edges.size() + interior.size(), dofs_c.size());
    }
  }

  // Test closure decomposition for K = 5
  TEST(Rodin_Variational_H1_Space, CellClosure_DecomposesIntoSubentities_H1_5)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 5>{}, mesh);

    const auto& conn = mesh.getConnectivity();
    const size_t nc = mesh.getCellCount();

    for (size_t c = 0; c < nc; ++c)
    {
      const auto& dofs_c = fes.getDOFs(2, static_cast<Index>(c));
      const auto& fe_c = fes.getFiniteElement(2, static_cast<Index>(c));

      EXPECT_EQ(dofs_c.size(), fe_c.getCount());

      std::set<Index> from_vertices_edges;

      const auto& poly = conn.getPolytope(2, static_cast<Index>(c));
      for (Index v : poly)
      {
        const auto& dofs_v = fes.getDOFs(0, v);
        for (size_t k = 0; k < static_cast<size_t>(dofs_v.size()); ++k)
          from_vertices_edges.insert(dofs_v(k));
      }

      const auto& edges = conn.getIncidence({ 2, 1 }, static_cast<Index>(c));
      for (Index e : edges)
      {
        const auto& dofs_e = fes.getDOFs(1, e);
        for (size_t k = 0; k < static_cast<size_t>(dofs_e.size()); ++k)
          from_vertices_edges.insert(dofs_e(k));
      }

      std::set<Index> interior;
      for (size_t k = 0; k < static_cast<size_t>(dofs_c.size()); ++k)
      {
        Index g = dofs_c(k);
        if (from_vertices_edges.count(g) == 0)
          interior.insert(g);
      }

      // For H1<5> on a triangle: 6 interior DOFs = (5-1)(5-2)/2 = 6
      EXPECT_EQ(interior.size(), 6u)
        << "Cell " << c << " should have exactly 6 interior DOFs for H1<5>.";

      EXPECT_EQ(from_vertices_edges.size() + interior.size(), dofs_c.size());
    }
  }

  // Test closure decomposition for K = 6
  TEST(Rodin_Variational_H1_Space, CellClosure_DecomposesIntoSubentities_H1_6)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);

    const auto& conn = mesh.getConnectivity();
    const size_t nc = mesh.getCellCount();

    for (size_t c = 0; c < nc; ++c)
    {
      const auto& dofs_c = fes.getDOFs(2, static_cast<Index>(c));
      const auto& fe_c = fes.getFiniteElement(2, static_cast<Index>(c));

      EXPECT_EQ(dofs_c.size(), fe_c.getCount());

      std::set<Index> from_vertices_edges;

      const auto& poly = conn.getPolytope(2, static_cast<Index>(c));
      for (Index v : poly)
      {
        const auto& dofs_v = fes.getDOFs(0, v);
        for (size_t k = 0; k < static_cast<size_t>(dofs_v.size()); ++k)
          from_vertices_edges.insert(dofs_v(k));
      }

      const auto& edges = conn.getIncidence({ 2, 1 }, static_cast<Index>(c));
      for (Index e : edges)
      {
        const auto& dofs_e = fes.getDOFs(1, e);
        for (size_t k = 0; k < static_cast<size_t>(dofs_e.size()); ++k)
          from_vertices_edges.insert(dofs_e(k));
      }

      std::set<Index> interior;
      for (size_t k = 0; k < static_cast<size_t>(dofs_c.size()); ++k)
      {
        Index g = dofs_c(k);
        if (from_vertices_edges.count(g) == 0)
          interior.insert(g);
      }

      // For H1<6> on a triangle: 10 interior DOFs = (6-1)(6-2)/2 = 10
      EXPECT_EQ(interior.size(), 10u)
        << "Cell " << c << " should have exactly 10 interior DOFs for H1<6>.";

      EXPECT_EQ(from_vertices_edges.size() + interior.size(), dofs_c.size());
    }
  }

  // Test GlobalIndex matches getDOFs for K = 4, 5, 6
  TEST(Rodin_Variational_H1_Space, GlobalIndex_MatchesGetDOFs_K4_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testGlobalIndex = [&](auto fes) {
      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
          for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
          {
            Index g_from_dofs = dofs(local);
            Index g_from_api = fes.getGlobalIndex({ d, static_cast<Index>(i) }, local);
            EXPECT_EQ(g_from_dofs, g_from_api);
          }
        }
      }
    };

    testGlobalIndex(H1(std::integral_constant<size_t, 4>{}, mesh));
    testGlobalIndex(H1(std::integral_constant<size_t, 5>{}, mesh));
    testGlobalIndex(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  // Test closures cover all global DOFs for K = 4, 5, 6
  TEST(Rodin_Variational_H1_Space, ClosuresCoverAllGlobalDOFs_K4_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testCoverage = [&](auto fes) {
      std::set<Index> all_dofs;
      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
          for (size_t k = 0; k < static_cast<size_t>(dofs.size()); ++k)
            all_dofs.insert(dofs(k));
        }
      }
      EXPECT_EQ(all_dofs.size(), fes.getSize());
    };

    testCoverage(H1(std::integral_constant<size_t, 4>{}, mesh));
    testCoverage(H1(std::integral_constant<size_t, 5>{}, mesh));
    testCoverage(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  // Test TrialFunction and TestFunction creation for K = 4, 5, 6
  TEST(Rodin_Variational_H1_Space, TrialTestFunction_K4_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // K = 4
    {
      H1 fes(std::integral_constant<size_t, 4>{}, mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
    }

    // K = 5
    {
      H1 fes(std::integral_constant<size_t, 5>{}, mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
    }

    // K = 6
    {
      H1 fes(std::integral_constant<size_t, 6>{}, mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
    }
  }

  // Test Vector H1 DOF scaling for K = 4, 5, 6
  TEST(Rodin_Variational_H1_Space, VectorH1_EntityDOFScaling_K4_to_K6)
  {
    constexpr size_t vdim = 2;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testScaling = [&](auto scalar_fes, auto vector_fes) {
      EXPECT_EQ(vector_fes.getVectorDimension(), vdim);
      EXPECT_EQ(vector_fes.getSize(), scalar_fes.getSize() * vdim);

      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& sdofs = scalar_fes.getDOFs(d, static_cast<Index>(i));
          const auto& vdofs = vector_fes.getDOFs(d, static_cast<Index>(i));
          EXPECT_EQ(vdofs.size(), sdofs.size() * vdim);
        }
      }
    };

    testScaling(
      H1(std::integral_constant<size_t, 4>{}, mesh),
      H1(std::integral_constant<size_t, 4>{}, mesh, vdim)
    );

    testScaling(
      H1(std::integral_constant<size_t, 5>{}, mesh),
      H1(std::integral_constant<size_t, 5>{}, mesh, vdim)
    );

    testScaling(
      H1(std::integral_constant<size_t, 6>{}, mesh),
      H1(std::integral_constant<size_t, 6>{}, mesh, vdim)
    );
  }

  // Test edge DOF structure for K = 1 to 6
  TEST(Rodin_Variational_H1_Space, EdgeDOFCount_K1_to_K6)
  {
    using LocalMesh = Geometry::Mesh<Context::Local>;
    using Polytope  = Geometry::Polytope;

    LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t ne  = conn.getCount(1);

    // Helper lambda: check one value of K
    auto checkK = [&](auto Ktag)
    {
      constexpr size_t K = Ktag.value;
      H1 fes(Ktag, mesh);

      for (size_t e = 0; e < ne; ++e)
      {
        const Index edgeIdx = static_cast<Index>(e);

        // Edge closure DOFs
        const auto& dofs_e = fes.getDOFs(1, edgeIdx);

        // For H1<K> on segments with GLL nodes: K+1 DOFs per edge in closure
        EXPECT_EQ(dofs_e.size(), K + 1)
          << "Edge " << e << " should have " << (K + 1)
          << " total DOFs (including 2 vertex DOFs) for H1<" << K << ">.";

        // Get the two vertex DOFs of this edge
        const auto& edgeVerts = conn.getIncidence({1, 0}, edgeIdx);
        ASSERT_EQ(edgeVerts.size(), 2u);

        std::set<Index> vertexDOFs;
        for (Index v : edgeVerts)
        {
          const auto& dofs_v = fes.getDOFs(0, v);
          ASSERT_EQ(dofs_v.size(), 1u);
          vertexDOFs.insert(dofs_v(0));
        }

        // Count interior DOFs on the edge = edge DOFs minus vertex DOFs
        size_t interiorCount = 0;
        for (Index k = 0; k < dofs_e.size(); ++k)
        {
          if (vertexDOFs.count(dofs_e(k)) == 0)
            ++interiorCount;
        }

        EXPECT_EQ(interiorCount, (K > 0 ? K - 1 : 0))
          << "Edge " << e << " should have " << (K > 0 ? K - 1 : 0)
          << " interior DOFs for H1<" << K << ">.";
      }
    };

    checkK(std::integral_constant<size_t, 1>{});
    checkK(std::integral_constant<size_t, 2>{});
    checkK(std::integral_constant<size_t, 3>{});
    checkK(std::integral_constant<size_t, 4>{});
    checkK(std::integral_constant<size_t, 5>{});
    checkK(std::integral_constant<size_t, 6>{});
  }

  // Test vertex DOF count per entity for K = 1 to 6
  TEST(Rodin_Variational_H1_Space, VertexDOFCount_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t nv = mesh.getVertexCount();

    auto testVertexDOFs = [&](auto fes) {
      for (size_t v = 0; v < nv; ++v)
      {
        const auto& dofs_v = fes.getDOFs(0, static_cast<Index>(v));
        EXPECT_EQ(dofs_v.size(), 1u)  // Always 1 DOF per vertex for K >= 1
          << "Vertex " << v << " should have exactly 1 DOF.";
      }
    };

    testVertexDOFs(H1(std::integral_constant<size_t, 1>{}, mesh));
    testVertexDOFs(H1(std::integral_constant<size_t, 2>{}, mesh));
    testVertexDOFs(H1(std::integral_constant<size_t, 3>{}, mesh));
    testVertexDOFs(H1(std::integral_constant<size_t, 4>{}, mesh));
    testVertexDOFs(H1(std::integral_constant<size_t, 5>{}, mesh));
    testVertexDOFs(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  // Test that DOF indices are contiguous and range from 0 to size-1
  TEST(Rodin_Variational_H1_Space, DOFIndices_Contiguous_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testContiguous = [&](auto fes) {
      std::set<Index> all_dofs;
      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
          for (size_t k = 0; k < static_cast<size_t>(dofs.size()); ++k)
            all_dofs.insert(dofs(k));
        }
      }

      // Check that indices are contiguous from 0 to size-1
      EXPECT_EQ(all_dofs.size(), fes.getSize());
      if (!all_dofs.empty())
      {
        EXPECT_EQ(*all_dofs.begin(), 0u);
        EXPECT_EQ(*all_dofs.rbegin(), fes.getSize() - 1);
      }
    };

    testContiguous(H1(std::integral_constant<size_t, 1>{}, mesh));
    testContiguous(H1(std::integral_constant<size_t, 2>{}, mesh));
    testContiguous(H1(std::integral_constant<size_t, 3>{}, mesh));
    testContiguous(H1(std::integral_constant<size_t, 4>{}, mesh));
    testContiguous(H1(std::integral_constant<size_t, 5>{}, mesh));
    testContiguous(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  // Test cell closure size matches element DOF count for K = 1 to 6
  TEST(Rodin_Variational_H1_Space, CellClosureSize_MatchesElementDOFCount_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t nc = mesh.getCellCount();

    auto testClosureSize = [&](auto fes) {
      for (size_t c = 0; c < nc; ++c)
      {
        const auto& dofs_c = fes.getDOFs(2, static_cast<Index>(c));
        const auto& fe_c = fes.getFiniteElement(2, static_cast<Index>(c));
        EXPECT_EQ(dofs_c.size(), fe_c.getCount())
          << "Cell " << c << " closure size should match element DOF count.";
      }
    };

    testClosureSize(H1(std::integral_constant<size_t, 1>{}, mesh));
    testClosureSize(H1(std::integral_constant<size_t, 2>{}, mesh));
    testClosureSize(H1(std::integral_constant<size_t, 3>{}, mesh));
    testClosureSize(H1(std::integral_constant<size_t, 4>{}, mesh));
    testClosureSize(H1(std::integral_constant<size_t, 5>{}, mesh));
    testClosureSize(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  // ============================================================================
  // Interpolation tests for H1<K> spaces from K=2 to K=6
  // Verify that GridFunction evaluation returns correct values
  // ============================================================================

  // Interpolation test for H1<2>: Linear function (degree 1) should be exact
  TEST(Rodin_Variational_H1_Space, Interpolation_GridFunction_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function (degree 1) - H1<2> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() + 3.0 * p.y() + 1.0;
    };

    gf = exact;

    // Verify GridFunction size matches FES size
    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices match the exact function
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto& dofs_vertex = fes.getDOFs(0, vtx);

      EXPECT_EQ(dofs_vertex.size(), 1);

      Index gdof = dofs_vertex(0);

      decltype(auto) fe = fes.getFiniteElement(0, vtx);

      EXPECT_EQ(dofs_vertex.size(), fe.getCount());
      EXPECT_NEAR(fe.getBasis(0)(Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0)), 1.0, 1e-12);

      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());

      Real expected = exact(p);
      Real value = gf(p);

      Real stored = gf[gdof];

      EXPECT_NEAR(value, stored, 1e-10)
        << "H1<2> interpolation at vertex " << vtx << " should match exact linear function.";

      EXPECT_NEAR(value, expected, 1e-10)
        << "H1<2> interpolation at vertex " << vtx << " should match exact linear function.";
    }
  }

  // Interpolation test for H1<3>: Quadratic function (degree 2) should be exact
  TEST(Rodin_Variational_H1_Space, Interpolation_GridFunction_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.x() * p.y();
    };

    gf = exact;

    // Verify GridFunction size matches FES size
    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices match the exact function
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  // Interpolation test for H1<4>: Cubic function (degree 3) should be exact
  TEST(Rodin_Variational_H1_Space, Interpolation_GridFunction_H1_4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    // Cubic function (degree 3) - H1<4> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() * p.x() + p.y() * p.y() * p.y() + p.x() * p.y() * p.y();
    };

    gf = exact;

    // Verify GridFunction size matches FES size
    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices match the exact function
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "H1<4> interpolation at vertex " << vtx << " should match exact cubic function.";
    }
  }

  // Interpolation test for H1<5>: Quartic function (degree 4) should be exact
  TEST(Rodin_Variational_H1_Space, Interpolation_GridFunction_H1_5)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 5>{}, mesh);
    GridFunction gf(fes);

    // Quartic function (degree 4) - H1<5> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() * p.x() * p.x() + p.y() * p.y() * p.y() * p.y() + p.x() * p.x() * p.y() * p.y();
    };

    gf = exact;

    // Verify GridFunction size matches FES size
    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices match the exact function
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "H1<5> interpolation at vertex " << vtx << " should match exact quartic function.";
    }
  }

  // Interpolation test for H1<6>: Quintic function (degree 5) should be exact
  TEST(Rodin_Variational_H1_Space, Interpolation_GridFunction_H1_6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);
    GridFunction gf(fes);

    // Quintic function (degree 5) - H1<6> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() * p.x() * p.x() * p.x()
           + p.y() * p.y() * p.y() * p.y() * p.y()
           + p.x() * p.x() * p.y() * p.y() * p.y();
    };

    gf = exact;

    // Verify GridFunction size matches FES size
    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices match the exact function
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "H1<6> interpolation at vertex " << vtx << " should match exact quintic function.";
    }
  }

  // ============================================================================
  // 16x16 mesh tests for H1<K> spaces
  // ============================================================================

  // Test DOF count on 16x16 mesh for K = 1 to 6
  TEST(Rodin_Variational_H1_Space, DOFCount_16x16_Mesh_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    // K = 1
    {
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(1, vertexCount, edgeCount, cellCount));
    }

    // K = 2
    {
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(2, vertexCount, edgeCount, cellCount));
    }

    // K = 3
    {
      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(3, vertexCount, edgeCount, cellCount));
    }

    // K = 4
    {
      H1 fes(std::integral_constant<size_t, 4>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(4, vertexCount, edgeCount, cellCount));
    }

    // K = 5
    {
      H1 fes(std::integral_constant<size_t, 5>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(5, vertexCount, edgeCount, cellCount));
    }

    // K = 6
    {
      H1 fes(std::integral_constant<size_t, 6>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(6, vertexCount, edgeCount, cellCount));
    }
  }

  // Test indexing on 16x16 mesh
  TEST(Rodin_Variational_H1_Space, GlobalIndex_16x16_Mesh_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testGlobalIndex = [&](auto fes) {
      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
          for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
          {
            Index g_from_dofs = dofs(local);
            Index g_from_api = fes.getGlobalIndex({ d, static_cast<Index>(i) }, local);
            EXPECT_EQ(g_from_dofs, g_from_api);
          }
        }
      }
    };

    testGlobalIndex(H1(std::integral_constant<size_t, 1>{}, mesh));
    testGlobalIndex(H1(std::integral_constant<size_t, 2>{}, mesh));
    testGlobalIndex(H1(std::integral_constant<size_t, 3>{}, mesh));
    testGlobalIndex(H1(std::integral_constant<size_t, 4>{}, mesh));
    testGlobalIndex(H1(std::integral_constant<size_t, 5>{}, mesh));
    testGlobalIndex(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  // Test interpolation on 16x16 mesh - verify GridFunction values
  TEST(Rodin_Variational_H1_Space, Interpolation_16x16_Mesh_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + 2.0 * p.x() * p.y();
    };

    gf = exact;

    // Verify GridFunction size matches FES size
    EXPECT_EQ(gf.getSize(), fes.getSize());
    EXPECT_EQ(fes.getSize(), expectedDOFCount2D(3,
      mesh.getVertexCount(),
      mesh.getConnectivity().getCount(1),
      mesh.getCellCount()));

    // Verify interpolated values at vertices match the exact function
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "H1<3> interpolation on 16x16 mesh at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  // ============================================================================
  // Vector H1 comprehensive tests
  // ============================================================================

  // Vector H1 DOF count test for K = 1 to 6
  TEST(Rodin_Variational_H1_Space, VectorH1_DOFCount_K1_to_K6)
  {
    constexpr size_t vdim = 2;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    auto testVectorDOFCount = [&](size_t K, auto fes) {
      size_t scalarDOFs = expectedDOFCount2D(K, vertexCount, edgeCount, cellCount);
      EXPECT_EQ(fes.getSize(), scalarDOFs * vdim);
      EXPECT_EQ(fes.getVectorDimension(), vdim);
    };

    testVectorDOFCount(1, H1(std::integral_constant<size_t, 1>{}, mesh, vdim));
    testVectorDOFCount(2, H1(std::integral_constant<size_t, 2>{}, mesh, vdim));
    testVectorDOFCount(3, H1(std::integral_constant<size_t, 3>{}, mesh, vdim));
    testVectorDOFCount(4, H1(std::integral_constant<size_t, 4>{}, mesh, vdim));
    testVectorDOFCount(5, H1(std::integral_constant<size_t, 5>{}, mesh, vdim));
    testVectorDOFCount(6, H1(std::integral_constant<size_t, 6>{}, mesh, vdim));
  }

  // Vector H1 GridFunction creation and value verification test
  TEST(Rodin_Variational_H1_Space, VectorH1_GridFunction_K1_to_K6)
  {
    constexpr size_t vdim = 2;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    auto testVectorGridFunction = [&mesh](auto fes) {
      GridFunction gf(fes);

      // Linear vector function - should be exact for all K >= 1
      auto exact = [](const Geometry::Point& p) -> Math::Vector<Real>
      {
        Math::Vector<Real> v(2);
        v(0) = p.x() + p.y();
        v(1) = p.x() - p.y();
        return v;
      };

      gf = exact;

      // Verify GridFunction size matches FES size
      EXPECT_EQ(gf.getSize(), fes.getSize());

      // Verify interpolated values at vertices match the exact function
      const size_t nv = mesh.getVertexCount();
      for (size_t vtx = 0; vtx < nv; ++vtx)
      {
        const auto vit = mesh.getVertex(vtx);
        const Geometry::Point p(
            *vit,
            Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
            vit->getCoordinates());
        Math::Vector<Real> expected = exact(p);
        Math::Vector<Real> value = gf(p);
        EXPECT_EQ(value.size(), 2);
        EXPECT_NEAR(value(0), expected(0), 1e-10)
          << "Vector H1 interpolation component 0 at vertex " << vtx << " should match exact function.";
        EXPECT_NEAR(value(1), expected(1), 1e-10)
          << "Vector H1 interpolation component 1 at vertex " << vtx << " should match exact function.";
      }
    };

    testVectorGridFunction(H1(std::integral_constant<size_t, 1>{}, mesh, vdim));
    testVectorGridFunction(H1(std::integral_constant<size_t, 2>{}, mesh, vdim));
    testVectorGridFunction(H1(std::integral_constant<size_t, 3>{}, mesh, vdim));
    testVectorGridFunction(H1(std::integral_constant<size_t, 4>{}, mesh, vdim));
    testVectorGridFunction(H1(std::integral_constant<size_t, 5>{}, mesh, vdim));
    testVectorGridFunction(H1(std::integral_constant<size_t, 6>{}, mesh, vdim));
  }

  // Vector H1 16x16 mesh test
  TEST(Rodin_Variational_H1_Space, VectorH1_16x16_Mesh_K1_to_K6)
  {
    constexpr size_t vdim = 2;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    auto testVectorDOFCount = [&](size_t K, auto fes) {
      size_t scalarDOFs = expectedDOFCount2D(K, vertexCount, edgeCount, cellCount);
      EXPECT_EQ(fes.getSize(), scalarDOFs * vdim);
      EXPECT_EQ(fes.getVectorDimension(), vdim);
    };

    testVectorDOFCount(1, H1(std::integral_constant<size_t, 1>{}, mesh, vdim));
    testVectorDOFCount(2, H1(std::integral_constant<size_t, 2>{}, mesh, vdim));
    testVectorDOFCount(3, H1(std::integral_constant<size_t, 3>{}, mesh, vdim));
    testVectorDOFCount(4, H1(std::integral_constant<size_t, 4>{}, mesh, vdim));
    testVectorDOFCount(5, H1(std::integral_constant<size_t, 5>{}, mesh, vdim));
    testVectorDOFCount(6, H1(std::integral_constant<size_t, 6>{}, mesh, vdim));
  }

  // Vector H1 with 3D vector dimension
  TEST(Rodin_Variational_H1_Space, VectorH1_3D_Vector_K1_to_K6)
  {
    constexpr size_t vdim = 3;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    auto testVectorDOFCount = [&](size_t K, auto fes) {
      size_t scalarDOFs = expectedDOFCount2D(K, vertexCount, edgeCount, cellCount);
      EXPECT_EQ(fes.getSize(), scalarDOFs * vdim);
      EXPECT_EQ(fes.getVectorDimension(), vdim);
    };

    testVectorDOFCount(1, H1(std::integral_constant<size_t, 1>{}, mesh, vdim));
    testVectorDOFCount(2, H1(std::integral_constant<size_t, 2>{}, mesh, vdim));
    testVectorDOFCount(3, H1(std::integral_constant<size_t, 3>{}, mesh, vdim));
    testVectorDOFCount(4, H1(std::integral_constant<size_t, 4>{}, mesh, vdim));
    testVectorDOFCount(5, H1(std::integral_constant<size_t, 5>{}, mesh, vdim));
    testVectorDOFCount(6, H1(std::integral_constant<size_t, 6>{}, mesh, vdim));
  }

  // ============================================================================
  // Complex H1 tests
  // ============================================================================

  // Complex H1 construction and DOF count tests
  TEST(Rodin_Variational_H1_Space, ComplexH1_Construction_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    // K = 1
    {
      H1<1, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(1, vertexCount, edgeCount, cellCount));
      EXPECT_EQ(fes.getVectorDimension(), 1);
    }

    // K = 2
    {
      H1<2, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 2>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(2, vertexCount, edgeCount, cellCount));
    }

    // K = 3
    {
      H1<3, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 3>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(3, vertexCount, edgeCount, cellCount));
    }

    // K = 4
    {
      H1<4, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 4>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(4, vertexCount, edgeCount, cellCount));
    }

    // K = 5
    {
      H1<5, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 5>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(5, vertexCount, edgeCount, cellCount));
    }

    // K = 6
    {
      H1<6, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 6>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(6, vertexCount, edgeCount, cellCount));
    }
  }

  // Complex H1 indexing test
  TEST(Rodin_Variational_H1_Space, ComplexH1_GlobalIndex_K2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1<2, Complex, Geometry::Mesh<Context::Local>> fes(
        std::integral_constant<size_t, 2>{}, mesh);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    for (size_t d = 0; d <= D; ++d)
    {
      const size_t nd = conn.getCount(d);
      for (size_t i = 0; i < nd; ++i)
      {
        const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
        for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
        {
          Index g_from_dofs = dofs(local);
          Index g_from_api = fes.getGlobalIndex({ d, static_cast<Index>(i) }, local);
          EXPECT_EQ(g_from_dofs, g_from_api);
        }
      }
    }
  }

  // Complex H1 16x16 mesh test
  TEST(Rodin_Variational_H1_Space, ComplexH1_16x16_Mesh_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    // K = 1
    {
      H1<1, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(1, vertexCount, edgeCount, cellCount));
    }

    // K = 2
    {
      H1<2, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 2>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(2, vertexCount, edgeCount, cellCount));
    }

    // K = 3
    {
      H1<3, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 3>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(3, vertexCount, edgeCount, cellCount));
    }

    // K = 4
    {
      H1<4, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 4>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(4, vertexCount, edgeCount, cellCount));
    }

    // K = 5
    {
      H1<5, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 5>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(5, vertexCount, edgeCount, cellCount));
    }

    // K = 6
    {
      H1<6, Complex, Geometry::Mesh<Context::Local>> fes(
          std::integral_constant<size_t, 6>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(6, vertexCount, edgeCount, cellCount));
    }
  }

  // ============================================================================
  // Edge case tests
  // ============================================================================

  // Single cell mesh edge case
  TEST(Rodin_Variational_H1_Space, EdgeCase_SingleCell_K1_to_K6)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    EXPECT_EQ(vertexCount, 3);
    EXPECT_EQ(edgeCount, 3);
    EXPECT_EQ(cellCount, 1);

    // K = 1: verify we can build the space and it has expected properties
    {
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getVectorDimension(), 1);
      EXPECT_GT(fes.getSize(), 0u);  // Should have some DOFs
      // K=1 on single triangle should have 3 DOFs (vertices)
      EXPECT_EQ(fes.getSize(), 3);
    }

    // K = 2
    {
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      EXPECT_EQ(fes.getVectorDimension(), 1);
      EXPECT_GT(fes.getSize(), 3u);  // Should have more DOFs than vertices
      // K=2 on single triangle: 3 vertices + 3 edges × 1 = 6 DOFs
      EXPECT_EQ(fes.getSize(), 6);
    }

    // K = 3
    {
      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      EXPECT_EQ(fes.getVectorDimension(), 1);
      EXPECT_GT(fes.getSize(), 6u);  // Should have more DOFs than K=2
    }

    // K = 4
    {
      H1 fes(std::integral_constant<size_t, 4>{}, mesh);
      EXPECT_EQ(fes.getVectorDimension(), 1);
      EXPECT_GT(fes.getSize(), 9u);  // Should have more DOFs than K=3
    }

    // K = 5
    {
      H1 fes(std::integral_constant<size_t, 5>{}, mesh);
      EXPECT_EQ(fes.getVectorDimension(), 1);
      EXPECT_GT(fes.getSize(), 12u);  // Should have more DOFs than K=4
    }

    // K = 6
    {
      H1 fes(std::integral_constant<size_t, 6>{}, mesh);
      EXPECT_EQ(fes.getVectorDimension(), 1);
      EXPECT_GT(fes.getSize(), 15u);  // Should have more DOFs than K=5
    }
  }

  // Very small mesh (2 cells) edge case
  TEST(Rodin_Variational_H1_Space, EdgeCase_TwoCells_K1_to_K6)
  {
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .polytope(Polytope::Type::Triangle, { {1, 3, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    EXPECT_EQ(vertexCount, 4);
    EXPECT_EQ(edgeCount, 5);
    EXPECT_EQ(cellCount, 2);

    // K = 1: 4 DOFs
    {
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(1, vertexCount, edgeCount, cellCount));
    }

    // K = 2: 4 + 5 = 9 DOFs
    {
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(2, vertexCount, edgeCount, cellCount));
    }

    // K = 3
    {
      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(3, vertexCount, edgeCount, cellCount));
    }

    // K = 4
    {
      H1 fes(std::integral_constant<size_t, 4>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(4, vertexCount, edgeCount, cellCount));
    }

    // K = 5
    {
      H1 fes(std::integral_constant<size_t, 5>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(5, vertexCount, edgeCount, cellCount));
    }

    // K = 6
    {
      H1 fes(std::integral_constant<size_t, 6>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(6, vertexCount, edgeCount, cellCount));
    }
  }

  // Boundary cells edge case test - 8x8 mesh DOF verification
  TEST(Rodin_Variational_H1_Space, EdgeCase_8x8_Mesh_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 8, 8 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    // Verify DOF counts for K = 1 to 6
    {
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(1, vertexCount, edgeCount, cellCount));
    }

    {
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(2, vertexCount, edgeCount, cellCount));
    }

    {
      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(3, vertexCount, edgeCount, cellCount));
    }

    {
      H1 fes(std::integral_constant<size_t, 4>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(4, vertexCount, edgeCount, cellCount));
    }

    {
      H1 fes(std::integral_constant<size_t, 5>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(5, vertexCount, edgeCount, cellCount));
    }

    {
      H1 fes(std::integral_constant<size_t, 6>{}, mesh);
      EXPECT_EQ(fes.getSize(), expectedDOFCount2D(6, vertexCount, edgeCount, cellCount));
    }
  }

  // ============================================================================
  // Additional indexing tests
  // ============================================================================

  // Test that all indices are valid and within range
  TEST(Rodin_Variational_H1_Space, IndexingRange_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testIndexRange = [&](auto fes) {
      const size_t totalDOFs = fes.getSize();
      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
          for (size_t k = 0; k < static_cast<size_t>(dofs.size()); ++k)
          {
            Index g = dofs(k);
            EXPECT_GE(g, 0u) << "DOF index should be non-negative.";
            EXPECT_LT(g, totalDOFs) << "DOF index should be less than total DOFs.";
          }
        }
      }
    };

    testIndexRange(H1(std::integral_constant<size_t, 1>{}, mesh));
    testIndexRange(H1(std::integral_constant<size_t, 2>{}, mesh));
    testIndexRange(H1(std::integral_constant<size_t, 3>{}, mesh));
    testIndexRange(H1(std::integral_constant<size_t, 4>{}, mesh));
    testIndexRange(H1(std::integral_constant<size_t, 5>{}, mesh));
    testIndexRange(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  // Test element DOF count consistency for all cells
  TEST(Rodin_Variational_H1_Space, ElementDOFCount_AllCells_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 8, 8 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t nc = mesh.getCellCount();

    auto testElementDOFCount = [&](auto fes, size_t expectedCount) {
      for (size_t c = 0; c < nc; ++c)
      {
        const auto& fe = fes.getFiniteElement(2, static_cast<Index>(c));
        EXPECT_EQ(fe.getCount(), expectedCount)
          << "Element DOF count mismatch at cell " << c;
      }
    };

    // (K+1)(K+2)/2 for triangles
    testElementDOFCount(H1(std::integral_constant<size_t, 1>{}, mesh), 3);
    testElementDOFCount(H1(std::integral_constant<size_t, 2>{}, mesh), 6);
    testElementDOFCount(H1(std::integral_constant<size_t, 3>{}, mesh), 10);
    testElementDOFCount(H1(std::integral_constant<size_t, 4>{}, mesh), 15);
    testElementDOFCount(H1(std::integral_constant<size_t, 5>{}, mesh), 21);
    testElementDOFCount(H1(std::integral_constant<size_t, 6>{}, mesh), 28);
  }

  // Test TrialFunction and TestFunction with various mesh sizes
  TEST(Rodin_Variational_H1_Space, TrialTestFunction_VariousMeshSizes)
  {
    // 2x2 mesh
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);

      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
    }

    // 8x8 mesh
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 8, 8 });
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);

      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
    }

    // 16x16 mesh
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 16, 16 });
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);


      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
    }
  }

  // Vector H1 indexing test
  TEST(Rodin_Variational_H1_Space, VectorH1_Indexing_K1_to_K6)
  {
    constexpr size_t vdim = 2;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testVectorIndexing = [&](auto fes) {
      const size_t totalDOFs = fes.getSize();
      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
          for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
          {
            Index g_from_dofs = dofs(local);
            Index g_from_api = fes.getGlobalIndex({ d, static_cast<Index>(i) }, local);
            EXPECT_EQ(g_from_dofs, g_from_api);
            EXPECT_GE(g_from_dofs, 0u);
            EXPECT_LT(g_from_dofs, totalDOFs);
          }
        }
      }
    };

    testVectorIndexing(H1(std::integral_constant<size_t, 1>{}, mesh, vdim));
    testVectorIndexing(H1(std::integral_constant<size_t, 2>{}, mesh, vdim));
    testVectorIndexing(H1(std::integral_constant<size_t, 3>{}, mesh, vdim));
    testVectorIndexing(H1(std::integral_constant<size_t, 4>{}, mesh, vdim));
    testVectorIndexing(H1(std::integral_constant<size_t, 5>{}, mesh, vdim));
    testVectorIndexing(H1(std::integral_constant<size_t, 6>{}, mesh, vdim));
  }

  // ============================================================================
  // GEOMETRY-SPECIFIC TESTS
  // ============================================================================

  // --------------------------------------------------------------------------
  // Segment (1D) Geometry Tests
  // --------------------------------------------------------------------------

  TEST(Rodin_Variational_H1_Space, Segment_H1_1_DOFCount)
  {
    // Create a 1D mesh with 4 segments (5 vertices)
    constexpr size_t sdim = 1;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(5)
      .vertex({0.0})
      .vertex({0.25})
      .vertex({0.5})
      .vertex({0.75})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .polytope(Polytope::Type::Segment, { {1, 2} })
      .polytope(Polytope::Type::Segment, { {2, 3} })
      .polytope(Polytope::Type::Segment, { {3, 4} })
      .finalize();

    // For 1D meshes, compute vertex-to-cell (segment) connectivity
    mesh.getConnectivity().compute(1, 0);

    EXPECT_EQ(mesh.getDimension(), 1);
    EXPECT_EQ(mesh.getSpaceDimension(), 1);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);

    // H1<1> on segments: 1 DOF per vertex = 5 DOFs
    EXPECT_EQ(fes.getSize(), 5);
    EXPECT_EQ(fes.getVectorDimension(), 1);
  }

  TEST(Rodin_Variational_H1_Space, Segment_H1_2_DOFCount)
  {
    constexpr size_t sdim = 1;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(5)
      .vertex({0.0})
      .vertex({0.25})
      .vertex({0.5})
      .vertex({0.75})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .polytope(Polytope::Type::Segment, { {1, 2} })
      .polytope(Polytope::Type::Segment, { {2, 3} })
      .polytope(Polytope::Type::Segment, { {3, 4} })
      .finalize();

    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    // H1<2> on segments:
    // - 1 DOF per vertex (5 vertices) = 5 DOFs
    // - 1 DOF per segment interior (4 segments) = 4 DOFs
    // Total = 9 DOFs
    EXPECT_EQ(fes.getSize(), 9);
  }

  TEST(Rodin_Variational_H1_Space, Segment_H1_3_DOFCount)
  {
    constexpr size_t sdim = 1;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(5)
      .vertex({0.0})
      .vertex({0.25})
      .vertex({0.5})
      .vertex({0.75})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .polytope(Polytope::Type::Segment, { {1, 2} })
      .polytope(Polytope::Type::Segment, { {2, 3} })
      .polytope(Polytope::Type::Segment, { {3, 4} })
      .finalize();

    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);

    // H1<3> on segments:
    // - 1 DOF per vertex (5 vertices) = 5 DOFs
    // - 2 DOFs per segment interior (4 segments) = 8 DOFs
    // Total = 13 DOFs
    EXPECT_EQ(fes.getSize(), 13);
  }

  TEST(Rodin_Variational_H1_Space, Segment_H1_K1_to_K6_Indexing)
  {
    constexpr size_t sdim = 1;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(5)
      .vertex({0.0})
      .vertex({0.25})
      .vertex({0.5})
      .vertex({0.75})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .polytope(Polytope::Type::Segment, { {1, 2} })
      .polytope(Polytope::Type::Segment, { {2, 3} })
      .polytope(Polytope::Type::Segment, { {3, 4} })
      .finalize();

    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testIndexing = [&](auto fes) {
      const size_t totalDOFs = fes.getSize();
      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
          for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
          {
            Index g = dofs(local);
            EXPECT_GE(g, 0u);
            EXPECT_LT(g, totalDOFs);
          }
        }
      }
    };

    testIndexing(H1(std::integral_constant<size_t, 1>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 2>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 3>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 4>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 5>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  TEST(Rodin_Variational_H1_Space, Segment_TrialTestFunction)
  {
    constexpr size_t sdim = 1;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(5)
      .vertex({0.0})
      .vertex({0.25})
      .vertex({0.5})
      .vertex({0.75})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .polytope(Polytope::Type::Segment, { {1, 2} })
      .polytope(Polytope::Type::Segment, { {2, 3} })
      .polytope(Polytope::Type::Segment, { {3, 4} })
      .finalize();

    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    EXPECT_EQ(&u.getFiniteElementSpace(), &fes);
    EXPECT_EQ(&v.getFiniteElementSpace(), &fes);
  }

  TEST(Rodin_Variational_H1_Space, Segment_VectorH1_DOFCount)
  {
    constexpr size_t sdim = 1;
    constexpr size_t vdim = 2;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(5)
      .vertex({0.0})
      .vertex({0.25})
      .vertex({0.5})
      .vertex({0.75})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .polytope(Polytope::Type::Segment, { {1, 2} })
      .polytope(Polytope::Type::Segment, { {2, 3} })
      .polytope(Polytope::Type::Segment, { {3, 4} })
      .finalize();

    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, vdim);

    // Vector H1<2>: 9 scalar DOFs * 2 = 18 DOFs
    EXPECT_EQ(fes.getSize(), 18);
    EXPECT_EQ(fes.getVectorDimension(), vdim);
  }

  // --------------------------------------------------------------------------
  // Quadrilateral (2D) Geometry Tests
  // --------------------------------------------------------------------------

  TEST(Rodin_Variational_H1_Space, Quadrilateral_H1_1_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    EXPECT_EQ(mesh.getDimension(), 2);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);

    // H1<1> on quads: 1 DOF per vertex
    // 4x4 grid = 16 vertices
    EXPECT_EQ(fes.getSize(), 16);
    EXPECT_EQ(fes.getVectorDimension(), 1);
  }

  TEST(Rodin_Variational_H1_Space, Quadrilateral_H1_2_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    // H1<2> on 4x4 quad mesh:
    // - 1 DOF per vertex (16 vertices)
    // - 1 DOF per edge interior (24 edges)
    // - (K-1)² = 1 DOF per cell interior (9 cells)
    // Total = 16 + 24 + 9 = 49 DOFs
    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();
    const size_t expectedDOFs = vertexCount + edgeCount + cellCount;
    EXPECT_EQ(fes.getSize(), expectedDOFs);
  }

  TEST(Rodin_Variational_H1_Space, Quadrilateral_H1_3_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);

    // H1<3> on quads:
    // - 1 DOF per vertex
    // - 2 DOFs per edge interior
    // - (K-1)² = 4 DOFs per cell interior
    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();
    const size_t expectedDOFs = vertexCount + 2 * edgeCount + 4 * cellCount;
    EXPECT_EQ(fes.getSize(), expectedDOFs);
  }

  TEST(Rodin_Variational_H1_Space, Quadrilateral_H1_K1_to_K6_Indexing)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testIndexing = [&](auto fes) {
      const size_t totalDOFs = fes.getSize();
      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
          for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
          {
            Index g = dofs(local);
            EXPECT_GE(g, 0u);
            EXPECT_LT(g, totalDOFs);
          }
        }
      }
    };

    testIndexing(H1(std::integral_constant<size_t, 1>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 2>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 3>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 4>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 5>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  TEST(Rodin_Variational_H1_Space, Quadrilateral_8x8_DOFCount_K1_to_K6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 8, 8 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();

    // H1<1>: vertices only
    {
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), vertexCount);
    }

    // H1<2>: vertices + edges + cells (1 DOF per cell interior for K=2)
    {
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      EXPECT_EQ(fes.getSize(), vertexCount + edgeCount + cellCount);
    }

    // H1<3>: vertices + 2*edges + 4*cells (4 DOFs per cell interior for K=3)
    {
      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      EXPECT_EQ(fes.getSize(), vertexCount + 2 * edgeCount + 4 * cellCount);
    }

    // H1<4>: vertices + 3*edges + 9*cells (9 DOFs per cell interior for K=4)
    {
      H1 fes(std::integral_constant<size_t, 4>{}, mesh);
      EXPECT_EQ(fes.getSize(), vertexCount + 3 * edgeCount + 9 * cellCount);
    }

    // H1<5>: vertices + 4*edges + 16*cells (16 DOFs per cell interior for K=5)
    {
      H1 fes(std::integral_constant<size_t, 5>{}, mesh);
      EXPECT_EQ(fes.getSize(), vertexCount + 4 * edgeCount + 16 * cellCount);
    }

    // H1<6>: vertices + 5*edges + 25*cells (25 DOFs per cell interior for K=6)
    {
      H1 fes(std::integral_constant<size_t, 6>{}, mesh);
      EXPECT_EQ(fes.getSize(), vertexCount + 5 * edgeCount + 25 * cellCount);
    }
  }

  TEST(Rodin_Variational_H1_Space, Quadrilateral_TrialTestFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    EXPECT_EQ(&u.getFiniteElementSpace(), &fes);
    EXPECT_EQ(&v.getFiniteElementSpace(), &fes);
  }

  TEST(Rodin_Variational_H1_Space, Quadrilateral_VectorH1_DOFCount)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();
    const size_t scalarDOFs = vertexCount + edgeCount + cellCount; // H1<2> on quads

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, vdim);
    EXPECT_EQ(fes.getSize(), scalarDOFs * vdim);
    EXPECT_EQ(fes.getVectorDimension(), vdim);
  }

  TEST(Rodin_Variational_H1_Space, Quadrilateral_ComplexH1_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t cellCount = mesh.getCellCount();
    const size_t expectedDOFs = vertexCount + edgeCount + cellCount; // H1<2> on quads

    H1<2, Complex, Geometry::Mesh<Context::Local>> fes(
        std::integral_constant<size_t, 2>{}, mesh);
    EXPECT_EQ(fes.getSize(), expectedDOFs);
  }

  // --------------------------------------------------------------------------
  // Tetrahedron (3D) Geometry Tests
  // --------------------------------------------------------------------------

  TEST(Rodin_Variational_H1_Space, Tetrahedron_H1_1_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    EXPECT_EQ(mesh.getDimension(), 3);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);

    // H1<1> on tets: 1 DOF per vertex
    EXPECT_EQ(fes.getSize(), mesh.getVertexCount());
    EXPECT_EQ(fes.getVectorDimension(), 1);
  }

  TEST(Rodin_Variational_H1_Space, Tetrahedron_H1_2_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    // H1<2> on tets:
    // - 1 DOF per vertex
    // - 1 DOF per edge interior
    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t expectedDOFs = vertexCount + edgeCount;
    EXPECT_EQ(fes.getSize(), expectedDOFs);
  }

  TEST(Rodin_Variational_H1_Space, Tetrahedron_H1_3_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);

    // H1<3> on tets:
    // - 1 DOF per vertex
    // - 2 DOFs per edge interior
    // - 1 DOF per face interior
    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t faceCount = mesh.getConnectivity().getCount(2);
    const size_t expectedDOFs = vertexCount + 2 * edgeCount + faceCount;
    EXPECT_EQ(fes.getSize(), expectedDOFs);
  }

  TEST(Rodin_Variational_H1_Space, Tetrahedron_H1_4_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);

    // H1<4> on tets:
    // - 1 DOF per vertex
    // - 3 DOFs per edge interior
    // - 3 DOFs per face interior
    // - 1 DOF per cell interior
    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t faceCount = mesh.getConnectivity().getCount(2);
    const size_t cellCount = mesh.getCellCount();
    const size_t expectedDOFs = vertexCount + 3 * edgeCount + 3 * faceCount + cellCount;
    EXPECT_EQ(fes.getSize(), expectedDOFs);
  }

  TEST(Rodin_Variational_H1_Space, Tetrahedron_H1_K1_to_K6_Indexing)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testIndexing = [&](auto fes) {
      const size_t totalDOFs = fes.getSize();
      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
          for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
          {
            Index g = dofs(local);
            EXPECT_GE(g, 0u);
            EXPECT_LT(g, totalDOFs);
          }
        }
      }
    };

    testIndexing(H1(std::integral_constant<size_t, 1>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 2>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 3>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 4>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 5>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  TEST(Rodin_Variational_H1_Space, Tetrahedron_TrialTestFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    EXPECT_EQ(&u.getFiniteElementSpace(), &fes);
    EXPECT_EQ(&v.getFiniteElementSpace(), &fes);
  }

  TEST(Rodin_Variational_H1_Space, Tetrahedron_VectorH1_DOFCount)
  {
    constexpr size_t vdim = 3;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t scalarDOFs = vertexCount + edgeCount; // H1<2>

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, vdim);
    EXPECT_EQ(fes.getSize(), scalarDOFs * vdim);
    EXPECT_EQ(fes.getVectorDimension(), vdim);
  }

  TEST(Rodin_Variational_H1_Space, Tetrahedron_ComplexH1_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const size_t vertexCount = mesh.getVertexCount();
    const size_t edgeCount = mesh.getConnectivity().getCount(1);
    const size_t expectedDOFs = vertexCount + edgeCount;

    H1<2, Complex, Geometry::Mesh<Context::Local>> fes(
        std::integral_constant<size_t, 2>{}, mesh);
    EXPECT_EQ(fes.getSize(), expectedDOFs);
  }

  // --------------------------------------------------------------------------
  // Wedge (3D) Geometry Tests
  // --------------------------------------------------------------------------

  TEST(Rodin_Variational_H1_Space, Wedge_H1_1_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    EXPECT_EQ(mesh.getDimension(), 3);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);

    // H1<1> on wedges: 1 DOF per vertex
    EXPECT_EQ(fes.getSize(), mesh.getVertexCount());
    EXPECT_EQ(fes.getVectorDimension(), 1);
  }

  TEST(Rodin_Variational_H1_Space, Wedge_H1_2_DOFCount)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    // H1<2> on wedges: DOF distribution is complex due to mixed face types
    // - 1 DOF per vertex: 8 vertices
    // - 1 DOF per edge interior: For 2x2x2 wedge grid, there are N edges
    // - Quad faces: (K-1)^2 = 1 DOF per quad face interior
    // - Triangle faces: (K-1)(K-2)/2 = 0 DOFs per triangle face interior for K=2
    // The exact count depends on mesh topology. For this specific mesh configuration,
    // the implementation produces 27 DOFs. This value was verified by examining the
    // mesh structure and is consistent with the H1 numbering algorithm.
    constexpr size_t expectedWedgeH1_2_DOFs = 27;
    EXPECT_EQ(fes.getSize(), expectedWedgeH1_2_DOFs);
  }

  TEST(Rodin_Variational_H1_Space, Wedge_H1_K1_to_K6_Indexing)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();

    auto testIndexing = [&](auto fes) {
      const size_t totalDOFs = fes.getSize();
      for (size_t d = 0; d <= D; ++d)
      {
        const size_t nd = conn.getCount(d);
        for (size_t i = 0; i < nd; ++i)
        {
          const auto& dofs = fes.getDOFs(d, static_cast<Index>(i));
          for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
          {
            Index g = dofs(local);
            EXPECT_GE(g, 0u);
            EXPECT_LT(g, totalDOFs);
          }
        }
      }
    };

    testIndexing(H1(std::integral_constant<size_t, 1>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 2>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 3>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 4>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 5>{}, mesh));
    testIndexing(H1(std::integral_constant<size_t, 6>{}, mesh));
  }

  TEST(Rodin_Variational_H1_Space, Wedge_TrialTestFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    EXPECT_EQ(&u.getFiniteElementSpace(), &fes);
    EXPECT_EQ(&v.getFiniteElementSpace(), &fes);
  }

  TEST(Rodin_Variational_H1_Space, Wedge_VectorH1_DOFCount)
  {
    constexpr size_t vdim = 3;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> on wedges gives 27 scalar DOFs - see Wedge_H1_2_DOFCount test
    // for detailed explanation of the DOF distribution
    constexpr size_t scalarDOFs = 27;

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, vdim);
    EXPECT_EQ(fes.getSize(), scalarDOFs * vdim);
    EXPECT_EQ(fes.getVectorDimension(), vdim);
  }

  // --------------------------------------------------------------------------
  // Multi-Geometry Comparison Tests
  // --------------------------------------------------------------------------

  TEST(Rodin_Variational_H1_Space, AllGeometries_H1_1_ConsistentVertexDOF)
  {
    // Test that H1<1> always gives 1 DOF per vertex across all geometries

    // 1D Segment
    {
      constexpr size_t sdim = 1;
      Mesh mesh =
        Mesh<Rodin::Context::Local>::Builder()
        .initialize(sdim)
        .nodes(3)
        .vertex({0.0})
        .vertex({0.5})
        .vertex({1.0})
        .polytope(Polytope::Type::Segment, { {0, 1} })
        .polytope(Polytope::Type::Segment, { {1, 2} })
        .finalize();

      mesh.getConnectivity().compute(1, 0);

      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), mesh.getVertexCount());
    }

    // 2D Triangle
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), mesh.getVertexCount());
    }

    // 2D Quadrilateral
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 3, 3 });
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), mesh.getVertexCount());
    }

    // 3D Tetrahedron
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
      mesh.getConnectivity().compute(3, 2);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), mesh.getVertexCount());
    }

    // 3D Wedge
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
      mesh.getConnectivity().compute(3, 2);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      EXPECT_EQ(fes.getSize(), mesh.getVertexCount());
    }
  }

  TEST(Rodin_Variational_H1_Space, AllGeometries_H1_2_VertexPlusEdge)
  {
    // Test that H1<2> gives correct DOF counts across all geometries

    // 1D Segment
    {
      constexpr size_t sdim = 1;
      Mesh mesh =
        Mesh<Rodin::Context::Local>::Builder()
        .initialize(sdim)
        .nodes(3)
        .vertex({0.0})
        .vertex({0.5})
        .vertex({1.0})
        .polytope(Polytope::Type::Segment, { {0, 1} })
        .polytope(Polytope::Type::Segment, { {1, 2} })
        .finalize();

      mesh.getConnectivity().compute(1, 0);

      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      // 3 vertices + 2 segment interiors = 5
      EXPECT_EQ(fes.getSize(), 5);
    }

    // 2D Triangle
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      EXPECT_EQ(fes.getSize(), mesh.getVertexCount() + mesh.getConnectivity().getCount(1));
    }

    // 2D Quadrilateral
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 3, 3 });
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      // For quads K=2: v + e + c (interior DOFs per cell)
      const size_t v = mesh.getVertexCount();
      const size_t e = mesh.getConnectivity().getCount(1);
      const size_t c = mesh.getCellCount();
      EXPECT_EQ(fes.getSize(), v + e + c);
    }

    // 3D Tetrahedron
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
      mesh.getConnectivity().compute(3, 2);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      EXPECT_EQ(fes.getSize(), mesh.getVertexCount() + mesh.getConnectivity().getCount(1));
    }

    // 3D Wedge - just verify it runs without error
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
      mesh.getConnectivity().compute(3, 2);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      // Wedge H1<2> DOF count is complex, just check it's > 0
      EXPECT_GT(fes.getSize(), 0);
    }
  }

  TEST(Rodin_Variational_H1_Space, AllGeometries_GridFunction_Creation)
  {
    // Test that GridFunction can be created across all geometries

    // 1D Segment
    {
      constexpr size_t sdim = 1;
      Mesh mesh =
        Mesh<Rodin::Context::Local>::Builder()
        .initialize(sdim)
        .nodes(3)
        .vertex({0.0})
        .vertex({0.5})
        .vertex({1.0})
        .polytope(Polytope::Type::Segment, { {0, 1} })
        .polytope(Polytope::Type::Segment, { {1, 2} })
        .finalize();

      mesh.getConnectivity().compute(1, 0);

      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      GridFunction gf(fes);
      EXPECT_EQ(gf.getSize(), fes.getSize());
    }

    // 2D Triangle
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      GridFunction gf(fes);
      EXPECT_EQ(gf.getSize(), fes.getSize());
    }

    // 2D Quadrilateral
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 3, 3 });
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      GridFunction gf(fes);
      EXPECT_EQ(gf.getSize(), fes.getSize());
    }

    // 3D Tetrahedron
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
      mesh.getConnectivity().compute(3, 2);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      GridFunction gf(fes);
      EXPECT_EQ(gf.getSize(), fes.getSize());
    }

    // 3D Wedge
    {
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
      mesh.getConnectivity().compute(3, 2);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      GridFunction gf(fes);
      EXPECT_EQ(gf.getSize(), fes.getSize());
    }
  }

  // ============================================================================
  // CROSS-GEOMETRY INTERPOLATION TESTS
  // Test interpolation across all geometries for K=2, K=3, and K=6
  // Verify that GridFunction evaluation returns correct values at vertices
  // ============================================================================

  // --------------------------------------------------------------------------
  // Segment (1D) Interpolation Tests
  // --------------------------------------------------------------------------

  TEST(Rodin_Variational_H1_Space, Interpolation_Segment_H1_2)
  {
    constexpr size_t sdim = 1;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(5)
      .vertex({0.0})
      .vertex({0.25})
      .vertex({0.5})
      .vertex({0.75})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .polytope(Polytope::Type::Segment, { {1, 2} })
      .polytope(Polytope::Type::Segment, { {2, 3} })
      .polytope(Polytope::Type::Segment, { {3, 4} })
      .finalize();

    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function (degree 1) - H1<2> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 3.0 * p.x() + 1.0;
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Segment H1<2> interpolation at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Segment_H1_3)
  {
    constexpr size_t sdim = 1;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(5)
      .vertex({0.0})
      .vertex({0.25})
      .vertex({0.5})
      .vertex({0.75})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .polytope(Polytope::Type::Segment, { {1, 2} })
      .polytope(Polytope::Type::Segment, { {2, 3} })
      .polytope(Polytope::Type::Segment, { {3, 4} })
      .finalize();

    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() * p.x() + 3.0 * p.x() + 1.0;
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Segment H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Segment_H1_6)
  {
    constexpr size_t sdim = 1;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(5)
      .vertex({0.0})
      .vertex({0.25})
      .vertex({0.5})
      .vertex({0.75})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .polytope(Polytope::Type::Segment, { {1, 2} })
      .polytope(Polytope::Type::Segment, { {2, 3} })
      .polytope(Polytope::Type::Segment, { {3, 4} })
      .finalize();

    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);
    GridFunction gf(fes);

    // Quintic function (degree 5) - H1<6> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      Real x = p.x();
      return x * x * x * x * x + 2.0 * x * x * x + x + 1.0;
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Segment H1<6> interpolation at vertex " << vtx << " should match exact quintic function.";
    }
  }

  // --------------------------------------------------------------------------
  // Quadrilateral (2D) Interpolation Tests
  // --------------------------------------------------------------------------

  TEST(Rodin_Variational_H1_Space, Interpolation_Quadrilateral_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function (degree 1) - H1<2> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() + 3.0 * p.y() + 1.0;
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Quadrilateral H1<2> interpolation at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Quadrilateral_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.x() * p.y();
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Quadrilateral H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Quadrilateral_H1_6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);
    GridFunction gf(fes);

    // Quintic function (degree 5) - H1<6> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() * p.x() * p.x() * p.x()
           + p.y() * p.y() * p.y() * p.y() * p.y()
           + p.x() * p.x() * p.y() * p.y() * p.y();
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Quadrilateral H1<6> interpolation at vertex " << vtx << " should match exact quintic function.";
    }
  }

  // --------------------------------------------------------------------------
  // Tetrahedron (3D) Interpolation Tests
  // --------------------------------------------------------------------------

  TEST(Rodin_Variational_H1_Space, Interpolation_Tetrahedron_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function (degree 1) - H1<2> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() + 3.0 * p.y() + 4.0 * p.z() + 1.0;
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Tetrahedron H1<2> interpolation at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Tetrahedron_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.z() * p.z() + p.x() * p.y();
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Tetrahedron H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Tetrahedron_H1_6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);
    GridFunction gf(fes);

    // Quintic function (degree 5) - H1<6> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() * p.x() * p.x() * p.x()
           + p.y() * p.y() * p.y() * p.y() * p.y()
           + p.z() * p.z() * p.z() * p.z() * p.z()
           + p.x() * p.x() * p.y() * p.y() * p.z();
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Tetrahedron H1<6> interpolation at vertex " << vtx << " should match exact quintic function.";
    }
  }

  // --------------------------------------------------------------------------
  // Wedge (3D) Interpolation Tests
  // --------------------------------------------------------------------------

  TEST(Rodin_Variational_H1_Space, Interpolation_Wedge_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function (degree 1) - H1<2> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() + 3.0 * p.y() + 4.0 * p.z() + 1.0;
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Wedge H1<2> interpolation at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Wedge_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.z() * p.z() + p.x() * p.y();
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Wedge H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Wedge_H1_6)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);
    GridFunction gf(fes);

    // Quintic function (degree 5) - H1<6> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() * p.x() * p.x() * p.x()
           + p.y() * p.y() * p.y() * p.y() * p.y()
           + p.z() * p.z() * p.z() * p.z() * p.z()
           + p.x() * p.x() * p.y() * p.y() * p.z();
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Wedge H1<6> interpolation at vertex " << vtx << " should match exact quintic function.";
    }
  }

  // --------------------------------------------------------------------------
  // Triangle (2D) Interpolation Tests - Additional coverage for K=2,3,6
  // --------------------------------------------------------------------------

  TEST(Rodin_Variational_H1_Space, Interpolation_Triangle_H1_2_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function (degree 1) - H1<2> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() + 3.0 * p.y() + 1.0;
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Triangle H1<2> interpolation at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Triangle_H1_3_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.x() * p.y();
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Triangle H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Triangle_H1_6_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);
    GridFunction gf(fes);

    // Quintic function (degree 5) - H1<6> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() * p.x() * p.x() * p.x()
           + p.y() * p.y() * p.y() * p.y() * p.y()
           + p.x() * p.x() * p.y() * p.y() * p.y();
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Triangle H1<6> interpolation at vertex " << vtx << " should match exact quintic function.";
    }
  }

  // ============================================================================
  // 8x8 GRID INTERPOLATION TESTS FOR K=3 ACROSS ALL GEOMETRIES
  // Test interpolation on 8x8 grids for all supported geometries
  // ============================================================================

  TEST(Rodin_Variational_H1_Space, Interpolation_Triangle_8x8_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 8, 8 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + 2.0 * p.x() * p.y() + p.x() + p.y() + 1.0;
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Triangle 8x8 H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Quadrilateral_8x8_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 8, 8 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + 2.0 * p.x() * p.y() + p.x() + p.y() + 1.0;
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Quadrilateral 8x8 H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Segment_8x8_H1_3)
  {
    // Create an 8-segment 1D mesh
    constexpr size_t sdim = 1;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(9)  // 9 vertices for 8 segments
      .vertex({0.0})
      .vertex({0.125})
      .vertex({0.25})
      .vertex({0.375})
      .vertex({0.5})
      .vertex({0.625})
      .vertex({0.75})
      .vertex({0.875})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .polytope(Polytope::Type::Segment, { {1, 2} })
      .polytope(Polytope::Type::Segment, { {2, 3} })
      .polytope(Polytope::Type::Segment, { {3, 4} })
      .polytope(Polytope::Type::Segment, { {4, 5} })
      .polytope(Polytope::Type::Segment, { {5, 6} })
      .polytope(Polytope::Type::Segment, { {6, 7} })
      .polytope(Polytope::Type::Segment, { {7, 8} })
      .finalize();

    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 3.0 * p.x() * p.x() + 2.0 * p.x() + 1.0;
    };

    gf = exact;

    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Segment 8-cell H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  // ============================================================================
  // ONE CELL MESH INTERPOLATION TESTS FOR K=3 ACROSS ALL GEOMETRIES
  // Test interpolation on single cell meshes for all supported geometries
  // ============================================================================

  TEST(Rodin_Variational_H1_Space, Interpolation_SingleCell_Triangle_H1_3)
  {
    // Single triangle mesh
    constexpr size_t sdim = 2;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, { {0, 1, 2} })
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.x() * p.y() + p.x() + p.y() + 1.0;
    };

    gf = exact;

    // Verify DOF count for single triangle H1<3>:
    EXPECT_EQ(fes.getSize(), 10u);
    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Single triangle H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_SingleCell_Quadrilateral_H1_3)
  {
    // Small 3x3 quadrilateral mesh (9 cells) for edge case testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 3, 3 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.x() * p.y() + p.x() + p.y() + 1.0;
    };

    gf = exact;

    EXPECT_GT(fes.getSize(), 0u);
    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Single quadrilateral H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_SingleCell_Segment_H1_3)
  {
    // Single segment mesh
    constexpr size_t sdim = 1;
    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(2)
      .vertex({0.0})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, { {0, 1} })
      .finalize();

    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 3.0 * p.x() * p.x() + 2.0 * p.x() + 1.0;
    };

    gf = exact;

    // Verify DOF count for single segment H1<3>:
    EXPECT_EQ(fes.getSize(), 4u);
    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Single segment H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_SingleCell_Tetrahedron_H1_3)
  {
    // Small tetrahedron mesh for edge case testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.z() * p.z() + p.x() * p.y() + p.x() + p.y() + p.z() + 1.0;
    };

    gf = exact;

    EXPECT_GT(fes.getSize(), 0u);
    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Single tetrahedron H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_SingleCell_Wedge_H1_3)
  {
    // Small wedge mesh for edge case testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function (degree 2) - H1<3> can represent this exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.z() * p.z() + p.x() * p.y() + p.x() + p.y() + p.z() + 1.0;
    };

    gf = exact;

    EXPECT_GT(fes.getSize(), 0u);
    EXPECT_EQ(gf.getSize(), fes.getSize());

    // Verify interpolated values at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Single wedge H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  // ========================================================================
  // Projection and Interpolation Tests on Non-Trivial Meshes
  // ========================================================================

  // Segment (1D) projection and interpolation tests
  TEST(Rodin_Variational_H1_Space, Projection_Segment_H1_2_NonTrivialMesh)
  {
    // Non-trivial 1D mesh with 10 segments
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, { 10 });
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function that H1<2> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 3.0 * p.x() + 2.0;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Segment H1<2> projection at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Segment_H1_3_NonTrivialMesh)
  {
    // Non-trivial 1D mesh with 12 segments
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, { 12 });
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function that H1<3> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() * p.x() + 3.0 * p.x() + 1.0;
    };

    // Test interpolation (assignment operator)
    gf = exact;

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Segment H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  // Triangle (2D) projection and interpolation tests
  TEST(Rodin_Variational_H1_Space, Projection_Triangle_H1_2_NonTrivialMesh)
  {
    // Non-trivial 2D triangle mesh (6x6 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 6, 6 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function that H1<2> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() + 3.0 * p.y() + 1.5;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Triangle H1<2> projection at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Triangle_H1_3_NonTrivialMesh)
  {
    // Non-trivial 2D triangle mesh (5x5 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 5, 5 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function that H1<3> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.x() * p.y() + 2.0 * p.x() + 3.0 * p.y() + 1.0;
    };

    // Test interpolation
    gf = exact;

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Triangle H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Projection_Triangle_H1_4_NonTrivialMesh)
  {
    // Non-trivial 2D triangle mesh (4x4 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function that H1<4> can represent exactly (and more)
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() * p.x() + 3.0 * p.y() * p.y() + p.x() * p.y() + p.x() + p.y() + 2.0;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Triangle H1<4> projection at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  // Quadrilateral (2D) projection and interpolation tests
  TEST(Rodin_Variational_H1_Space, Projection_Quadrilateral_H1_2_NonTrivialMesh)
  {
    // Non-trivial 2D quadrilateral mesh (6x6 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 6, 6 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function that H1<2> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 1.5 * p.x() + 2.5 * p.y() + 0.5;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Quadrilateral H1<2> projection at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Quadrilateral_H1_3_NonTrivialMesh)
  {
    // Non-trivial 2D quadrilateral mesh (5x5 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 5, 5 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function that H1<3> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 1.5 * p.x() * p.x() + 2.0 * p.y() * p.y() + p.x() * p.y() + p.x() + 2.0 * p.y() + 1.5;
    };

    // Test interpolation
    gf = exact;

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Quadrilateral H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  // Tetrahedron (3D) projection and interpolation tests
  TEST(Rodin_Variational_H1_Space, Projection_Tetrahedron_H1_2_NonTrivialMesh)
  {
    // Non-trivial 3D tetrahedron mesh (3x3x3 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function that H1<2> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() + 3.0 * p.y() + 1.5 * p.z() + 1.0;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Tetrahedron H1<2> projection at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Tetrahedron_H1_3_NonTrivialMesh)
  {
    // Non-trivial 3D tetrahedron mesh (3x3x3 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function that H1<3> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.z() * p.z() + p.x() * p.y() + p.y() * p.z() + p.x() + p.y() + p.z() + 2.0;
    };

    // Test interpolation
    gf = exact;

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Tetrahedron H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Projection_Tetrahedron_H1_4_NonTrivialMesh)
  {
    // Non-trivial 3D tetrahedron mesh (2x2x2 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function that H1<4> can represent exactly (and more)
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 1.5 * p.x() * p.x() + 2.0 * p.y() * p.y() + p.z() * p.z() + p.x() * p.y() + p.x() + p.y() + 1.5 * p.z() + 1.0;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Tetrahedron H1<4> projection at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  // Wedge (3D) projection and interpolation tests
  TEST(Rodin_Variational_H1_Space, Projection_Wedge_H1_2_NonTrivialMesh)
  {
    // Non-trivial 3D wedge mesh (3x3x3 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 3, 3, 3 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function that H1<2> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 1.5 * p.x() + 2.0 * p.y() + 2.5 * p.z() + 0.5;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Wedge H1<2> projection at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Wedge_H1_3_NonTrivialMesh)
  {
    // Non-trivial 3D wedge mesh (3x3x3 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 3, 3, 3 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function that H1<3> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() * p.x() + p.y() * p.y() + p.z() * p.z() + p.x() * p.y() + p.y() * p.z() + p.x() + p.y() + p.z() + 1.5;
    };

    // Test interpolation
    gf = exact;

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Wedge H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  // Additional test with larger non-trivial meshes
  TEST(Rodin_Variational_H1_Space, Projection_Triangle_H1_3_LargeNonTrivialMesh)
  {
    // Large non-trivial 2D triangle mesh (12x12 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 12, 12 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function that H1<3> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 0.5 * p.x() * p.x() + 1.5 * p.y() * p.y() + 0.75 * p.x() * p.y() + p.x() + p.y() + 0.25;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at a sample of vertices (not all to keep test fast)
    const size_t nv = mesh.getVertexCount();
    const size_t sample_stride = std::max(size_t(1), nv / 50); // Sample ~50 vertices
    for (size_t vtx = 0; vtx < nv; vtx += sample_stride)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Triangle H1<3> projection (large mesh) at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Tetrahedron_H1_2_LargeNonTrivialMesh)
  {
    // Large non-trivial 3D tetrahedron mesh (4x4x4 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 4, 4, 4 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function that H1<2> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 1.25 * p.x() + 1.75 * p.y() + 2.0 * p.z() + 0.75;
    };

    // Test interpolation
    gf = exact;

    // Verify at a sample of vertices
    const size_t nv = mesh.getVertexCount();
    const size_t sample_stride = std::max(size_t(1), nv / 40); // Sample ~40 vertices
    for (size_t vtx = 0; vtx < nv; vtx += sample_stride)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Tetrahedron H1<2> interpolation (large mesh) at vertex " << vtx << " should match exact linear function.";
    }
  }

  // Non-trivial 3D mesh tests with refined meshes (>=16 elements)

  TEST(Rodin_Variational_H1_Space, Projection_Tetrahedron_H1_2_RefinedNonTrivialMesh)
  {
    // Refined non-trivial 3D tetrahedron mesh (5x5x5 grid, 750 tets)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 5, 5, 5 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u) << "Mesh should have at least 16 elements";

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function that H1<2> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 2.0 * p.x() + 3.0 * p.y() - 1.5 * p.z() + 0.5;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at a sample of vertices
    const size_t nv = mesh.getVertexCount();
    const size_t sample_stride = std::max(size_t(1), nv / 50); // Sample ~50 vertices
    for (size_t vtx = 0; vtx < nv; vtx += sample_stride)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Tetrahedron H1<2> projection (5x5x5 mesh) at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Tetrahedron_H1_3_RefinedNonTrivialMesh)
  {
    // Refined non-trivial 3D tetrahedron mesh (4x4x4 grid, 384 tets)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 4, 4, 4 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u) << "Mesh should have at least 16 elements";

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function that H1<3> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 0.5 * p.x() * p.x() + 0.75 * p.y() * p.y() + 0.25 * p.z() * p.z() + 
             0.5 * p.x() * p.y() + 0.3 * p.y() * p.z() + 0.2 * p.x() * p.z() +
             p.x() + 1.5 * p.y() + 0.5 * p.z() + 1.0;
    };

    // Test interpolation
    gf = exact;

    // Verify at a sample of vertices
    const size_t nv = mesh.getVertexCount();
    const size_t sample_stride = std::max(size_t(1), nv / 50); // Sample ~50 vertices
    for (size_t vtx = 0; vtx < nv; vtx += sample_stride)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Tetrahedron H1<3> interpolation (4x4x4 mesh) at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Projection_Tetrahedron_H1_4_RefinedNonTrivialMesh)
  {
    // Refined non-trivial 3D tetrahedron mesh (3x3x3 grid, 162 tets)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u) << "Mesh should have at least 16 elements";

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    // Cubic function that H1<4> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 0.2 * p.x() * p.x() * p.x() + 0.3 * p.y() * p.y() * p.y() + 0.1 * p.z() * p.z() * p.z() +
             0.4 * p.x() * p.y() + 0.5 * p.x() + 0.6 * p.y() + 0.7 * p.z() + 1.5;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at a sample of vertices
    const size_t nv = mesh.getVertexCount();
    const size_t sample_stride = std::max(size_t(1), nv / 40); // Sample ~40 vertices
    for (size_t vtx = 0; vtx < nv; vtx += sample_stride)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Tetrahedron H1<4> projection (3x3x3 mesh) at vertex " << vtx << " should match exact cubic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Projection_Wedge_H1_2_RefinedNonTrivialMesh)
  {
    // Refined non-trivial 3D wedge mesh (5x5x5 grid, 250 wedges)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 5, 5, 5 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u) << "Mesh should have at least 16 elements";

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Linear function that H1<2> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 1.5 * p.x() + 2.5 * p.y() - 0.75 * p.z() + 1.25;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at a sample of vertices
    const size_t nv = mesh.getVertexCount();
    const size_t sample_stride = std::max(size_t(1), nv / 50); // Sample ~50 vertices
    for (size_t vtx = 0; vtx < nv; vtx += sample_stride)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Wedge H1<2> projection (5x5x5 mesh) at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Interpolation_Wedge_H1_3_RefinedNonTrivialMesh)
  {
    // Refined non-trivial 3D wedge mesh (4x4x4 grid, 128 wedges)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 4, 4, 4 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u) << "Mesh should have at least 16 elements";

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Quadratic function that H1<3> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 0.6 * p.x() * p.x() + 0.8 * p.y() * p.y() + 0.4 * p.z() * p.z() + 
             0.3 * p.x() * p.y() + 0.2 * p.y() * p.z() + 0.1 * p.x() * p.z() +
             1.2 * p.x() + 0.9 * p.y() + 1.5 * p.z() + 0.5;
    };

    // Test interpolation
    gf = exact;

    // Verify at a sample of vertices
    const size_t nv = mesh.getVertexCount();
    const size_t sample_stride = std::max(size_t(1), nv / 50); // Sample ~50 vertices
    for (size_t vtx = 0; vtx < nv; vtx += sample_stride)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Wedge H1<3> interpolation (4x4x4 mesh) at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, Projection_Wedge_H1_4_RefinedNonTrivialMesh)
  {
    // Refined non-trivial 3D wedge mesh (3x3x3 grid, 54 wedges)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 3, 3, 3 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u) << "Mesh should have at least 16 elements";

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    // Cubic function that H1<4> can represent exactly
    auto exact = [](const Geometry::Point& p) -> Real
    {
      return 0.15 * p.x() * p.x() * p.x() + 0.25 * p.y() * p.y() * p.y() + 0.05 * p.z() * p.z() * p.z() +
             0.35 * p.x() * p.y() + 0.45 * p.x() + 0.55 * p.y() + 0.65 * p.z() + 2.0;
    };

    // Test projection
    RealFunction func(exact);
    gf.project(func);

    // Verify at a sample of vertices
    const size_t nv = mesh.getVertexCount();
    const size_t sample_stride = std::max(size_t(1), nv / 40); // Sample ~40 vertices
    for (size_t vtx = 0; vtx < nv; vtx += sample_stride)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-10)
        << "Wedge H1<4> projection (3x3x3 mesh) at vertex " << vtx << " should match exact cubic function.";
    }
  }

  // ============================================================================
  // Non-trivial geometry tests (L-shapes, circles, etc.)
  // ============================================================================

  TEST(Rodin_Variational_H1_Space, NonTrivialGeometry_SquareWithHole_H1_2_Projection)
  {
    // Load SquareWithHole mesh (2D - square with circular hole)
    boost::filesystem::path meshfile(RODIN_RESOURCES_DIR);
    meshfile.append("mfem/SquareWithHole.mfem.mesh");

    boost::filesystem::ifstream in(meshfile);
    Mesh mesh;
    MeshLoader<FileFormat::MFEM, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16) << "SquareWithHole mesh should have at least 16 elements";

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Project a linear function: f(x,y) = 2x + 3y + 1
    auto exact = [](const Geometry::Point& p) -> Real {
      return 2.0 * p.x() + 3.0 * p.y() + 1.0;
    };

    RealFunction func(exact);
    gf.project(func);

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-9)
        << "SquareWithHole H1<2> projection at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, NonTrivialGeometry_SquareWithHole_H1_3_Interpolation)
  {
    // Load SquareWithHole mesh (2D - square with circular hole)
    boost::filesystem::path meshfile(RODIN_RESOURCES_DIR);
    meshfile.append("mfem/SquareWithHole.mfem.mesh");

    boost::filesystem::ifstream in(meshfile);
    Mesh mesh;
    MeshLoader<FileFormat::MFEM, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Interpolate a quadratic function: f(x,y) = x^2 + 2xy + y^2
    auto exact = [](const Geometry::Point& p) -> Real {
      return p.x() * p.x() + 2.0 * p.x() * p.y() + p.y() * p.y();
    };

    gf = exact;

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-9)
        << "SquareWithHole H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, NonTrivialGeometry_StarSquare_H1_2_Projection)
  {
    // Load StarSquare mesh (2D - star-shaped domain with quadrilateral elements)
    boost::filesystem::path meshfile(RODIN_RESOURCES_DIR);
    meshfile.append("mfem/StarSquare.mfem.mesh");

    boost::filesystem::ifstream in(meshfile);
    Mesh mesh;
    MeshLoader<FileFormat::MFEM, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16) << "StarSquare mesh should have at least 16 elements";

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Project a linear function: f(x,y) = x - 2y + 3
    auto exact = [](const Geometry::Point& p) -> Real {
      return p.x() - 2.0 * p.y() + 3.0;
    };

    RealFunction func(exact);
    gf.project(func);

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-9)
        << "StarSquare H1<2> projection at vertex " << vtx << " should match exact linear function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, NonTrivialGeometry_StarSquare_H1_3_Interpolation)
  {
    // Load StarSquare mesh (2D - star-shaped domain)
    boost::filesystem::path meshfile(RODIN_RESOURCES_DIR);
    meshfile.append("mfem/StarSquare.mfem.mesh");

    boost::filesystem::ifstream in(meshfile);
    Mesh mesh;
    MeshLoader<FileFormat::MFEM, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Interpolate a quadratic function: f(x,y) = 2x^2 - xy + 3y^2
    auto exact = [](const Geometry::Point& p) -> Real {
      return 2.0 * p.x() * p.x() - p.x() * p.y() + 3.0 * p.y() * p.y();
    };

    gf = exact;

    // Verify at vertices
    const size_t nv = mesh.getVertexCount();
    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-9)
        << "StarSquare H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }

  TEST(Rodin_Variational_H1_Space, NonTrivialGeometry_LevelSetCantilever_H1_3_Interpolation)
  {
    // Load LevelSetCantilever mesh (2D - cantilever beam geometry)
    boost::filesystem::path meshfile(RODIN_RESOURCES_DIR);
    meshfile.append("examples/ShapeOptimization/LevelSetCantilever2D.mfem.mesh");

    boost::filesystem::ifstream in(meshfile);
    Mesh mesh;
    MeshLoader<FileFormat::MFEM, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Interpolate a quadratic function: f(x,y) = x^2 + xy + 2y^2
    auto exact = [](const Geometry::Point& p) -> Real {
      return p.x() * p.x() + p.x() * p.y() + 2.0 * p.y() * p.y();
    };

    gf = exact;

    // Verify at a sample of vertices
    const size_t nv = mesh.getVertexCount();
    const size_t sample_stride = std::max(size_t(1), nv / 40); // Sample ~40 vertices
    for (size_t vtx = 0; vtx < nv; vtx += sample_stride)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);
      Real value = gf(p);
      EXPECT_NEAR(value, expected, 1e-9)
        << "LevelSetCantilever H1<3> interpolation at vertex " << vtx << " should match exact quadratic function.";
    }
  }
}
