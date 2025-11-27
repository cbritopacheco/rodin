#include <gtest/gtest.h>
#include <type_traits>
#include "Rodin/Test/Random.h"

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
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

    mesh.getConnectivity().compute(1, 2);

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

    mesh.getConnectivity().compute(1, 2);

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

    mesh.getConnectivity().compute(1, 2);

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
    mesh.getConnectivity().compute(1, 2);

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
    mesh.getConnectivity().compute(1, 2);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    const auto& fe = fes.getFiniteElement(2, 0);
    EXPECT_EQ(fe.getGeometry(), Polytope::Type::Triangle);
    // H1Element<2> on Triangle has 6 DOFs
    EXPECT_EQ(fe.getCount(), 6);
    EXPECT_EQ(fe.getOrder(), 2);
  }

  // Test DOF ownership (owned DOFs per entity) for H1<2>
  // Test DOF structure (closure on cells) for H1<2>
  TEST(Rodin_Variational_H1_Space, GetDOFs_H1_2)
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

    // Build edge–cell incidence so that H1 can enumerate edge entities
    mesh.getConnectivity().compute(1, 2);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);

    const auto& conn = mesh.getConnectivity();
    const size_t nv  = mesh.getVertexCount();
    const size_t ne  = conn.getCount(1);
    const size_t nc  = mesh.getCellCount();

    // 1) Vertices (d = 0): H1<2> has 1 DOF per vertex
    for (size_t v = 0; v < nv; ++v)
    {
      const auto& dofs_v = fes.getDOFs(0, v);
      EXPECT_EQ(dofs_v.size(), 1)
        << "Vertex " << v << " should have exactly 1 DOF for H1<2>.";
    }

    // 2) Edges (d = 1): H1<2> has (K-1) = 1 interior DOF per edge
    for (size_t e = 0; e < ne; ++e)
    {
      const auto& dofs_e = fes.getDOFs(1, e);
      EXPECT_EQ(dofs_e.size(), 1)
        << "Edge " << e << " should have exactly 1 interior DOF for H1<2>.";
    }

    // 3) Cells (d = 2): H1<2> triangle has closure = 3 vertex DOFs + 3 edge DOFs = 6 DOFs
    for (size_t c = 0; c < nc; ++c)
    {
      const auto& dofs_c = fes.getDOFs(2, c);
      const auto& fe_c   = fes.getFiniteElement(2, c);

      // Check that cell closure size matches element DOF count
      EXPECT_EQ(dofs_c.size(), fe_c.getCount())
        << "Cell " << c << " should have closure size equal to element DOF count.";

      EXPECT_EQ(dofs_c.size(), 6)
        << "Triangle " << c << " should have 6 DOFs in its closure for H1<2>.";
    }

    // 4) Consistency check: for one cell, its closure is exactly the union
    //    of its vertex DOFs and edge DOFs (no interior DOFs for K=2).
    {
      const size_t c = 0; // first cell
      const auto& dofs_c = fes.getDOFs(2, c);

      std::set<Index> closure_from_subentities;

      // vertices of cell c
      const auto& poly = conn.getPolytope(2, c);
      for (Index v : poly)
      {
        const auto& dofs_v = fes.getDOFs(0, v);
        for (size_t k = 0; k < static_cast<size_t>(dofs_v.size()); ++k)
          closure_from_subentities.insert(dofs_v(k));
      }

      // edges of cell c
      const auto& edges = conn.getIncidence({ 2, 1 }, c);
      for (Index e : edges)
      {
        const auto& dofs_e = fes.getDOFs(1, e);
        for (size_t k = 0; k < static_cast<size_t>(dofs_e.size()); ++k)
          closure_from_subentities.insert(dofs_e(k));
      }

      // No interior DOFs for H1<2> on triangles, so closure_from_subentities
      // should match dofs_c exactly.
      EXPECT_EQ(closure_from_subentities.size(), dofs_c.size());
      for (size_t k = 0; k < static_cast<size_t>(dofs_c.size()); ++k)
      {
        EXPECT_EQ(closure_from_subentities.count(dofs_c(k)), 1u)
          << "Cell closure DOF " << dofs_c(k)
          << " should come from a vertex or an edge DOF.";
      }
    }

    // 5) Global count still matches fes.getSize()
    //    (sum of owned DOFs: vertices + edges + cell interior)
    size_t total_owned = 0;

    // owned at vertices
    for (size_t v = 0; v < nv; ++v)
      total_owned += fes.getDOFs(0, v).size();

    // owned at edges
    for (size_t e = 0; e < ne; ++e)
      total_owned += fes.getDOFs(1, e).size();

    // owned at cells: for K=2 on triangles, 0 interior DOFs
    // nothing to add; getDOFs(2,c) is closure, not owned

    EXPECT_EQ(total_owned, fes.getSize());
  }

  // Test that H1<1> size matches P1 size
  TEST(Rodin_Variational_H1_Space, H1_1_MatchesP1_Size)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

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

    mesh.getConnectivity().compute(1, 2);

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
    mesh.getConnectivity().compute(1, 2);

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
    mesh.getConnectivity().compute(1, 2);

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
    mesh.getConnectivity().compute(1, 2);

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
    mesh.getConnectivity().compute(1, 2);

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

    mesh.getConnectivity().compute(1, 2);

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
    mesh.getConnectivity().compute(1, 2);

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

  // Behavioral test: interpolate a simple function and check values at nodes for H1<1>
  TEST(Rodin_Variational_H1_Space, Interpolation_Behavior_H1_1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(1, 2);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);

    // Trial/Test functions just to ensure assembly/interpolation pipeline works
    TrialFunction u(fes);
    TestFunction  v(fes);

    // GridFunction interpolation test: u(x, y) = x + 2y
    GridFunction gf(fes);

    auto exact = [](const Geometry::Point& p) -> Real
    {
      return p.x() + 2.0 * p.y();
    };

    gf = exact; // relies on existing interpolation operator in your framework

    // Check that DOF values match the exact function at vertex coordinates
    const size_t nv  = mesh.getVertexCount();

    for (size_t vtx = 0; vtx < nv; ++vtx)
    {
      const auto vit = mesh.getVertex(vtx);
      const Geometry::Point p(
          *vit,
          Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
          vit->getCoordinates());
      Real expected = exact(p);

      Real value = gf(p); // assuming scalar gf is indexed by global DOF (vertex) for H1<1>

      EXPECT_NEAR(value, expected, 1e-12)
        << "Interpolated value at vertex " << vtx
        << " should match exact function.";
    }
  }
}
