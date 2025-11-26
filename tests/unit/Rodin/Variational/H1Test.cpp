#include <gtest/gtest.h>
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

    H1<1, Real> fes(mesh);

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

    H1<2, Real> fes(mesh);

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

    H1<3, Real> fes(mesh);

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

    H1<2, Real> fes(mesh);

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

    H1<2, Real> fes(mesh);

    const auto& fe = fes.getFiniteElement(2, 0);
    EXPECT_EQ(fe.getGeometry(), Polytope::Type::Triangle);
    // H1Element<2> on Triangle has 6 DOFs
    EXPECT_EQ(fe.getCount(), 6);
    EXPECT_EQ(fe.getOrder(), 2);
  }

  // Test DOF retrieval
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

    mesh.getConnectivity().compute(1, 2);

    H1<2, Real> fes(mesh);

    // Get DOFs for first cell
    const auto& dofs = fes.getDOFs(2, 0);
    // Triangle with H1<2> should have 6 DOFs (3 vertices + 3 edges)
    EXPECT_EQ(dofs.size(), 6);
  }

  // Test that H1<1> size matches P1 size
  TEST(Rodin_Variational_H1_Space, H1_1_MatchesP1_Size)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    H1<1, Real> h1_fes(mesh);
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

    H1<2, Math::Vector<Real>> fes(mesh, vdim);

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

    H1<1, Real> fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
  }

  // Test that H1<2> can be used with TrialFunction and TestFunction
  TEST(Rodin_Variational_H1_Space, H1_2_TrialTestFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    H1<2, Real> fes(mesh);
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
    H1<1, Real> h1_1(mesh);
    EXPECT_EQ(h1_1.getSize(), vertexCount);

    // H1<2>: vertices + edges
    H1<2, Real> h1_2(mesh);
    EXPECT_EQ(h1_2.getSize(), vertexCount + edgeCount);

    // H1<3>: vertices + 2*edges + faces
    H1<3, Real> h1_3(mesh);
    EXPECT_EQ(h1_3.getSize(), vertexCount + 2 * edgeCount + cellCount);
  }
}
