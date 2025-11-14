#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include <complex>
#include "Rodin/Variational/Pk.h"
#include "Rodin/Geometry.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  // Test Pk<2> (P2) space on simple mesh
  TEST(Rodin_Variational_PkSpace, SanityTest_P2_Triangle_Mesh)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});
    
    // Create P2 space
    Pk<2, Real> Vh(mesh);
    
    // Check that the space was created successfully
    EXPECT_GT(Vh.getSize(), 0);
    EXPECT_EQ(Vh.getVectorDimension(), 1);
    
    // Check that we can get the finite element
    const auto& elem = Vh.getFiniteElement(mesh.getDimension(), 0);
    EXPECT_EQ(elem.getCount(), 6);  // P2 triangle has 6 DOFs
  }

  // Test Pk<1> (P1) space matches P1 DOF count
  TEST(Rodin_Variational_PkSpace, Consistency_P1_Matches_Pk1)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {3, 3});
    
    // Create Pk<1> space
    Pk<1, Real> Vh_pk(mesh);
    
    // P1 space has one DOF per vertex
    EXPECT_EQ(Vh_pk.getSize(), mesh.getVertexCount());
    EXPECT_EQ(Vh_pk.getVectorDimension(), 1);
  }

  // Test Pk<0> (P0) space matches P0 DOF count
  TEST(Rodin_Variational_PkSpace, Consistency_P0_Matches_Pk0)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {4, 4});
    
    // Create Pk<0> space
    Pk<0, Real> Vh_pk(mesh);
    
    // P0 space has one DOF per cell
    EXPECT_EQ(Vh_pk.getSize(), mesh.getCellCount());
    EXPECT_EQ(Vh_pk.getVectorDimension(), 1);
  }

  // Test vector-valued Pk<2> space
  TEST(Rodin_Variational_PkSpace, SanityTest_VectorP2_Triangle_Mesh)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});
    
    // Create vector P2 space with dimension 2
    Pk<2, Math::Vector<Real>> Vh(mesh, 2);
    
    // Check that the space was created successfully
    EXPECT_GT(Vh.getSize(), 0);
    EXPECT_EQ(Vh.getVectorDimension(), 2);
    
    // Check that we can get the finite element
    const auto& elem = Vh.getFiniteElement(mesh.getDimension(), 0);
    EXPECT_EQ(elem.getCount(), 12);  // P2 triangle has 6 nodes × 2 components = 12 DOFs
  }

  // Test Pk<3> (P3) space
  TEST(Rodin_Variational_PkSpace, SanityTest_P3_Triangle_Mesh)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});
    
    // Create P3 space
    Pk<3, Real> Vh(mesh);
    
    // Check that the space was created successfully
    EXPECT_GT(Vh.getSize(), 0);
    EXPECT_EQ(Vh.getVectorDimension(), 1);
    
    // Check that we can get the finite element for a triangle
    const auto& elem = Vh.getFiniteElement(mesh.getDimension(), 0);
    EXPECT_EQ(elem.getCount(), 10);  // P3 triangle has 10 DOFs
  }

  // Test Pk<4> space
  TEST(Rodin_Variational_PkSpace, SanityTest_P4_Triangle_Mesh)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});
    
    // Create P4 space
    Pk<4, Real> Vh(mesh);
    
    // Check that the space was created successfully
    EXPECT_GT(Vh.getSize(), 0);
    EXPECT_EQ(Vh.getVectorDimension(), 1);
    
    // Check that we can get the finite element for a triangle
    const auto& elem = Vh.getFiniteElement(mesh.getDimension(), 0);
    EXPECT_EQ(elem.getCount(), 15);  // P4 triangle has 15 DOFs
  }

  // Test getting DOFs for a cell
  TEST(Rodin_Variational_PkSpace, GetDOFs_P2_Triangle)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});
    
    // Create P2 space
    Pk<2, Real> Vh(mesh);
    
    // Get DOFs for first cell
    const size_t dim = mesh.getDimension();
    const auto& dofs = Vh.getDOFs(dim, 0);
    
    // P2 triangle has 6 DOFs
    EXPECT_EQ(dofs.size(), 6);
    
    // DOF indices should be valid
    for (Index i = 0; i < dofs.size(); ++i)
    {
      EXPECT_GE(dofs(i), 0);
      EXPECT_LT(static_cast<size_t>(dofs(i)), Vh.getSize());
    }
  }

  // Test copy constructor
  TEST(Rodin_Variational_PkSpace, CopyConstructor_P2)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});
    
    // Create P2 space
    Pk<2, Real> Vh1(mesh);
    
    // Copy construct
    Pk<2, Real> Vh2(Vh1);
    
    // Should have same properties
    EXPECT_EQ(Vh1.getSize(), Vh2.getSize());
    EXPECT_EQ(Vh1.getVectorDimension(), Vh2.getVectorDimension());
  }

  // Test move constructor
  TEST(Rodin_Variational_PkSpace, MoveConstructor_P2)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});
    
    // Create P2 space
    Pk<2, Real> Vh1(mesh);
    size_t originalSize = Vh1.getSize();
    
    // Move construct
    Pk<2, Real> Vh2(std::move(Vh1));
    
    // Vh2 should have the original properties
    EXPECT_EQ(Vh2.getSize(), originalSize);
    EXPECT_EQ(Vh2.getVectorDimension(), 1);
  }

  // Test complex-valued Pk space
  TEST(Rodin_Variational_PkSpace, SanityTest_ComplexP2_Triangle)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});
    
    // Create complex P2 space
    Pk<2, Complex> Vh(mesh);
    
    // Check that the space was created successfully
    EXPECT_GT(Vh.getSize(), 0);
    EXPECT_EQ(Vh.getVectorDimension(), 1);
  }

  // Test Pk<5> and Pk<6> spaces
  TEST(Rodin_Variational_PkSpace, HighOrder_P5_P6_Triangle)
  {
    // Create a simple triangle mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {2, 2});
    
    // Test P5
    {
      Pk<5, Real> Vh(mesh);
      EXPECT_GT(Vh.getSize(), 0);
      const auto& elem = Vh.getFiniteElement(mesh.getDimension(), 0);
      EXPECT_EQ(elem.getCount(), 21);  // P5 triangle has 21 DOFs
    }
    
    // Test P6
    {
      Pk<6, Real> Vh(mesh);
      EXPECT_GT(Vh.getSize(), 0);
      const auto& elem = Vh.getFiniteElement(mesh.getDimension(), 0);
      EXPECT_EQ(elem.getCount(), 28);  // P6 triangle has 28 DOFs
    }
  }
}
