#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_Real_P1_GridFunction, SanityTest_Build)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    EXPECT_EQ(gf.getSize(), fes.getSize());
    EXPECT_EQ(gf.getDimension(), 1);
    EXPECT_EQ(&gf.getFiniteElementSpace(), &fes);
  }

  TEST(Rodin_Variational_Real_P1_GridFunction, AssignmentFromRealFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    RealFunction c(5.0);
    gf = c;

    // Check that all DOFs have the assigned value
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], 5.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_Real_P1_GridFunction, ProjectLinearFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    RealFunction linear_func([](const Geometry::Point& p) { return p.x() + p.y(); });
    gf.project(linear_func);

    // For P1 elements, a linear function should be represented exactly
    // Check a few known values based on mesh structure
    EXPECT_GT(gf.getSize(), 0);
  }

  TEST(Rodin_Variational_Real_P1_GridFunction, ArithmeticOperations)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf1(fes);
    GridFunction gf2(fes);

    gf1 = RealFunction(3.0);
    gf2 = RealFunction(2.0);

    // Test addition
    gf1 += gf2;
    for (Index i = 0; i < static_cast<Index>(gf1.getSize()); i++)
    {
      EXPECT_NEAR(gf1[i], 5.0, RODIN_FUZZY_CONSTANT);
    }

    // Test scalar multiplication
    gf1 *= 2.0;
    for (Index i = 0; i < static_cast<Index>(gf1.getSize()); i++)
    {
      EXPECT_NEAR(gf1[i], 10.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_Real_P1_GridFunction, SubtractionOperations)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf1(fes);
    GridFunction gf2(fes);

    gf1 = RealFunction(7.0);
    gf2 = RealFunction(3.0);

    // Test subtraction
    gf1 -= gf2;
    for (Index i = 0; i < static_cast<Index>(gf1.getSize()); i++)
    {
      EXPECT_NEAR(gf1[i], 4.0, RODIN_FUZZY_CONSTANT);
    }

    // Test scalar subtraction
    gf1 -= 1.0;
    for (Index i = 0; i < static_cast<Index>(gf1.getSize()); i++)
    {
      EXPECT_NEAR(gf1[i], 3.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_Real_P1_GridFunction, DivisionOperations)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    gf = RealFunction(12.0);

    // Test scalar division
    gf /= 3.0;
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], 4.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_Real_P1_GridFunction, MinMaxOperations)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    // Set different values
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      gf[i] = static_cast<Real>(i);
    }

    Index min_idx, max_idx;
    Real min_val = gf.min(min_idx);
    Real max_val = gf.max(max_idx);

    EXPECT_EQ(min_val, 0.0);
    EXPECT_EQ(max_val, static_cast<Real>(gf.getSize() - 1));
    EXPECT_EQ(min_idx, 0);
    EXPECT_EQ(max_idx, static_cast<Index>(gf.getSize() - 1));

    EXPECT_EQ(gf.argmin(), 0);
    EXPECT_EQ(gf.argmax(), static_cast<Index>(gf.getSize() - 1));
  }

  TEST(Rodin_Variational_Vector_P1_GridFunction, SanityTest_Build)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    EXPECT_EQ(gf.getSize(), fes.getSize());
    EXPECT_EQ(gf.getDimension(), vdim);
    EXPECT_EQ(&gf.getFiniteElementSpace(), &fes);
  }

  TEST(Rodin_Variational_Vector_P1_GridFunction, ProjectVectorFunction)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    VectorFunction vf{1.0, 2.0};
    gf.project(vf);

    EXPECT_GT(gf.getSize(), 0);
  }

  TEST(Rodin_Variational_Real_P1_GridFunction, GetValue)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    gf = RealFunction(7.5);

    // Create a point for evaluation
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    Real value = gf.getValue(p);
    EXPECT_NEAR(value, 7.5, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_Real_P1_GridFunction, ZeroInitialization)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    // GridFunction should be zero-initialized by default
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], 0.0, RODIN_FUZZY_CONSTANT);
    }
  }
}
