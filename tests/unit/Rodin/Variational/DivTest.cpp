#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_Div, VectorGridFunction_Construction)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    auto div_gf = Div(gf);

    EXPECT_EQ(&div_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_Div, VectorShapeFunction_Construction)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);

    auto div_u = Div(u);
    auto div_v = Div(v);
  }

  TEST(Rodin_Variational_Div, ConstantVectorField)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    // Project a constant vector field
    VectorFunction constant_field{5.0, 7.0};
    gf.project(constant_field);

    auto div_gf = Div(gf);

    // Divergence of a constant vector field should be zero
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    Real div_value = div_gf.getValue(p);
    EXPECT_NEAR(div_value, 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_Div, LinearVectorField)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });  // Use finer mesh
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    // Project a linear vector field: (x, y) - divergence should be 2
    VectorFunction linear_field{
      [](const Geometry::Point& p) { return p.x(); },
      [](const Geometry::Point& p) { return p.y(); }
    };
    gf.project(linear_field);

    auto div_gf = Div(gf);

    // For P1 elements, divergence should be approximately constant = 2
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    Real div_value = div_gf.getValue(p);
    EXPECT_NEAR(div_value, 2.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_Div, Copy)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    auto div_gf = Div(gf);
    auto copied = div_gf.copy();

    EXPECT_NE(copied, nullptr);

    delete copied;
  }

  TEST(Rodin_Variational_Div, UsageInLinearForm)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes_scalar(mesh);
    P1 fes_vector(mesh, vdim);

    TestFunction v_scalar(fes_scalar);
    GridFunction u_vector(fes_vector);

    LinearForm lf(v_scalar);

    // ∫ div(u) * v dx - this is a mixed formulation
    lf = Integral(Div(u_vector), v_scalar);

    // This should be constructible
    EXPECT_FALSE(lf.getIntegrators().empty());
  }

  TEST(Rodin_Variational_Div, UsageInBilinearForm)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);

    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);

    // ∫ div(u) * div(v) dx - part of mixed elasticity formulations
    bf = Integral(Div(u), Div(v));
    bf.assemble();

    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
  }

  TEST(Rodin_Variational_Div, 3D_VectorField)
  {
    constexpr size_t vdim = 3;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });  // 2D mesh with 3D vector field
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    // Project a 3D constant vector field
    VectorFunction constant_field{1.0, 2.0, 3.0};
    gf.project(constant_field);

    auto div_gf = Div(gf);

    // Test that divergence can be constructed
  }

  TEST(Rodin_Variational_Div, ZeroVectorField)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    // Zero vector field (default initialization)
    auto div_gf = Div(gf);

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    Real div_value = div_gf.getValue(p);
    EXPECT_NEAR(div_value, 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_Div, MultipleEvaluations)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    // Project a divergence-free field: (-y, x)
    VectorFunction divergence_free{
      [](const Geometry::Point& p) { return -p.y(); },
      [](const Geometry::Point& p) { return p.x(); }
    };
    gf.project(divergence_free);

    auto div_gf = Div(gf);

    // Test multiple evaluation points
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;

    std::vector<Math::Vector<Real>> test_coords = {
      Math::Vector<Real>{{0.0, 0.0}},
      Math::Vector<Real>{{1.0, 0.0}},
      Math::Vector<Real>{{0.0, 1.0}},
      Math::Vector<Real>{{0.5, 0.5}},
      Math::Vector<Real>{{0.25, 0.75}}
    };

    for (const auto& rc : test_coords)
    {
      Point p(polytope, rc);
      Real div_value = div_gf.getValue(p);
  
      // Divergence should be zero everywhere for this field
      EXPECT_NEAR(div_value, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_Div, GetOperand)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    auto div_gf = Div(gf);
    const auto& operand = div_gf.getOperand();

    EXPECT_EQ(&operand, &gf);
  }

  TEST(Rodin_Variational_Div, CopyConstructor)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    gf = VectorFunction{3.0, 4.0};

    auto div_gf = Div(gf);
    auto div_copy(div_gf);

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    Real div_value_orig = div_gf.getValue(p);
    Real div_value_copy = div_copy.getValue(p);

    EXPECT_NEAR(div_value_orig, div_value_copy, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_Div, MoveConstructor)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    GridFunction gf(fes);

    gf = VectorFunction{1.5, 2.5};

    auto div_gf = Div(gf);

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);

    Real original_value = div_gf.getValue(p);

    auto div_moved(std::move(div_gf));
    EXPECT_NEAR(div_moved.getValue(p), original_value, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_Div, ShapeFunction_getDOFs_Triangle_P1)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    P1 Vh(mesh, vdim);
    ShapeFunction u(Vh);
    auto div_u = Div(u);

    auto cellIt = mesh.getCell(0);
    size_t dofs = div_u.getDOFs(*cellIt);

    // P1 triangle with 2 components has 3 DOFs/component × 2 components = 6 total DOFs
    EXPECT_EQ(dofs, 6);
  }

  TEST(Rodin_Variational_Div, ShapeFunction_getBasis_Triangle_P1)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    P1 Vh(mesh, vdim);
    ShapeFunction u(Vh);
    auto div_u = Div(u);

    auto cellIt = mesh.getCell(0);
    size_t dofs = div_u.getDOFs(*cellIt);

    Math::SpatialPoint<Real> sp(2);
    sp(0) = 0.3;
    sp(1) = 0.4;
    Point p(*cellIt, sp);
    div_u.setLeaf(p);

    // Check that getBasis returns valid scalar values for all local DOFs
    for (size_t local = 0; local < dofs; local++)
    {
      Real basis_val = div_u.getBasis(local);
      // Check value is finite
      EXPECT_TRUE(std::isfinite(basis_val));
    }
  }

  TEST(Rodin_Variational_Div, ShapeFunction_getDOFs_Tetrahedron_P1)
  {
    constexpr size_t vdim = 3;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    P1 Vh(mesh, vdim);
    ShapeFunction u(Vh);
    auto div_u = Div(u);

    auto cellIt = mesh.getCell(0);
    size_t dofs = div_u.getDOFs(*cellIt);

    // P1 tetrahedron with 3 components has 4 DOFs/component × 3 components = 12 total DOFs
    EXPECT_EQ(dofs, 12);
  }

  TEST(Rodin_Variational_Div, ShapeFunction_getBasis_Tetrahedron_P1)
  {
    constexpr size_t vdim = 3;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    P1 Vh(mesh, vdim);
    ShapeFunction u(Vh);
    auto div_u = Div(u);

    auto cellIt = mesh.getCell(0);
    size_t dofs = div_u.getDOFs(*cellIt);

    Math::SpatialPoint<Real> sp(3);
    sp(0) = 0.2;
    sp(1) = 0.2;
    sp(2) = 0.2;
    Point p(*cellIt, sp);
    div_u.setLeaf(p);

    // Check that getBasis returns valid scalar values for all local DOFs
    for (size_t local = 0; local < dofs; local++)
    {
      Real basis_val = div_u.getBasis(local);
      // Check value is finite
      EXPECT_TRUE(std::isfinite(basis_val));
    }
  }
}
