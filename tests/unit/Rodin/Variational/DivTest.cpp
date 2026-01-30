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
    TrialFunction u(Vh);
    auto div_u = Div(u);

    auto cellIt = mesh.getCell(0);
    size_t dofs = div_u.getDOFs(*cellIt);

    // P1 triangle with 2 components has 3 DOFs/component × 2 components = 6 total DOFs
    EXPECT_EQ(dofs, 6);
  }

  TEST(Rodin_Variational_Div, ShapeFunction_getDOFs_Tetrahedron_P1)
  {
    constexpr size_t vdim = 3;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    P1 Vh(mesh, vdim);
    TrialFunction u(Vh);
    auto div_u = Div(u);

    auto cellIt = mesh.getCell(0);
    size_t dofs = div_u.getDOFs(*cellIt);

    // P1 tetrahedron with 3 components has 4 DOFs/component × 3 components = 12 total DOFs
    EXPECT_EQ(dofs, 12);
  }

  TEST(Rodin_Variational_Div, RandomCoordinates_LinearVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Create vector P1 space for the vector field
    P1 fes_vector(mesh, mesh.getSpaceDimension());
    GridFunction gf(fes_vector);

    // Project linear vector field u = (2x, 3y) with div(u) = 2 + 3 = 5
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << 2.0 * p.x(), 3.0 * p.y();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(2, linear_lambda);
    gf.project(linear_func);

    auto div_gf = Div(gf);

    // Test at 20 random points across different cells
    RandomFloat gen(0.0, 1.0);
    for (int test = 0; test < 20; test++)
    {
      Index cellIdx = gen() * (mesh.getCellCount() - 1);
      auto it = mesh.getPolytope(mesh.getDimension(), cellIdx);
      const auto& polytope = *it;

      Real x = gen();
      Real y = gen();
      if (x + y > 1.0) {
        x = 1.0 - x;
        y = 1.0 - y;
      }
      const Math::Vector<Real> rc{{x, y}};
      Point p(polytope, rc);

      Real div_value = div_gf.getValue(p);
      // Expected div = 2 + 3 = 5
      EXPECT_NEAR(div_value, 5.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_Div, RandomCoordinates_DivergenceFreeField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 5, 5 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Create vector P1 space for the vector field
    P1 fes_vector(mesh, mesh.getSpaceDimension());
    GridFunction gf(fes_vector);

    // Project divergence-free field u = (y, -x) with div(u) = 0
    auto divfree_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.y(), -p.x();
      return v;
    };
    VectorFunction<decltype(divfree_lambda)> divfree_func(2, divfree_lambda);
    gf.project(divfree_func);

    auto div_gf = Div(gf);

    // Test at 15 random points across different cells
    RandomFloat gen(0.0, 1.0);
    for (int test = 0; test < 15; test++)
    {
      Index cellIdx = gen() * (mesh.getCellCount() - 1);
      auto it = mesh.getPolytope(mesh.getDimension(), cellIdx);
      const auto& polytope = *it;

      Real x = gen();
      Real y = gen();
      if (x + y > 1.0) {
        x = 1.0 - x;
        y = 1.0 - y;
      }
      const Math::Vector<Real> rc{{x, y}};
      Point p(polytope, rc);

      Real div_value = div_gf.getValue(p);
      // Expected div = 0
      EXPECT_NEAR(div_value, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_Div, RandomCoordinates_Tetrahedron_LinearVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Create vector P1 space for the 3D vector field
    P1 fes_vector(mesh, mesh.getSpaceDimension());
    GridFunction gf(fes_vector);

    // Project linear vector field u = (x, 2y, 3z) with div(u) = 1 + 2 + 3 = 6
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(3);
      v << p.x(), 2.0 * p.y(), 3.0 * p.z();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(3, linear_lambda);
    gf.project(linear_func);

    auto div_gf = Div(gf);

    // Test at 15 random points across different cells
    RandomFloat gen(0.0, 1.0);
    for (int test = 0; test < 15; test++)
    {
      Index cellIdx = gen() * (mesh.getCellCount() - 1);
      auto it = mesh.getPolytope(mesh.getDimension(), cellIdx);
      const auto& polytope = *it;

      // Generate random barycentric coordinates for tetrahedron
      Real r1 = gen();
      Real r2 = gen();
      Real r3 = gen();
      Real sum = r1 + r2 + r3;
      if (sum > 1.0) {
        r1 /= sum;
        r2 /= sum;
        r3 /= sum;
      }
      const Math::Vector<Real> rc{{r1, r2, r3}};
      Point p(polytope, rc);

      Real div_value = div_gf.getValue(p);
      // Expected div = 1 + 2 + 3 = 6
      EXPECT_NEAR(div_value, 6.0, RODIN_FUZZY_CONSTANT);
    }
  }
}
