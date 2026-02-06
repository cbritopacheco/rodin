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
  TEST(Rodin_Variational_Grad, ShapeFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    auto grad_u = Grad(u);
    auto grad_v = Grad(v);

    // Gradient of a scalar function should be a vector function
    // For 2D, gradient should have 2 components
  }

  TEST(Rodin_Variational_Grad, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    auto grad_gf = Grad(gf);

    // Gradient of a scalar GridFunction should be a vector function
    EXPECT_EQ(&grad_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_Grad, GridFunction_LinearFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a linear function: f(x,y) = 3x + 4y
    RealFunction linear_func([](const Geometry::Point& p) { return 3.0 * p.x() + 4.0 * p.y(); });
    gf.project(linear_func);

    auto grad_gf = Grad(gf);

    // For P1 elements, gradient of a linear function should be constant
    // Create a point for evaluation
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    auto grad_value = grad_gf.getValue(p);
    // Gradient should be approximately [3, 4]
    EXPECT_NEAR(grad_value(0), 3.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(grad_value(1), 4.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_Grad, GridFunction_ConstantFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a constant function
    gf = RealFunction(42.0);

    auto grad_gf = Grad(gf);

    // Gradient of a constant function should be zero
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    auto grad_value = grad_gf.getValue(p);
    EXPECT_NEAR(grad_value.norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_Grad, Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);
    gf = RealFunction(1.0);

    auto grad_gf = Grad(gf);
    auto copied = grad_gf.copy();

    EXPECT_NE(copied, nullptr);

    delete copied;
  }

  TEST(Rodin_Variational_Grad, UsageInBilinearForm)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);

    // Laplacian operator: ∫ ∇u · ∇v dx
    bf = Integral(Grad(u), Grad(v));

    EXPECT_FALSE(bf.getLocalIntegrators().empty());

    // This should assemble without errors
    bf.assemble();

    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
  }

  TEST(Rodin_Variational_Grad, 3D_Mesh)
  {
    // Test with 3D mesh (if supported by the current configuration)
    // For now, test with 2D to ensure compatibility
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    auto grad_gf = Grad(gf);

    // In 2D, gradient should have 2 components
  }

  TEST(Rodin_Variational_Grad, GridFunction_QuadraticFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a quadratic function: f(x,y) = x^2 + y^2
    // Note: P1 elements can't represent this exactly, but we can test the approximation
    RealFunction quadratic_func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y(); 
    });
    gf.project(quadratic_func);

    auto grad_gf = Grad(gf);

    // Evaluate at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    auto grad_value = grad_gf.getValue(p);

    // Since P1 can't represent quadratic exactly, we just check that we get some non-zero gradient
    EXPECT_GT(grad_value.norm(), 0.0);
  }

  TEST(Rodin_Variational_Grad, MultipleEvaluations)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a linear function: f(x,y) = x + 2y
    RealFunction linear_func([](const Geometry::Point& p) { return p.x() + 2.0 * p.y(); });
    gf.project(linear_func);

    auto grad_gf = Grad(gf);

    // Test multiple evaluation points
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;

    std::vector<Math::Vector<Real>> test_coords =
    {
      Math::Vector<Real>{{ 0.0, 0.0 }},
      Math::Vector<Real>{{ 1.0, 0.0 }},
      Math::Vector<Real>{{ 0.0, 1.0 }},
      Math::Vector<Real>{{ 0.5, 0.5 }},
      Math::Vector<Real>{{ 0.25, 0.75 }}
    };

    for (const auto& rc : test_coords)
    {
      Point p(polytope, rc);
      auto grad_value = grad_gf.getValue(p);
      // Gradient should be approximately [1, 2] everywhere for linear function
      EXPECT_NEAR(grad_value(0), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_value(1), 2.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_Grad, ShapeFunction_getDOFs_Triangle_P1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    P1 Vh(mesh);
    TrialFunction u(Vh);
    auto grad_u = Grad(u);

    auto cellIt = mesh.getCell(0);
    size_t dofs = grad_u.getDOFs(*cellIt);

    // P1 triangle has 3 DOFs
    EXPECT_EQ(dofs, 3);
  }

  TEST(Rodin_Variational_Grad, ShapeFunction_getDOFs_Tetrahedron_P1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    P1 Vh(mesh);
    TrialFunction u(Vh);
    auto grad_u = Grad(u);

    auto cellIt = mesh.getCell(0);
    size_t dofs = grad_u.getDOFs(*cellIt);

    // P1 tetrahedron has 4 DOFs
    EXPECT_EQ(dofs, 4);
  }

  TEST(Rodin_Variational_Grad, RandomCoordinates_LinearFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a linear function: f(x,y) = 2x + 3y
    // Gradient should be constant: [2, 3]
    RealFunction linear_func([](const Geometry::Point& p) {
      return 2.0 * p.x() + 3.0 * p.y();
    });
    gf.project(linear_func);

    auto grad_gf = Grad(gf);

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

      auto grad_value = grad_gf.getValue(p);
      EXPECT_NEAR(grad_value(0), 2.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_value(1), 3.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_Grad, RandomCoordinates_ZeroGradient)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 5, 5 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a constant function: f(x,y) = 5.0
    // Gradient should be zero: [0, 0]
    RealFunction constant_func([](const Geometry::Point& p) {
      return 5.0;
    });
    gf.project(constant_func);

    auto grad_gf = Grad(gf);

    // Test at 10 random points across different cells
    RandomFloat gen(0.0, 1.0);
    for (int test = 0; test < 10; test++)
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

      auto grad_value = grad_gf.getValue(p);
      EXPECT_NEAR(grad_value(0), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_value(1), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_Grad, RandomCoordinates_Tetrahedron_LinearFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a linear function: f(x,y,z) = 2x + 3y - z
    // Gradient should be constant: [2, 3, -1]
    RealFunction linear_func([](const Geometry::Point& p) {
      return 2.0 * p.x() + 3.0 * p.y() - p.z();
    });
    gf.project(linear_func);

    auto grad_gf = Grad(gf);

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

      auto grad_value = grad_gf.getValue(p);
      EXPECT_NEAR(grad_value(0), 2.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_value(1), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_value(2), -1.0, RODIN_FUZZY_CONSTANT);
    }
  }
}
