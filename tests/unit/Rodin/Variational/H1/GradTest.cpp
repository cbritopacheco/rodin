#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/H1/Grad.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_H1_Grad, ShapeFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    auto grad_u = Grad(u);
    auto grad_v = Grad(v);

    // Gradient of a scalar function should be a vector function
    // For 2D, gradient should have 2 components
  }

  TEST(Rodin_Variational_H1_Grad, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    auto grad_gf = Grad(gf);

    // Gradient of a scalar GridFunction should be a vector function
    EXPECT_EQ(&grad_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Grad, GridFunction_LinearFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Project a linear function: f(x,y) = 3x + 4y
    RealFunction linear_func([](const Geometry::Point& p) { return 3.0 * p.x() + 4.0 * p.y(); });
    gf.project(linear_func);

    auto grad_gf = Grad(gf);

    // For linear functions, gradient should be constant
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

  TEST(Rodin_Variational_H1_Grad, GridFunction_ConstantFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
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

  TEST(Rodin_Variational_H1_Grad, Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);
    gf = RealFunction(1.0);

    auto grad_gf = Grad(gf);
    auto copied = grad_gf.copy();

    EXPECT_NE(copied, nullptr);

    delete copied;
  }

  TEST(Rodin_Variational_H1_Grad, UsageInBilinearForm)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
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

  TEST(Rodin_Variational_H1_Grad, GridFunction_QuadraticFunction_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> (quadratic) can represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Project a quadratic function: f(x,y) = x^2 + y^2
    RealFunction quadratic_func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y(); 
    });
    gf.project(quadratic_func);

    auto grad_gf = Grad(gf);

    // Evaluate at a point in the center of the domain
    auto it = mesh.getPolytope(mesh.getDimension(), mesh.getCellCount() / 2);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    auto grad_value = grad_gf.getValue(p);

    // For f(x,y) = x^2 + y^2, grad f = [2x, 2y]
    // At the point, we expect approximately [2*x_phys, 2*y_phys]
    const auto& phys_coords = p.getPhysicalCoordinates();
    EXPECT_NEAR(grad_value(0), 2.0 * phys_coords(0), 1e-10);
    EXPECT_NEAR(grad_value(1), 2.0 * phys_coords(1), 1e-10);
  }

  TEST(Rodin_Variational_H1_Grad, GridFunction_QuadraticFunction_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<3> (cubic) can also represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Project a quadratic function: f(x,y) = x^2 + 2xy + y^2
    RealFunction quadratic_func([](const Geometry::Point& p) { 
      return p.x() * p.x() + 2.0 * p.x() * p.y() + p.y() * p.y(); 
    });
    gf.project(quadratic_func);

    auto grad_gf = Grad(gf);

    // Evaluate at several points
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    std::vector<Math::Vector<Real>> test_coords =
    {
      Math::Vector<Real>{{ 0.25, 0.25 }},
      Math::Vector<Real>{{ 0.5, 0.5 }},
      Math::Vector<Real>{{ 0.75, 0.25 }}
    };

    for (const auto& rc : test_coords)
    {
      Point p(polytope, rc);
      auto grad_value = grad_gf.getValue(p);
      
      // For f(x,y) = x^2 + 2xy + y^2, grad f = [2x + 2y, 2x + 2y]
      const auto& phys_coords = p.getPhysicalCoordinates();
      Real expected_grad_x = 2.0 * phys_coords(0) + 2.0 * phys_coords(1);
      Real expected_grad_y = 2.0 * phys_coords(0) + 2.0 * phys_coords(1);
      
      EXPECT_NEAR(grad_value(0), expected_grad_x, 1e-10);
      EXPECT_NEAR(grad_value(1), expected_grad_y, 1e-10);
    }
  }

  TEST(Rodin_Variational_H1_Grad, MultipleEvaluations)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
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

  TEST(Rodin_Variational_H1_Grad, DifferentPolynomialDegrees)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Linear function for testing
    RealFunction linear_func([](const Geometry::Point& p) { return 3.0 * p.x() + 4.0 * p.y(); });

    // Test with H1<1> (linear)
    {
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      GridFunction gf(fes);
      gf.project(linear_func);
      auto grad_gf = Grad(gf);
      
      auto it = mesh.getPolytope(mesh.getDimension(), 0);
      const auto& polytope = *it;
      Point p(polytope, Math::Vector<Real>{{0.5, 0.5}});
      auto grad_value = grad_gf.getValue(p);
      
      EXPECT_NEAR(grad_value(0), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_value(1), 4.0, RODIN_FUZZY_CONSTANT);
    }

    // Test with H1<2> (quadratic)
    {
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      GridFunction gf(fes);
      gf.project(linear_func);
      auto grad_gf = Grad(gf);
      
      auto it = mesh.getPolytope(mesh.getDimension(), 0);
      const auto& polytope = *it;
      Point p(polytope, Math::Vector<Real>{{0.5, 0.5}});
      auto grad_value = grad_gf.getValue(p);
      
      EXPECT_NEAR(grad_value(0), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_value(1), 4.0, RODIN_FUZZY_CONSTANT);
    }

    // Test with H1<3> (cubic)
    {
      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      GridFunction gf(fes);
      gf.project(linear_func);
      auto grad_gf = Grad(gf);
      
      auto it = mesh.getPolytope(mesh.getDimension(), 0);
      const auto& polytope = *it;
      Point p(polytope, Math::Vector<Real>{{0.5, 0.5}});
      auto grad_value = grad_gf.getValue(p);
      
      EXPECT_NEAR(grad_value(0), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_value(1), 4.0, RODIN_FUZZY_CONSTANT);
    }
  }
}
