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
    const Math::Vector<Real> rc{{1.0/3.0, 1.0/3.0}};
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
      Math::Vector<Real>{{ 0.2, 0.2 }},
      Math::Vector<Real>{{ 0.6, 0.2 }},
      Math::Vector<Real>{{ 0.2, 0.6 }}
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
      Math::Vector<Real>{{ 0.2, 0.2 }},
      Math::Vector<Real>{{ 0.6, 0.2 }},
      Math::Vector<Real>{{ 0.2, 0.6 }},
      Math::Vector<Real>{{ 1.0/3.0, 1.0/3.0 }},
      Math::Vector<Real>{{ 0.1, 0.7 }} // sum=0.8
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

  TEST(Rodin_Variational_H1_Grad_Tet, ShapeFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    auto grad_u = Grad(u);
    auto grad_v = Grad(v);

    // No explicit asserts here: compilation + construction is the test.
    (void)grad_u;
    (void)grad_v;
  }

  TEST(Rodin_Variational_H1_Grad_Tet, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    auto grad_gf = Grad(gf);

    EXPECT_EQ(&grad_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Grad_Tet, GridFunction_ConstantFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    gf = RealFunction(42.0);

    auto grad_gf = Grad(gf);

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& cell = *it;

    // Interior reference point (sum < 1)
    const Math::Vector<Real> rc{{ 0.2, 0.2, 0.2 }};
    Point p(cell, rc);

    auto grad_value = grad_gf.getValue(p);
    EXPECT_NEAR(grad_value.norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Grad_Tet, GridFunction_LinearFunction_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> can represent linear functions exactly.
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // f(x,y,z) = 3x + 4y + 5z
    RealFunction linear_func([](const Geometry::Point& p)
    {
      return 3.0 * p.x() + 4.0 * p.y() + 5.0 * p.z();
    });
    gf.project(linear_func);

    auto grad_gf = Grad(gf);

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& cell = *it;

    const Math::Vector<Real> rc{{ 0.25, 0.25, 0.25 }};
    Point p(cell, rc);

    auto grad_value = grad_gf.getValue(p);

    EXPECT_NEAR(grad_value(0), 3.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(grad_value(1), 4.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(grad_value(2), 5.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Grad_Tet, GridFunction_QuadraticFunction_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> can represent quadratics exactly.
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // f(x,y,z) = x^2 + y^2 + z^2
    RealFunction quad([](const Geometry::Point& p)
    {
      return p.x() * p.x() + p.y() * p.y() + p.z() * p.z();
    });
    gf.project(quad);

    auto grad_gf = Grad(gf);

    auto it = mesh.getPolytope(mesh.getDimension(), mesh.getCellCount() / 2);
    const auto& cell = *it;

    const Math::Vector<Real> rc{{ 0.2, 0.3, 0.1 }}; // sum = 0.6 < 1
    Point p(cell, rc);

    auto grad_value = grad_gf.getValue(p);

    const auto& X = p.getPhysicalCoordinates();
    EXPECT_NEAR(grad_value(0), 2.0 * X(0), 1e-10);
    EXPECT_NEAR(grad_value(1), 2.0 * X(1), 1e-10);
    EXPECT_NEAR(grad_value(2), 2.0 * X(2), 1e-10);
  }

  TEST(Rodin_Variational_H1_Grad_Tet, GridFunction_QuadraticFunction_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<3> can also represent quadratics exactly.
    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // f(x,y,z) = x^2 + 2xy + y^2 + 3z^2
    RealFunction quad([](const Geometry::Point& p)
    {
      return p.x()*p.x() + 2.0*p.x()*p.y() + p.y()*p.y() + 3.0*p.z()*p.z();
    });
    gf.project(quad);

    auto grad_gf = Grad(gf);

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& cell = *it;

    std::vector<Math::Vector<Real>> test_coords =
    {
      Math::Vector<Real>{{ 0.1, 0.1, 0.1 }}, // sum=0.3
      Math::Vector<Real>{{ 0.2, 0.3, 0.1 }}, // sum=0.6
      Math::Vector<Real>{{ 0.3, 0.1, 0.2 }}  // sum=0.6
    };

    for (const auto& rc : test_coords)
    {
      Point p(cell, rc);
      auto grad_value = grad_gf.getValue(p);

      const auto& X = p.getPhysicalCoordinates();

      // ∇f = [2x + 2y, 2x + 2y, 6z]
      const Real ex = 2.0 * X(0) + 2.0 * X(1);
      const Real ey = 2.0 * X(0) + 2.0 * X(1);
      const Real ez = 6.0 * X(2);

      EXPECT_NEAR(grad_value(0), ex, 1e-10);
      EXPECT_NEAR(grad_value(1), ey, 1e-10);
      EXPECT_NEAR(grad_value(2), ez, 1e-10);
    }
  }

  TEST(Rodin_Variational_H1_Grad_Tet, MultipleEvaluations)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // f(x,y,z) = x + 2y + 3z  -> ∇f = [1,2,3]
    RealFunction linear_func([](const Geometry::Point& p)
    {
      return p.x() + 2.0 * p.y() + 3.0 * p.z();
    });
    gf.project(linear_func);

    auto grad_gf = Grad(gf);

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& cell = *it;

    std::vector<Math::Vector<Real>> test_coords =
    {
      Math::Vector<Real>{{ 0.1, 0.1, 0.1 }},
      Math::Vector<Real>{{ 0.2, 0.2, 0.1 }},
      Math::Vector<Real>{{ 0.3, 0.1, 0.1 }},
      Math::Vector<Real>{{ 0.2, 0.3, 0.1 }},
      Math::Vector<Real>{{ 0.1, 0.2, 0.3 }}
    };

    for (const auto& rc : test_coords)
    {
      ASSERT_LT(rc(0) + rc(1) + rc(2), 1.0); // stay inside reference tet
      Point p(cell, rc);
      auto grad_value = grad_gf.getValue(p);

      EXPECT_NEAR(grad_value(0), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_value(1), 2.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_value(2), 3.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Grad_Tet, DifferentPolynomialDegrees)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // f(x,y,z) = 3x + 4y + 5z
    RealFunction linear_func([](const Geometry::Point& p)
    {
      return 3.0 * p.x() + 4.0 * p.y() + 5.0 * p.z();
    });

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& cell = *it;
    Point p(cell, Math::Vector<Real>{{ 0.25, 0.25, 0.25 }});

    // H1<1>
    {
      H1 fes(std::integral_constant<size_t, 1>{}, mesh);
      GridFunction gf(fes);
      gf.project(linear_func);
      auto grad_gf = Grad(gf);
      auto g = grad_gf.getValue(p);

      EXPECT_NEAR(g(0), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(g(1), 4.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(g(2), 5.0, RODIN_FUZZY_CONSTANT);
    }

    // H1<2>
    {
      H1 fes(std::integral_constant<size_t, 2>{}, mesh);
      GridFunction gf(fes);
      gf.project(linear_func);
      auto grad_gf = Grad(gf);
      auto g = grad_gf.getValue(p);

      EXPECT_NEAR(g(0), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(g(1), 4.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(g(2), 5.0, RODIN_FUZZY_CONSTANT);
    }

    // H1<3>
    {
      H1 fes(std::integral_constant<size_t, 3>{}, mesh);
      GridFunction gf(fes);
      gf.project(linear_func);
      auto grad_gf = Grad(gf);
      auto g = grad_gf.getValue(p);

      EXPECT_NEAR(g(0), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(g(1), 4.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(g(2), 5.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Grad, RandomCoordinates_LinearFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
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

  TEST(Rodin_Variational_H1_Grad, RandomCoordinates_QuadraticFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 5, 5 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> can represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Project a quadratic function: f(x,y) = x^2 + y^2
    // Gradient should be: [2x, 2y]
    RealFunction quadratic_func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y(); 
    });
    gf.project(quadratic_func);

    auto grad_gf = Grad(gf);

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

      auto grad_value = grad_gf.getValue(p);
      const auto& phys_coords = p.getPhysicalCoordinates();
      
      // Expected gradient: [2x, 2y]
      EXPECT_NEAR(grad_value(0), 2.0 * phys_coords(0), 1e-10);
      EXPECT_NEAR(grad_value(1), 2.0 * phys_coords(1), 1e-10);
    }
  }

  TEST(Rodin_Variational_H1_Grad, RandomCoordinates_H1_3_QuadraticFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<3> can also represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Project: f(x,y) = x^2 - xy + 2y^2
    // Gradient: [2x - y, -x + 4y]
    RealFunction quadratic_func([](const Geometry::Point& p) { 
      return p.x() * p.x() - p.x() * p.y() + 2.0 * p.y() * p.y(); 
    });
    gf.project(quadratic_func);

    auto grad_gf = Grad(gf);

    // Test at 25 random points
    RandomFloat gen(0.0, 1.0);
    for (int test = 0; test < 25; test++)
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
      const auto& phys_coords = p.getPhysicalCoordinates();
      
      Real px = phys_coords(0);
      Real py = phys_coords(1);
      
      // Expected gradient: [2x - y, -x + 4y]
      EXPECT_NEAR(grad_value(0), 2.0 * px - py, 1e-10);
      EXPECT_NEAR(grad_value(1), -px + 4.0 * py, 1e-10);
    }
  }

  TEST(Rodin_Variational_H1_Grad, RandomCoordinates_ZeroGradient)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Project a constant function
    gf = RealFunction(7.5);

    auto grad_gf = Grad(gf);

    // Test at 10 random points - gradient should always be zero
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
      EXPECT_NEAR(grad_value.norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Grad, RandomCoordinates_H1_5_QuarticFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<5> can represent quartic (degree 4) functions exactly
    H1 fes(std::integral_constant<size_t, 5>{}, mesh);
    GridFunction gf(fes);

    // Project a quartic function: f(x,y) = x^4 + x^2*y^2 + y^4
    // Gradient: [4x^3 + 2xy^2, 2x^2*y + 4y^3]
    RealFunction quartic_func([](const Geometry::Point& p) { 
      Real x = p.x();
      Real y = p.y();
      return x*x*x*x + x*x*y*y + y*y*y*y; 
    });
    gf.project(quartic_func);

    auto grad_gf = Grad(gf);

    // Test at 30 random points across different cells
    RandomFloat gen(0.0, 1.0);
    for (int test = 0; test < 30; test++)
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
      const auto& phys_coords = p.getPhysicalCoordinates();
      
      Real px = phys_coords(0);
      Real py = phys_coords(1);
      
      // Expected gradient: [4x^3 + 2xy^2, 2x^2*y + 4y^3]
      Real expected_grad_x = 4.0 * px*px*px + 2.0 * px * py*py;
      Real expected_grad_y = 2.0 * px*px * py + 4.0 * py*py*py;
      
      EXPECT_NEAR(grad_value(0), expected_grad_x, 1e-9);
      EXPECT_NEAR(grad_value(1), expected_grad_y, 1e-9);
    }
  }

  TEST(Rodin_Variational_H1_Grad, ProjectGradOntoGridFunction_LinearFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // Create scalar H1<2> space for the original function
    H1 fes_scalar(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes_scalar);

    // Project linear function f(x,y) = 2x + 3y with gradient [2, 3]
    RealFunction linear_func([](const Geometry::Point& p) { 
      return 2.0 * p.x() + 3.0 * p.y(); 
    });
    gf.project(linear_func);

    // Create vector H1<2> space for the gradient (dimension = 2)
    H1 fes_vector(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf_grad(fes_vector);

    // Project Grad(gf) onto gf_grad
    gf_grad.project(Grad(gf));

    // Test at multiple random points
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

      // Evaluate the projected gradient grid function
      auto grad_value = gf_grad.getValue(p);
      
      // Should match the constant gradient [2, 3]
      EXPECT_NEAR(grad_value(0), 2.0, 0.1);  // Relaxed tolerance for FE projection
      EXPECT_NEAR(grad_value(1), 3.0, 0.1);
    }
  }

  TEST(Rodin_Variational_H1_Grad, ProjectGradOntoGridFunction_QuadraticFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // Create H1<3> space for the quadratic function
    H1 fes_scalar(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes_scalar);

    // Project quadratic function f(x,y) = x^2 + y^2 with gradient [2x, 2y]
    RealFunction quadratic_func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y(); 
    });
    gf.project(quadratic_func);

    // Create vector H1<3> space for the gradient
    H1 fes_vector(std::integral_constant<size_t, 3>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf_grad(fes_vector);

    // Project Grad(gf) onto gf_grad
    gf_grad = Grad(gf);  // Using operator= instead of project()

    // Test at multiple random points
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

      // Evaluate the projected gradient grid function
      auto grad_value = gf_grad.getValue(p);
      const auto& phys_coords = p.getPhysicalCoordinates();
      
      // Should match gradient [2x, 2y]
      EXPECT_NEAR(grad_value(0), 2.0 * phys_coords(0), 1.0);  // Larger tolerance for quadratic FE projection
      EXPECT_NEAR(grad_value(1), 2.0 * phys_coords(1), 1.0);
    }
  }

  // ============================================================================
  // Quadrilateral tests
  // ============================================================================

  TEST(Rodin_Variational_H1_Grad_Quad, ShapeFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
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

  TEST(Rodin_Variational_H1_Grad_Quad, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    auto grad_gf = Grad(gf);

    // Gradient of a scalar GridFunction should be a vector function
    EXPECT_EQ(&grad_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Grad_Quad, GridFunction_LinearFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Project a linear function: f(x,y) = x + y
    RealFunction linear_func([](const Geometry::Point& p) { return p.x() + p.y(); });
    gf.project(linear_func);

    auto grad_gf = Grad(gf);

    // For linear functions, gradient should be constant
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.4}};
    Point p(polytope, rc);

    auto grad_value = grad_gf.getValue(p);
    // Gradient should be approximately [1, 1]
    EXPECT_NEAR(grad_value(0), 1.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(grad_value(1), 1.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Grad_Quad, GridFunction_QuadraticFunction_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
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

  TEST(Rodin_Variational_H1_Grad_Quad, RandomCoordinates_QuadraticFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> can represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Project a quadratic function: f(x,y) = x^2 - xy + 2y^2
    // Gradient: [2x - y, -x + 4y]
    RealFunction quadratic_func([](const Geometry::Point& p) { 
      return p.x() * p.x() - p.x() * p.y() + 2.0 * p.y() * p.y(); 
    });
    gf.project(quadratic_func);

    auto grad_gf = Grad(gf);

    // Test at 15 random points across different cells
    RandomFloat gen(0.0, 1.0);
    for (int test = 0; test < 15; test++)
    {
      Index cellIdx = gen() * (mesh.getCellCount() - 1);
      auto it = mesh.getPolytope(mesh.getDimension(), cellIdx);
      const auto& polytope = *it;
      
      Real x = gen();
      Real y = gen();
      const Math::Vector<Real> rc{{x, y}};
      Point p(polytope, rc);

      auto grad_value = grad_gf.getValue(p);
      const auto& phys_coords = p.getPhysicalCoordinates();
      
      Real px = phys_coords(0);
      Real py = phys_coords(1);
      
      // Expected gradient: [2x - y, -x + 4y]
      EXPECT_NEAR(grad_value(0), 2.0 * px - py, 1e-10);
      EXPECT_NEAR(grad_value(1), -px + 4.0 * py, 1e-10);
    }
  }

  // Tests for setPoint and getBasis methods
  TEST(Rodin_Variational_H1_Grad, ShapeFunction_setPoint_Triangle)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    auto grad_u = Grad(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Create evaluation point
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);

    // Test setPoint
    grad_u.setPoint(p);
    
    // Verify getPoint returns the correct point
    const auto& retrieved_point = grad_u.getPoint();
    EXPECT_EQ(&retrieved_point, &p);
    EXPECT_EQ(&retrieved_point.getPolytope(), &polytope);
  }

  TEST(Rodin_Variational_H1_Grad, ShapeFunction_getBasis_Triangle_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    auto grad_u = Grad(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Create evaluation point at centroid
    const Math::Vector<Real> rc{{1.0/3.0, 1.0/3.0}};
    Point p(polytope, rc);

    // Set the point
    grad_u.setPoint(p);
    
    // Test getBasis for each local DOF
    const size_t num_dofs = grad_u.getDOFs(polytope);
    EXPECT_EQ(num_dofs, 6); // H1<2> on triangle has 6 DOFs
    
    // Verify that getBasis returns vector values for each local DOF
    for (size_t local = 0; local < num_dofs; local++)
    {
      auto basis_grad = grad_u.getBasis(local);
      // Should be a 2D vector for 2D mesh
      EXPECT_EQ(basis_grad.size(), 2);
      // Values should be finite
      EXPECT_TRUE(std::isfinite(basis_grad(0)));
      EXPECT_TRUE(std::isfinite(basis_grad(1)));
    }
  }

  TEST(Rodin_Variational_H1_Grad, ShapeFunction_getBasis_Tetrahedron_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    auto grad_u = Grad(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Create evaluation point
    const Math::Vector<Real> rc{{0.2, 0.2, 0.2}};
    Point p(polytope, rc);

    // Set the point
    grad_u.setPoint(p);
    
    // Test getBasis for each local DOF
    const size_t num_dofs = grad_u.getDOFs(polytope);
    EXPECT_EQ(num_dofs, 10); // H1<2> on tetrahedron has 10 DOFs
    
    // Verify that getBasis returns vector values for each local DOF
    for (size_t local = 0; local < num_dofs; local++)
    {
      auto basis_grad = grad_u.getBasis(local);
      // Should be a 3D vector for 3D mesh
      EXPECT_EQ(basis_grad.size(), 3);
      // Values should be finite
      EXPECT_TRUE(std::isfinite(basis_grad(0)));
      EXPECT_TRUE(std::isfinite(basis_grad(1)));
      EXPECT_TRUE(std::isfinite(basis_grad(2)));
    }
  }

  TEST(Rodin_Variational_H1_Grad, ShapeFunction_getBasis_Quadrilateral_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    auto grad_u = Grad(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Create evaluation point
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    // Set the point
    grad_u.setPoint(p);
    
    // Test getBasis for each local DOF
    const size_t num_dofs = grad_u.getDOFs(polytope);
    EXPECT_EQ(num_dofs, 9); // H1<2> on quadrilateral has 9 DOFs
    
    // Verify that getBasis returns vector values for each local DOF
    for (size_t local = 0; local < num_dofs; local++)
    {
      auto basis_grad = grad_u.getBasis(local);
      // Should be a 2D vector for 2D mesh
      EXPECT_EQ(basis_grad.size(), 2);
      // Values should be finite
      EXPECT_TRUE(std::isfinite(basis_grad(0)));
      EXPECT_TRUE(std::isfinite(basis_grad(1)));
    }
  }

  TEST(Rodin_Variational_H1_Grad, ShapeFunction_setPoint_MultipleCalls)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    auto grad_u = Grad(u);

    // Get first element
    auto it1 = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope1 = *it1;
    const Math::Vector<Real> rc1{{0.2, 0.3}};
    Point p1(polytope1, rc1);

    // Set first point
    grad_u.setPoint(p1);
    EXPECT_EQ(&grad_u.getPoint(), &p1);
    
    // Get gradients at first point
    auto basis_grad1_0 = grad_u.getBasis(0);

    // Get second element
    auto it2 = mesh.getPolytope(mesh.getDimension(), 1);
    const auto& polytope2 = *it2;
    const Math::Vector<Real> rc2{{0.4, 0.5}};
    Point p2(polytope2, rc2);

    // Set second point
    grad_u.setPoint(p2);
    EXPECT_EQ(&grad_u.getPoint(), &p2);
    
    // Get gradients at second point (should be different if elements differ)
    auto basis_grad2_0 = grad_u.getBasis(0);
    
    // The gradients should potentially be different for different elements
    // Just verify both are valid
    EXPECT_TRUE(std::isfinite(basis_grad1_0(0)));
    EXPECT_TRUE(std::isfinite(basis_grad2_0(0)));
  }

  TEST(Rodin_Variational_H1_Grad, ShapeFunction_getBasis_PartitionOfUnity_Triangle)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Use H1<1> for simpler partition of unity test
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TrialFunction u(fes);
    auto grad_u = Grad(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Test at centroid
    const Math::Vector<Real> rc{{1.0/3.0, 1.0/3.0}};
    Point p(polytope, rc);

    grad_u.setPoint(p);
    
    const size_t num_dofs = grad_u.getDOFs(polytope);
    EXPECT_EQ(num_dofs, 3); // H1<1> on triangle has 3 DOFs (vertices)
    
    // Sum of gradients should equal zero (constant function gradient is zero)
    Math::Vector<Real> grad_sum(2);
    grad_sum.setZero();
    
    for (size_t local = 0; local < num_dofs; local++)
    {
      auto basis_grad = grad_u.getBasis(local);
      grad_sum += basis_grad;
    }
    
    // Sum of basis function gradients should be zero
    EXPECT_NEAR(grad_sum.norm(), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_H1_Grad, ShapeFunction_getBasis_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    TrialFunction u(fes);
    auto grad_u = Grad(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    const Math::Vector<Real> rc{{0.3, 0.4}};
    Point p(polytope, rc);

    grad_u.setPoint(p);
    
    const size_t num_dofs = grad_u.getDOFs(polytope);
    EXPECT_EQ(num_dofs, 10); // H1<3> on triangle has 10 DOFs
    
    // Verify all basis gradients are valid
    for (size_t local = 0; local < num_dofs; local++)
    {
      auto basis_grad = grad_u.getBasis(local);
      EXPECT_EQ(basis_grad.size(), 2);
      EXPECT_TRUE(std::isfinite(basis_grad(0)));
      EXPECT_TRUE(std::isfinite(basis_grad(1)));
    }
  }

  TEST(Rodin_Variational_H1_Grad, ShapeFunction_getBasis_RandomPoints)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    auto grad_u = Grad(u);

    // Test with 10 random evaluation points
    const size_t num_random_points = 10;
    const size_t num_cells = mesh.getCellCount();
    
    for (size_t i = 0; i < num_random_points; i++)
    {
      // Random cell
      const size_t cell_idx = randomSize(0, num_cells - 1);
      auto it = mesh.getPolytope(mesh.getDimension(), cell_idx);
      const auto& polytope = *it;
      
      // Random reference coordinates (barycentric for triangle)
      const Real r1 = randomReal(0.0, 1.0);
      const Real r2 = randomReal(0.0, 1.0 - r1);
      const Math::Vector<Real> rc{{r1, r2}};
      Point p(polytope, rc);

      // Set point and verify getBasis works
      grad_u.setPoint(p);
      EXPECT_EQ(&grad_u.getPoint(), &p);
      
      const size_t num_dofs = grad_u.getDOFs(polytope);
      for (size_t local = 0; local < num_dofs; local++)
      {
        auto basis_grad = grad_u.getBasis(local);
        EXPECT_EQ(basis_grad.size(), 2);
        EXPECT_TRUE(std::isfinite(basis_grad(0)));
        EXPECT_TRUE(std::isfinite(basis_grad(1)));
        
        // Gradient values should be reasonable (not too large)
        EXPECT_LT(std::abs(basis_grad(0)), 1000.0);
        EXPECT_LT(std::abs(basis_grad(1)), 1000.0);
      }
    }
  }
}
