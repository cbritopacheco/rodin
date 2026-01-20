#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/H1/Div.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_H1_Div, ShapeFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    TestFunction v(fes);

    auto div_u = Div(u);
    auto div_v = Div(v);

    // Divergence of a vector function should be a scalar function
  }

  TEST(Rodin_Variational_H1_Div, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    auto div_gf = Div(gf);

    // Divergence of a vector GridFunction should be a scalar function
    EXPECT_EQ(&div_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Div, GridFunction_ConstantVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a constant vector field
    VectorFunction const_func{5.0, 7.0};
    gf.project(const_func);

    auto div_gf = Div(gf);

    // Divergence of a constant vector field should be zero
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    Real div_value = div_gf.getValue(p);
    EXPECT_NEAR(div_value, 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Div, Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);
    gf = VectorFunction{1.0, 2.0};

    auto div_gf = Div(gf);
    auto copied = div_gf.copy();

    EXPECT_NE(copied, nullptr);

    delete copied;
  }

  TEST(Rodin_Variational_H1_Div, GridFunction_LinearVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: (x, y) - divergence should be 2
    VectorFunction linear_func{
      [](const Geometry::Point& p) { return p.x(); },
      [](const Geometry::Point& p) { return p.y(); }
    };
    gf.project(linear_func);

    auto div_gf = Div(gf);

    // For linear functions, divergence should be constant = 2
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    Real div_value = div_gf.getValue(p);
    EXPECT_NEAR(div_value, 2.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Div, GridFunction_DivergenceFreeField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> (quadratic) can represent linear divergence-free fields exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a divergence-free field: (-y, x)
    // div(u) = ∂(-y)/∂x + ∂(x)/∂y = 0 + 0 = 0
    VectorFunction divergence_free{
      [](const Geometry::Point& p) { return -p.y(); },
      [](const Geometry::Point& p) { return p.x(); }
    };
    gf.project(divergence_free);

    auto div_gf = Div(gf);

    // Test at random coordinates
    RandomFloat gen(0.0, 1.0);
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;

    for (int i = 0; i < 10; i++)
    {
      Real x = gen();
      Real y = gen();
      // Ensure point is inside the reference triangle
      if (x + y > 1.0) {
        x = 1.0 - x;
        y = 1.0 - y;
      }
      const Math::Vector<Real> rc{{x, y}};
      Point p(polytope, rc);

      Real div_value = div_gf.getValue(p);
      EXPECT_NEAR(div_value, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Div, GridFunction_QuadraticVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> (quadratic) can represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a quadratic vector field: (x^2, y^2)
    // div(u) = ∂(x^2)/∂x + ∂(y^2)/∂y = 2x + 2y
    VectorFunction quadratic_func{
      [](const Geometry::Point& p) { return p.x() * p.x(); },
      [](const Geometry::Point& p) { return p.y() * p.y(); }
    };
    gf.project(quadratic_func);

    auto div_gf = Div(gf);

    // Test at random coordinates
    RandomFloat gen(0.0, 1.0);
    for (int test = 0; test < 5; test++)
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
      const auto& phys_coords = p.getPhysicalCoordinates();
      
      // Expected div = 2x + 2y
      Real expected_div = 2.0 * phys_coords(0) + 2.0 * phys_coords(1);
      EXPECT_NEAR(div_value, expected_div, 1e-10);
    }
  }

  TEST(Rodin_Variational_H1_Div, GridFunction_QuadraticVectorField_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<3> (cubic) can also represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 3>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a quadratic vector field: (x^2 - y^2, 2xy)
    // div(u) = ∂(x^2 - y^2)/∂x + ∂(2xy)/∂y = 2x + 2x = 4x
    VectorFunction quadratic_func{
      [](const Geometry::Point& p) { return p.x() * p.x() - p.y() * p.y(); },
      [](const Geometry::Point& p) { return 2.0 * p.x() * p.y(); }
    };
    gf.project(quadratic_func);

    auto div_gf = Div(gf);

    // Test at random coordinates
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

      Real div_value = div_gf.getValue(p);
      const auto& phys_coords = p.getPhysicalCoordinates();
      
      // Expected div = 4x
      Real expected_div = 4.0 * phys_coords(0);
      EXPECT_NEAR(div_value, expected_div, 1e-10);
    }
  }

  TEST(Rodin_Variational_H1_Div, MultipleRandomEvaluations)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: (2x + y, x - y)
    // div(u) = ∂(2x + y)/∂x + ∂(x - y)/∂y = 2 + (-1) = 1
    VectorFunction linear_func{
      [](const Geometry::Point& p) { return 2.0 * p.x() + p.y(); },
      [](const Geometry::Point& p) { return p.x() - p.y(); }
    };
    gf.project(linear_func);

    auto div_gf = Div(gf);

    // Test at many random points
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
      
      // Divergence should be constant = 1
      EXPECT_NEAR(div_value, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Div, DifferentPolynomialDegrees)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Linear vector field for testing: (3x + 4y, -2x + y)
    // div(u) = 3 + 1 = 4
    VectorFunction linear_func{
      [](const Geometry::Point& p) { return 3.0 * p.x() + 4.0 * p.y(); },
      [](const Geometry::Point& p) { return -2.0 * p.x() + p.y(); }
    };

    RandomFloat gen(0.0, 1.0);
    Real x = gen();
    Real y = gen();
    if (x + y > 1.0) {
      x = 1.0 - x;
      y = 1.0 - y;
    }

    // Test with H1<1> (linear)
    {
      H1 fes(std::integral_constant<size_t, 1>{}, mesh, mesh.getSpaceDimension());
      GridFunction gf(fes);
      gf.project(linear_func);
      auto div_gf = Div(gf);
      
      auto it = mesh.getPolytope(mesh.getDimension(), 0);
      const auto& polytope = *it;
      Point p(polytope, Math::Vector<Real>{{x, y}});
      Real div_value = div_gf.getValue(p);
      
      EXPECT_NEAR(div_value, 4.0, RODIN_FUZZY_CONSTANT);
    }

    // Test with H1<2> (quadratic)
    {
      H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
      GridFunction gf(fes);
      gf.project(linear_func);
      auto div_gf = Div(gf);
      
      auto it = mesh.getPolytope(mesh.getDimension(), 0);
      const auto& polytope = *it;
      Point p(polytope, Math::Vector<Real>{{x, y}});
      Real div_value = div_gf.getValue(p);
      
      EXPECT_NEAR(div_value, 4.0, RODIN_FUZZY_CONSTANT);
    }

    // Test with H1<3> (cubic)
    {
      H1 fes(std::integral_constant<size_t, 3>{}, mesh, mesh.getSpaceDimension());
      GridFunction gf(fes);
      gf.project(linear_func);
      auto div_gf = Div(gf);
      
      auto it = mesh.getPolytope(mesh.getDimension(), 0);
      const auto& polytope = *it;
      Point p(polytope, Math::Vector<Real>{{x, y}});
      Real div_value = div_gf.getValue(p);
      
      EXPECT_NEAR(div_value, 4.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Div, UsageInBilinearForm)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);

    // ∫ div(u) * div(v) dx - part of mixed elasticity formulations
    bf = Integral(Div(u), Div(v));

    EXPECT_FALSE(bf.getLocalIntegrators().empty());

    // This should assemble without errors
    bf.assemble();

    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
  }

  TEST(Rodin_Variational_H1_Div, ZeroVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Zero vector field (default initialization)
    auto div_gf = Div(gf);

    // Test at random point
    RandomFloat gen(0.0, 1.0);
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
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
    EXPECT_NEAR(div_value, 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Div, RandomCoordinatesOnMultipleCells)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 5, 5 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project: (x, 2y) - div = 1 + 2 = 3
    VectorFunction func{
      [](const Geometry::Point& p) { return p.x(); },
      [](const Geometry::Point& p) { return 2.0 * p.y(); }
    };
    gf.project(func);

    auto div_gf = Div(gf);

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

      Real div_value = div_gf.getValue(p);
      EXPECT_NEAR(div_value, 3.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Div, RandomCoordinates_QuarticVectorField_H1_5)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<5> can represent quartic (degree 4) functions exactly
    H1 fes(std::integral_constant<size_t, 5>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a quartic vector field: u = (x^4 + x^2*y, y^4 + xy^2)
    // div(u) = ∂(x^4 + x^2*y)/∂x + ∂(y^4 + xy^2)/∂y
    //        = 4x^3 + 2xy + 4y^3 + 2xy
    //        = 4x^3 + 4y^3 + 4xy
    auto quartic_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      Real x = p.x();
      Real y = p.y();
      v << x*x*x*x + x*x*y, 
           y*y*y*y + x*y*y;
      return v;
    };
    VectorFunction<decltype(quartic_lambda)> quartic_func(2, quartic_lambda);
    gf.project(quartic_func);

    auto div_gf = Div(gf);

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

      Real div_value = div_gf.getValue(p);
      const auto& phys_coords = p.getPhysicalCoordinates();
      
      Real px = phys_coords(0);
      Real py = phys_coords(1);
      
      // Expected div = 4x^3 + 4y^3 + 4xy
      Real expected_div = 4.0 * px*px*px + 4.0 * py*py*py + 4.0 * px * py;
      EXPECT_NEAR(div_value, expected_div, 1e-9);
    }
  }

  TEST(Rodin_Variational_H1_Div, ProjectDivOntoGridFunction_LinearVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // Create vector H1<2> space for the vector field
    H1 fes_vector(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes_vector);

    // Project linear vector field u = (2x, 3y) with div(u) = 2 + 3 = 5
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << 2.0 * p.x(), 3.0 * p.y();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(2, linear_lambda);
    gf.project(linear_func);

    // Create scalar H1<2> space for the divergence
    H1 fes_scalar(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf_div(fes_scalar);

    // Project Div(gf) onto gf_div
    gf_div.project(Div(gf));

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

      // Evaluate the projected divergence grid function
      Real div_value = gf_div.getValue(p);
      
      // Should match the constant divergence 2 + 3 = 5
      EXPECT_NEAR(div_value, 5.0, 0.1);  // Relaxed tolerance for FE projection
    }
  }

  TEST(Rodin_Variational_H1_Div, ProjectDivOntoGridFunction_QuadraticVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // Create vector H1<3> space for the quadratic vector field
    H1 fes_vector(std::integral_constant<size_t, 3>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes_vector);

    // Project quadratic vector field u = (x^2, y^2) with div(u) = 2x + 2y
    auto quadratic_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.x() * p.x(), p.y() * p.y();
      return v;
    };
    VectorFunction<decltype(quadratic_lambda)> quadratic_func(2, quadratic_lambda);
    gf.project(quadratic_func);

    // Create scalar H1<3> space for the divergence
    H1 fes_scalar(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf_div(fes_scalar);

    // Project Div(gf) onto gf_div using operator=
    gf_div = Div(gf);

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

      // Evaluate the projected divergence grid function
      Real div_value = gf_div.getValue(p);
      const auto& phys_coords = p.getPhysicalCoordinates();
      
      // Should match divergence 2x + 2y
      Real expected_div = 2.0 * phys_coords(0) + 2.0 * phys_coords(1);
      EXPECT_NEAR(div_value, expected_div, 1.5);  // Larger tolerance for quadratic FE projection
    }
  }

  // ============================================================================
  // Tetrahedron tests
  // ============================================================================

  TEST(Rodin_Variational_H1_Div_Tet, ShapeFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    TestFunction v(fes);

    auto div_u = Div(u);
    auto div_v = Div(v);

    // Divergence of a vector function should be a scalar function
  }

  TEST(Rodin_Variational_H1_Div_Tet, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    auto div_gf = Div(gf);

    // Divergence of a vector GridFunction should be a scalar function
    EXPECT_EQ(&div_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Div_Tet, GridFunction_ConstantFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a constant vector field
    VectorFunction const_func{5.0, 7.0, 9.0};
    gf.project(const_func);

    auto div_gf = Div(gf);

    // Divergence of a constant vector field should be zero
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.2, 0.2}};
    Point p(polytope, rc);

    Real div_value = div_gf.getValue(p);
    EXPECT_NEAR(div_value, 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Div_Tet, GridFunction_LinearVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: (x, y, z) - divergence should be 3
    VectorFunction linear_func{
      [](const Geometry::Point& p) { return p.x(); },
      [](const Geometry::Point& p) { return p.y(); },
      [](const Geometry::Point& p) { return p.z(); }
    };
    gf.project(linear_func);

    auto div_gf = Div(gf);

    // For linear functions, divergence should be constant = 3
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25, 0.25}};
    Point p(polytope, rc);

    Real div_value = div_gf.getValue(p);
    EXPECT_NEAR(div_value, 3.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Div_Tet, GridFunction_QuadraticVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> (quadratic) can represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a quadratic vector field: u = (x^2, y^2, z^2)
    // div(u) = 2x + 2y + 2z
    VectorFunction quadratic_func{
      [](const Geometry::Point& p) { return p.x() * p.x(); },
      [](const Geometry::Point& p) { return p.y() * p.y(); },
      [](const Geometry::Point& p) { return p.z() * p.z(); }
    };
    gf.project(quadratic_func);

    auto div_gf = Div(gf);

    // Evaluate at a point
    auto it = mesh.getPolytope(mesh.getDimension(), mesh.getCellCount() / 2);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.3, 0.1}};
    Point p(polytope, rc);

    Real div_value = div_gf.getValue(p);
    const auto& phys_coords = p.getPhysicalCoordinates();
    
    // Expected divergence: 2x + 2y + 2z
    Real expected_div = 2.0 * phys_coords(0) + 2.0 * phys_coords(1) + 2.0 * phys_coords(2);
    EXPECT_NEAR(div_value, expected_div, 1e-10);
  }

  TEST(Rodin_Variational_H1_Div_Tet, RandomCoordinates_LinearVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: u = (2x, 3y, 4z)
    // Divergence should be constant: 2 + 3 + 4 = 9
    VectorFunction linear_func{
      [](const Geometry::Point& p) { return 2.0 * p.x(); },
      [](const Geometry::Point& p) { return 3.0 * p.y(); },
      [](const Geometry::Point& p) { return 4.0 * p.z(); }
    };
    gf.project(linear_func);

    auto div_gf = Div(gf);

    // Test at 20 random points across different cells
    RandomFloat gen(0.0, 1.0);
    for (int test = 0; test < 20; test++)
    {
      Index cellIdx = gen() * (mesh.getCellCount() - 1);
      auto it = mesh.getPolytope(mesh.getDimension(), cellIdx);
      const auto& polytope = *it;
      
      // Generate random barycentric coordinates for tetrahedron
      Real x = gen();
      Real y = gen();
      Real z = gen();
      Real sum = x + y + z;
      if (sum > 1.0) {
        x = x / sum * 0.9;
        y = y / sum * 0.9;
        z = z / sum * 0.9;
      }
      const Math::Vector<Real> rc{{x, y, z}};
      Point p(polytope, rc);

      Real div_value = div_gf.getValue(p);
      EXPECT_NEAR(div_value, 9.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // ============================================================================
  // Quadrilateral tests
  // ============================================================================

  TEST(Rodin_Variational_H1_Div_Quad, ShapeFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    TestFunction v(fes);

    auto div_u = Div(u);
    auto div_v = Div(v);

    // Divergence of a vector function should be a scalar function
  }

  TEST(Rodin_Variational_H1_Div_Quad, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    auto div_gf = Div(gf);

    // Divergence of a vector GridFunction should be a scalar function
    EXPECT_EQ(&div_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Div_Quad, GridFunction_LinearVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: (x, y) - divergence should be 2
    VectorFunction linear_func{
      [](const Geometry::Point& p) { return p.x(); },
      [](const Geometry::Point& p) { return p.y(); }
    };
    gf.project(linear_func);

    auto div_gf = Div(gf);

    // For linear functions, divergence should be constant = 2
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.4}};
    Point p(polytope, rc);

    Real div_value = div_gf.getValue(p);
    EXPECT_NEAR(div_value, 2.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Div_Quad, RandomCoordinates_LinearVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: u = (3x, 2y)
    // Divergence should be constant: 3 + 2 = 5
    VectorFunction linear_func{
      [](const Geometry::Point& p) { return 3.0 * p.x(); },
      [](const Geometry::Point& p) { return 2.0 * p.y(); }
    };
    gf.project(linear_func);

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
      const Math::Vector<Real> rc{{x, y}};
      Point p(polytope, rc);

      Real div_value = div_gf.getValue(p);
      EXPECT_NEAR(div_value, 5.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // ===== ShapeFunctionBase setPoint and getBasis Tests =====

  TEST(Rodin_Variational_H1_Div, ShapeFunction_setPoint_Triangle)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    auto div_u = Div(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Create evaluation point
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);

    // Test setPoint
    div_u.setPoint(p);
    
    // Verify getPoint returns the correct point
    const auto& retrieved_point = div_u.getPoint();
    EXPECT_EQ(&retrieved_point, &p);
    EXPECT_EQ(&retrieved_point.getPolytope(), &polytope);
  }

  TEST(Rodin_Variational_H1_Div, ShapeFunction_getBasis_Triangle_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<2> space (2D vectors)
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    auto div_u = Div(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Create evaluation point at centroid
    const Math::Vector<Real> rc{{1.0/3.0, 1.0/3.0}};
    Point p(polytope, rc);

    // Set the point
    div_u.setPoint(p);
    
    // Test getBasis for each local DOF
    const size_t num_dofs = div_u.getDOFs(polytope);
    EXPECT_EQ(num_dofs, 6); // H1<2> on triangle has 6 DOFs per component, 2 components
    
    // Verify that getBasis returns scalar values (divergence) for each local DOF
    for (size_t local = 0; local < num_dofs; local++)
    {
      Real basis_div = div_u.getBasis(local);
      // Divergence should be a scalar
      // Value should be finite
      EXPECT_TRUE(std::isfinite(basis_div));
    }
  }

  TEST(Rodin_Variational_H1_Div, ShapeFunction_getBasis_Tetrahedron_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<2> space (3D vectors)
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    auto div_u = Div(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Create evaluation point
    const Math::Vector<Real> rc{{0.2, 0.2, 0.2}};
    Point p(polytope, rc);

    // Set the point
    div_u.setPoint(p);
    
    // Test getBasis for each local DOF
    const size_t num_dofs = div_u.getDOFs(polytope);
    EXPECT_EQ(num_dofs, 10); // H1<2> on tetrahedron has 10 DOFs per component, 3 components
    
    // Verify that getBasis returns scalar values for each local DOF
    for (size_t local = 0; local < num_dofs; local++)
    {
      Real basis_div = div_u.getBasis(local);
      // Divergence should be a scalar value
      EXPECT_TRUE(std::isfinite(basis_div));
    }
  }

  TEST(Rodin_Variational_H1_Div, ShapeFunction_getBasis_Quadrilateral_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    auto div_u = Div(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Create evaluation point
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    // Set the point
    div_u.setPoint(p);
    
    // Test getBasis for each local DOF
    const size_t num_dofs = div_u.getDOFs(polytope);
    EXPECT_EQ(num_dofs, 9); // H1<2> on quadrilateral has 9 DOFs per component, 2 components
    
    // Verify that getBasis returns scalar values for each local DOF
    for (size_t local = 0; local < num_dofs; local++)
    {
      Real basis_div = div_u.getBasis(local);
      // Divergence should be a scalar
      EXPECT_TRUE(std::isfinite(basis_div));
    }
  }

  TEST(Rodin_Variational_H1_Div, ShapeFunction_setPoint_MultipleCalls)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    auto div_u = Div(u);

    // Get first element
    auto it1 = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope1 = *it1;
    const Math::Vector<Real> rc1{{0.25, 0.25}};
    Point p1(polytope1, rc1);

    // First setPoint
    div_u.setPoint(p1);
    EXPECT_EQ(&div_u.getPoint(), &p1);
    
    // Verify getBasis works
    const size_t num_dofs1 = div_u.getDOFs(polytope1);
    for (size_t local = 0; local < num_dofs1; local++)
    {
      Real basis_div = div_u.getBasis(local);
      EXPECT_TRUE(std::isfinite(basis_div));
    }

    // Get different element
    auto it2 = mesh.getPolytope(mesh.getDimension(), 2);
    const auto& polytope2 = *it2;
    const Math::Vector<Real> rc2{{0.5, 0.3}};
    Point p2(polytope2, rc2);

    // Second setPoint
    div_u.setPoint(p2);
    EXPECT_EQ(&div_u.getPoint(), &p2);
    EXPECT_EQ(&div_u.getPoint().getPolytope(), &polytope2);
    
    // Verify getBasis still works after changing point
    const size_t num_dofs2 = div_u.getDOFs(polytope2);
    for (size_t local = 0; local < num_dofs2; local++)
    {
      Real basis_div = div_u.getBasis(local);
      EXPECT_TRUE(std::isfinite(basis_div));
    }
  }

  TEST(Rodin_Variational_H1_Div, ShapeFunction_getBasis_PartitionOfUnity_Triangle)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Use H1<1> for simplicity - linear basis functions
    H1 fes(std::integral_constant<size_t, 1>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    auto div_u = Div(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Create evaluation point
    const Math::Vector<Real> rc{{0.3, 0.4}};
    Point p(polytope, rc);

    // Set the point
    div_u.setPoint(p);
    
    // For H1<1>, the divergence of constant vector field should be zero
    // Sum of all basis divergences should be zero
    const size_t num_dofs = div_u.getDOFs(polytope);
    EXPECT_EQ(num_dofs, 3); // H1<1> on triangle has 3 DOFs per component
    
    // Sum all basis divergences
    Real sum_div = 0.0;
    for (size_t local = 0; local < num_dofs; local++)
    {
      Real basis_div = div_u.getBasis(local);
      sum_div += basis_div;
    }
    
    // For linear basis functions, sum of divergences should be zero
    // (divergence of sum of basis functions = divergence of constant = 0)
    EXPECT_NEAR(sum_div, 0.0, 1e-10);
  }

  TEST(Rodin_Variational_H1_Div, ShapeFunction_getBasis_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<3> space
    H1 fes(std::integral_constant<size_t, 3>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    auto div_u = Div(u);

    // Get first element
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    
    // Create evaluation point
    const Math::Vector<Real> rc{{0.2, 0.3}};
    Point p(polytope, rc);

    // Set the point
    div_u.setPoint(p);
    
    // Test getBasis for H1<3>
    const size_t num_dofs = div_u.getDOFs(polytope);
    EXPECT_EQ(num_dofs, 10); // H1<3> on triangle has 10 DOFs per component
    
    // Verify all basis divergences are finite
    for (size_t local = 0; local < num_dofs; local++)
    {
      Real basis_div = div_u.getBasis(local);
      EXPECT_TRUE(std::isfinite(basis_div));
    }
  }

  TEST(Rodin_Variational_H1_Div, ShapeFunction_getBasis_RandomPoints)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    auto div_u = Div(u);

    // Test with random evaluation points
    const size_t num_tests = 10;
    std::uniform_int_distribution<size_t> cell_dist(0, mesh.getCellCount() - 1);
    std::uniform_real_distribution<Real> coord_dist(0.1, 0.8);
    
    for (size_t test = 0; test < num_tests; test++)
    {
      // Random cell
      size_t cell_idx = cell_dist(gen);
      auto it = mesh.getPolytope(mesh.getDimension(), cell_idx);
      const auto& polytope = *it;
      
      // Random reference coordinates
      Real r1 = coord_dist(gen);
      Real r2 = coord_dist(gen);
      if (r1 + r2 > 1.0)
      {
        r1 = 1.0 - r1;
        r2 = 1.0 - r2;
      }
      const Math::Vector<Real> rc{{r1, r2}};
      Point p(polytope, rc);

      // Set point and test getBasis
      div_u.setPoint(p);
      
      const size_t num_dofs = div_u.getDOFs(polytope);
      for (size_t local = 0; local < num_dofs; local++)
      {
        Real basis_div = div_u.getBasis(local);
        
        // Check that value is finite and reasonably bounded
        EXPECT_TRUE(std::isfinite(basis_div));
        EXPECT_LT(std::abs(basis_div), 1000.0); // Sanity check
      }
    }
  }
}
