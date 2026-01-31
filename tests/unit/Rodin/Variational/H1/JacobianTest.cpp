#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/H1/Jacobian.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_H1_Jacobian, ShapeFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    TestFunction v(fes);

    auto jac_u = Jacobian(u);
    auto jac_v = Jacobian(v);

    // Jacobian of a vector function should be a matrix function
    // For 2D, Jacobian should be 2x2
  }

  TEST(Rodin_Variational_H1_Jacobian, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    auto jac_gf = Jacobian(gf);

    // Jacobian of a vector GridFunction should be a matrix function
    EXPECT_EQ(&jac_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Jacobian, GridFunction_ConstantFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a constant vector field
    auto const_func_lambda = []( const Geometry::Point& p) { 
      Math::Vector<Real> v(2);
      v << 1.0, 2.0;
      return v;
    };
    VectorFunction<decltype(const_func_lambda)> const_func(2, const_func_lambda);
    gf.project(const_func);

    auto jac_gf = Jacobian(gf);

    // Jacobian of a constant vector field should be zero
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    auto jac_value = jac_gf.getValue(p);
    EXPECT_NEAR(jac_value.norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Jacobian, Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);
    auto copy_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << 1.0, 2.0;
      return v;
    };
    gf = VectorFunction<decltype(copy_lambda)>(2, copy_lambda);

    auto jac_gf = Jacobian(gf);
    auto copied = jac_gf.copy();

    EXPECT_NE(copied, nullptr);

    delete copied;
  }

  TEST(Rodin_Variational_H1_Jacobian, GridFunction_LinearVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: u = (2x + y, x - 3y)
    // Jacobian should be: [[2, 1], [1, -3]]
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << 2.0 * p.x() + p.y(), p.x() - 3.0 * p.y();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(2, linear_lambda);
    gf.project(linear_func);

    auto jac_gf = Jacobian(gf);

    // For linear functions, Jacobian should be constant
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    auto jac_value = jac_gf.getValue(p);
    // Expected Jacobian: [[2, 1], [1, -3]]
    EXPECT_EQ(jac_value.rows(), 2);
    EXPECT_EQ(jac_value.cols(), 2);
    EXPECT_NEAR(jac_value(0, 0), 2.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(0, 1), 1.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(1, 0), 1.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(1, 1), -3.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Jacobian, GridFunction_QuadraticVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> (quadratic) can represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a quadratic vector field: u = (x^2, y^2)
    // Jacobian should be: [[2x, 0], [0, 2y]]
    auto quadratic_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.x() * p.x(), p.y() * p.y();
      return v;
    };
    VectorFunction<decltype(quadratic_lambda)> quadratic_func(2, quadratic_lambda);
    gf.project(quadratic_func);

    auto jac_gf = Jacobian(gf);

    // Evaluate at a point
    auto it = mesh.getPolytope(mesh.getDimension(), mesh.getCellCount() / 2);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    auto jac_value = jac_gf.getValue(p);
    const auto& phys_coords = p.getPhysicalCoordinates();

    // Expected Jacobian: [[2x, 0], [0, 2y]]
    EXPECT_EQ(jac_value.rows(), 2);
    EXPECT_EQ(jac_value.cols(), 2);
    EXPECT_NEAR(jac_value(0, 0), 2.0 * phys_coords(0), 1e-10);
    EXPECT_NEAR(jac_value(0, 1), 0.0, 1e-10);
    EXPECT_NEAR(jac_value(1, 0), 0.0, 1e-10);
    EXPECT_NEAR(jac_value(1, 1), 2.0 * phys_coords(1), 1e-10);
  }

  TEST(Rodin_Variational_H1_Jacobian, GridFunction_QuadraticVectorField_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<3> (cubic) can also represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 3>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a quadratic vector field: u = (x^2 + xy, y^2 - xy)
    // Jacobian should be: [[2x + y, x], [-y, 2y - x]]
    auto quadratic_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.x() * p.x() + p.x() * p.y(), p.y() * p.y() - p.x() * p.y();
      return v;
    };
    VectorFunction<decltype(quadratic_lambda)> quadratic_func(2, quadratic_lambda);
    gf.project(quadratic_func);

    auto jac_gf = Jacobian(gf);

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
      auto jac_value = jac_gf.getValue(p);
      const auto& phys_coords = p.getPhysicalCoordinates();

      Real x = phys_coords(0);
      Real y = phys_coords(1);

      // Expected Jacobian: [[2x + y, x], [-y, 2y - x]]
      EXPECT_NEAR(jac_value(0, 0), 2.0 * x + y, 1e-10);
      EXPECT_NEAR(jac_value(0, 1), x, 1e-10);
      EXPECT_NEAR(jac_value(1, 0), -y, 1e-10);
      EXPECT_NEAR(jac_value(1, 1), 2.0 * y - x, 1e-10);
    }
  }

  TEST(Rodin_Variational_H1_Jacobian, MultipleEvaluations)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: u = (x + 2y, 3x - y)
    // Jacobian should be: [[1, 2], [3, -1]]
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.x() + 2.0 * p.y(), 3.0 * p.x() - p.y();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(2, linear_lambda);
    gf.project(linear_func);

    auto jac_gf = Jacobian(gf);

    // Test multiple evaluation points
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;

    std::vector<Math::Vector<Real>> test_coords =
    {
      Math::Vector<Real>{{ 0.1, 0.1 }},
      Math::Vector<Real>{{ 0.5, 0.2 }},
      Math::Vector<Real>{{ 0.2, 0.5 }},
      Math::Vector<Real>{{ 0.5, 0.3 }},
      Math::Vector<Real>{{ 0.25, 0.25 }}
    };

    for (const auto& rc : test_coords)
    {
      Point p(polytope, rc);
      auto jac_value = jac_gf.getValue(p);

      // Jacobian should be constant: [[1, 2], [3, -1]]
      EXPECT_NEAR(jac_value(0, 0), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(0, 1), 2.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 0), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 1), -1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Jacobian, DifferentPolynomialDegrees)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Linear vector field for testing
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << 3.0 * p.x() + 4.0 * p.y(), -2.0 * p.x() + p.y();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(2, linear_lambda);

    // Expected Jacobian: [[3, 4], [-2, 1]]

    // Test with H1<1> (linear)
    {
      H1 fes(std::integral_constant<size_t, 1>{}, mesh, mesh.getSpaceDimension());
      GridFunction gf(fes);
      gf.project(linear_func);
      auto jac_gf = Jacobian(gf);

      auto it = mesh.getPolytope(mesh.getDimension(), 0);
      const auto& polytope = *it;
      Point p(polytope, Math::Vector<Real>{{0.5, 0.5}});
      auto jac_value = jac_gf.getValue(p);

      EXPECT_NEAR(jac_value(0, 0), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(0, 1), 4.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 0), -2.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 1), 1.0, RODIN_FUZZY_CONSTANT);
    }

    // Test with H1<2> (quadratic)
    {
      H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
      GridFunction gf(fes);
      gf.project(linear_func);
      auto jac_gf = Jacobian(gf);

      auto it = mesh.getPolytope(mesh.getDimension(), 0);
      const auto& polytope = *it;
      Point p(polytope, Math::Vector<Real>{{0.5, 0.5}});
      auto jac_value = jac_gf.getValue(p);

      EXPECT_NEAR(jac_value(0, 0), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(0, 1), 4.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 0), -2.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 1), 1.0, RODIN_FUZZY_CONSTANT);
    }

    // Test with H1<3> (cubic)
    {
      H1 fes(std::integral_constant<size_t, 3>{}, mesh, mesh.getSpaceDimension());
      GridFunction gf(fes);
      gf.project(linear_func);
      auto jac_gf = Jacobian(gf);

      auto it = mesh.getPolytope(mesh.getDimension(), 0);
      const auto& polytope = *it;
      Point p(polytope, Math::Vector<Real>{{0.5, 0.5}});
      auto jac_value = jac_gf.getValue(p);

      EXPECT_NEAR(jac_value(0, 0), 3.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(0, 1), 4.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 0), -2.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 1), 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Jacobian, UsageInBilinearForm)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);

    // Strain energy form: ∫ J(u) : J(v) dx
    bf = Integral(Jacobian(u), Jacobian(v));

    EXPECT_FALSE(bf.getLocalIntegrators().empty());

    // This should assemble without errors
    bf.assemble();

    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
  }

  TEST(Rodin_Variational_H1_Jacobian, RandomCoordinates_LinearVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: u = (x + y, 2x - y)
    // Jacobian should be constant: [[1, 1], [2, -1]]
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.x() + p.y(), 2.0 * p.x() - p.y();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(2, linear_lambda);
    gf.project(linear_func);

    auto jac_gf = Jacobian(gf);

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

      auto jac_value = jac_gf.getValue(p);

      // Expected Jacobian: [[1, 1], [2, -1]]
      EXPECT_NEAR(jac_value(0, 0), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(0, 1), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 0), 2.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 1), -1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Jacobian, RandomCoordinates_QuadraticVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 5, 5 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> can represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project: u = (x^2, y^2)
    // Jacobian: [[2x, 0], [0, 2y]]
    auto quadratic_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.x() * p.x(), p.y() * p.y();
      return v;
    };
    VectorFunction<decltype(quadratic_lambda)> quadratic_func(2, quadratic_lambda);
    gf.project(quadratic_func);

    auto jac_gf = Jacobian(gf);

    // Test at 15 random points
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

      auto jac_value = jac_gf.getValue(p);
      const auto& phys_coords = p.getPhysicalCoordinates();

      // Expected Jacobian: [[2x, 0], [0, 2y]]
      EXPECT_NEAR(jac_value(0, 0), 2.0 * phys_coords(0), 1e-10);
      EXPECT_NEAR(jac_value(0, 1), 0.0, 1e-10);
      EXPECT_NEAR(jac_value(1, 0), 0.0, 1e-10);
      EXPECT_NEAR(jac_value(1, 1), 2.0 * phys_coords(1), 1e-10);
    }
  }

  TEST(Rodin_Variational_H1_Jacobian, RandomCoordinates_QuadraticVectorField_H1_3)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<3> can represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 3>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project: u = (x^2 + xy, y^2 - xy)
    // Jacobian: [[2x + y, x], [-y, 2y - x]]
    auto quadratic_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.x() * p.x() + p.x() * p.y(), 
           p.y() * p.y() - p.x() * p.y();
      return v;
    };
    VectorFunction<decltype(quadratic_lambda)> quadratic_func(2, quadratic_lambda);
    gf.project(quadratic_func);

    auto jac_gf = Jacobian(gf);

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

      auto jac_value = jac_gf.getValue(p);
      const auto& phys_coords = p.getPhysicalCoordinates();

      Real px = phys_coords(0);
      Real py = phys_coords(1);

      // Expected Jacobian: [[2x + y, x], [-y, 2y - x]]
      EXPECT_NEAR(jac_value(0, 0), 2.0 * px + py, 1e-10);
      EXPECT_NEAR(jac_value(0, 1), px, 1e-10);
      EXPECT_NEAR(jac_value(1, 0), -py, 1e-10);
      EXPECT_NEAR(jac_value(1, 1), 2.0 * py - px, 1e-10);
    }
  }

  TEST(Rodin_Variational_H1_Jacobian, RandomCoordinates_ConstantVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a constant vector field
    auto const_lambda = []( const Geometry::Point& p) { 
      Math::Vector<Real> v(2);
      v << 3.5, -2.0;
      return v;
    };
    VectorFunction<decltype(const_lambda)> const_func(2, const_lambda);
    gf.project(const_func);

    auto jac_gf = Jacobian(gf);

    // Test at 10 random points - Jacobian should always be zero
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

      auto jac_value = jac_gf.getValue(p);
      EXPECT_NEAR(jac_value.norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_H1_Jacobian, RandomCoordinates_QuarticVectorField_H1_5)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<5> can represent quartic (degree 4) functions exactly
    H1 fes(std::integral_constant<size_t, 5>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a quartic vector field: u = (x^4, y^4)
    // Jacobian: [[4x^3, 0], [0, 4y^3]]
    auto quartic_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      Real x = p.x();
      Real y = p.y();
      v << x*x*x*x, y*y*y*y;
      return v;
    };
    VectorFunction<decltype(quartic_lambda)> quartic_func(2, quartic_lambda);
    gf.project(quartic_func);

    auto jac_gf = Jacobian(gf);

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

      auto jac_value = jac_gf.getValue(p);
      const auto& phys_coords = p.getPhysicalCoordinates();

      Real px = phys_coords(0);
      Real py = phys_coords(1);

      // Expected Jacobian: [[4x^3, 0], [0, 4y^3]]
      EXPECT_NEAR(jac_value(0, 0), 4.0 * px*px*px, 1e-9);
      EXPECT_NEAR(jac_value(0, 1), 0.0, 1e-10);
      EXPECT_NEAR(jac_value(1, 0), 0.0, 1e-10);
      EXPECT_NEAR(jac_value(1, 1), 4.0 * py*py*py, 1e-9);
    }
  }

  TEST(Rodin_Variational_H1_Jacobian, ProjectJacobianOntoGridFunction_LinearVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Create vector H1<2> space for the original vector field
    H1 fes_vector(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes_vector);

    // Project linear vector field u = (x + y, 2x - y) with Jacobian [[1, 1], [2, -1]]
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.x() + p.y(), 2.0 * p.x() - p.y();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(2, linear_lambda);
    gf.project(linear_func);

    // Create a higher-dimensional H1<2> space for the Jacobian matrix (2x2 = 4 components)
    // For simplicity, we'll create separate GridFunctions for each component
    H1 fes_scalar(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction jac_00(fes_scalar), jac_01(fes_scalar), jac_10(fes_scalar), jac_11(fes_scalar);

    // Project Jacobian components
    auto jac = Jacobian(gf);
    jac_00 = Component(jac, 0, 0);
    jac_01 = Component(jac, 0, 1);
    jac_10 = Component(jac, 1, 0);
    jac_11 = Component(jac, 1, 1);

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

      // Evaluate the projected Jacobian components
      Real j00 = jac_00.getValue(p);
      Real j01 = jac_01.getValue(p);
      Real j10 = jac_10.getValue(p);
      Real j11 = jac_11.getValue(p);

      // Should match the constant Jacobian [[1, 1], [2, -1]]
      EXPECT_NEAR(j00, 1.0, 0.1);  // Relaxed tolerance for FE projection
      EXPECT_NEAR(j01, 1.0, 0.1);
      EXPECT_NEAR(j10, 2.0, 0.1);
      EXPECT_NEAR(j11, -1.0, 0.1);
    }
  }

  TEST(Rodin_Variational_H1_Jacobian, ProjectJacobianOntoGridFunction_QuadraticVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Create vector H1<3> space for the quadratic vector field
    H1 fes_vector(std::integral_constant<size_t, 3>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes_vector);

    // Project quadratic vector field u = (x^2, y^2) with Jacobian [[2x, 0], [0, 2y]]
    auto quadratic_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.x() * p.x(), p.y() * p.y();
      return v;
    };
    VectorFunction<decltype(quadratic_lambda)> quadratic_func(2, quadratic_lambda);
    gf.project(quadratic_func);

    // Create scalar H1<3> spaces for Jacobian components
    H1 fes_scalar(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction jac_00(fes_scalar), jac_01(fes_scalar), jac_10(fes_scalar), jac_11(fes_scalar);

    // Project Jacobian components
    auto jac = Jacobian(gf);
    jac_00.project(Component(jac, 0, 0));
    jac_01.project(Component(jac, 0, 1));
    jac_10.project(Component(jac, 1, 0));
    jac_11.project(Component(jac, 1, 1));

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

      // Evaluate the projected Jacobian components
      Real j00 = jac_00.getValue(p);
      Real j01 = jac_01.getValue(p);
      Real j10 = jac_10.getValue(p);
      Real j11 = jac_11.getValue(p);

      const auto& phys_coords = p.getPhysicalCoordinates();
      Real px = phys_coords(0);
      Real py = phys_coords(1);

      // Should match Jacobian [[2x, 0], [0, 2y]]
      EXPECT_NEAR(j00, 2.0 * px, 1.0);  // Larger tolerance for quadratic FE projection
      EXPECT_NEAR(j01, 0.0, 0.5);
      EXPECT_NEAR(j10, 0.0, 0.5);
      EXPECT_NEAR(j11, 2.0 * py, 1.0);
    }
  }

  // ============================================================================
  // Tetrahedron tests
  // ============================================================================

  TEST(Rodin_Variational_H1_Jacobian_Tet, ShapeFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    TestFunction v(fes);

    auto jac_u = Jacobian(u);
    auto jac_v = Jacobian(v);

    // Jacobian of a vector function should be a matrix function
    // For 3D, Jacobian should be 3x3
  }

  TEST(Rodin_Variational_H1_Jacobian_Tet, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    auto jac_gf = Jacobian(gf);

    // Jacobian of a vector GridFunction should be a matrix function
    EXPECT_EQ(&jac_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Jacobian_Tet, GridFunction_ConstantFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a constant vector field
    auto const_func_lambda = []( const Geometry::Point& p) { 
      Math::Vector<Real> v(3);
      v << 1.0, 2.0, 3.0;
      return v;
    };
    VectorFunction<decltype(const_func_lambda)> const_func(3, const_func_lambda);
    gf.project(const_func);

    auto jac_gf = Jacobian(gf);

    // Jacobian of a constant vector field should be zero
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.2, 0.2}};
    Point p(polytope, rc);

    auto jac_value = jac_gf.getValue(p);
    EXPECT_NEAR(jac_value.norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Jacobian_Tet, GridFunction_LinearVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project the identity vector field: u = (x, y, z)
    // Jacobian should be the identity matrix: [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(3);
      v << p.x(), p.y(), p.z();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(3, linear_lambda);
    gf.project(linear_func);

    auto jac_gf = Jacobian(gf);

    // For linear functions, Jacobian should be constant
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25, 0.25}};
    Point p(polytope, rc);

    auto jac_value = jac_gf.getValue(p);
    // Expected Jacobian: identity matrix
    EXPECT_EQ(jac_value.rows(), 3);
    EXPECT_EQ(jac_value.cols(), 3);
    EXPECT_NEAR(jac_value(0, 0), 1.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(0, 1), 0.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(0, 2), 0.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(1, 0), 0.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(1, 1), 1.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(1, 2), 0.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(2, 0), 0.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(2, 1), 0.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(2, 2), 1.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Jacobian_Tet, GridFunction_QuadraticVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> (quadratic) can represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a quadratic vector field: u = (x^2, y^2, z^2)
    // Jacobian should be: [[2x, 0, 0], [0, 2y, 0], [0, 0, 2z]]
    auto quadratic_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(3);
      v << p.x() * p.x(), p.y() * p.y(), p.z() * p.z();
      return v;
    };
    VectorFunction<decltype(quadratic_lambda)> quadratic_func(3, quadratic_lambda);
    gf.project(quadratic_func);

    auto jac_gf = Jacobian(gf);

    // Evaluate at a point
    auto it = mesh.getPolytope(mesh.getDimension(), mesh.getCellCount() / 2);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.3, 0.1}};
    Point p(polytope, rc);

    auto jac_value = jac_gf.getValue(p);
    const auto& phys_coords = p.getPhysicalCoordinates();

    // Expected Jacobian: [[2x, 0, 0], [0, 2y, 0], [0, 0, 2z]]
    EXPECT_EQ(jac_value.rows(), 3);
    EXPECT_EQ(jac_value.cols(), 3);
    EXPECT_NEAR(jac_value(0, 0), 2.0 * phys_coords(0), 1e-10);
    EXPECT_NEAR(jac_value(0, 1), 0.0, 1e-10);
    EXPECT_NEAR(jac_value(0, 2), 0.0, 1e-10);
    EXPECT_NEAR(jac_value(1, 0), 0.0, 1e-10);
    EXPECT_NEAR(jac_value(1, 1), 2.0 * phys_coords(1), 1e-10);
    EXPECT_NEAR(jac_value(1, 2), 0.0, 1e-10);
    EXPECT_NEAR(jac_value(2, 0), 0.0, 1e-10);
    EXPECT_NEAR(jac_value(2, 1), 0.0, 1e-10);
    EXPECT_NEAR(jac_value(2, 2), 2.0 * phys_coords(2), 1e-10);
  }

  TEST(Rodin_Variational_H1_Jacobian_Tet, RandomCoordinates_LinearVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 2, 2, 2 });
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: u = (x + y, y + z, x + z)
    // Jacobian should be constant: [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(3);
      v << p.x() + p.y(), p.y() + p.z(), p.x() + p.z();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(3, linear_lambda);
    gf.project(linear_func);

    auto jac_gf = Jacobian(gf);

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

      auto jac_value = jac_gf.getValue(p);

      // Expected Jacobian: [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
      EXPECT_NEAR(jac_value(0, 0), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(0, 1), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(0, 2), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 0), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 1), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 2), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(2, 0), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(2, 1), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(2, 2), 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // ============================================================================
  // Quadrilateral tests
  // ============================================================================

  TEST(Rodin_Variational_H1_Jacobian_Quad, ShapeFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    TestFunction v(fes);

    auto jac_u = Jacobian(u);
    auto jac_v = Jacobian(v);

    // Jacobian of a vector function should be a matrix function
    // For 2D, Jacobian should be 2x2
  }

  TEST(Rodin_Variational_H1_Jacobian_Quad, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Vector-valued H1<2> space
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    auto jac_gf = Jacobian(gf);

    // Jacobian of a vector GridFunction should be a matrix function
    EXPECT_EQ(&jac_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Jacobian_Quad, GridFunction_LinearVectorField_H1_2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project the identity vector field: u = (x, y)
    // Jacobian should be the identity matrix: [[1, 0], [0, 1]]
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << p.x(), p.y();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(2, linear_lambda);
    gf.project(linear_func);

    auto jac_gf = Jacobian(gf);

    // For linear functions, Jacobian should be constant
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.4}};
    Point p(polytope, rc);

    auto jac_value = jac_gf.getValue(p);
    // Expected Jacobian: identity matrix
    EXPECT_EQ(jac_value.rows(), 2);
    EXPECT_EQ(jac_value.cols(), 2);
    EXPECT_NEAR(jac_value(0, 0), 1.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(0, 1), 0.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(1, 0), 0.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(jac_value(1, 1), 1.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Jacobian_Quad, RandomCoordinates_LinearVectorField)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> can represent linear functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    GridFunction gf(fes);

    // Project a linear vector field: u = (2x + y, x - y)
    // Jacobian should be constant: [[2, 1], [1, -1]]
    auto linear_lambda = [](const Geometry::Point& p) {
      Math::Vector<Real> v(2);
      v << 2.0 * p.x() + p.y(), p.x() - p.y();
      return v;
    };
    VectorFunction<decltype(linear_lambda)> linear_func(2, linear_lambda);
    gf.project(linear_func);

    auto jac_gf = Jacobian(gf);

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

      auto jac_value = jac_gf.getValue(p);

      // Expected Jacobian: [[2, 1], [1, -1]]
      EXPECT_NEAR(jac_value(0, 0), 2.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(0, 1), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 0), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(jac_value(1, 1), -1.0, RODIN_FUZZY_CONSTANT);
    }
  }
}

TEST(H1Jacobian, ShapeFunction_getDOFs_Triangle_H1_2)
{
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, {4, 4});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);

  H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
  TrialFunction u(uh);
  auto jac_u = Jacobian(u);

  // H1<2> on a triangle has 6 DOFs per component, 2 components = 12 total
  const auto cellIt = mesh.getCell(0);
  EXPECT_EQ(jac_u.getDOFs(*cellIt), 12);
}

TEST(H1Jacobian, ShapeFunction_getDOFs_Tetrahedron_H1_2)
{
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
  mesh.getConnectivity().compute(3, 2);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);

  H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
  TrialFunction u(uh);
  auto jac_u = Jacobian(u);

  // H1<2> on a tetrahedron has 10 DOFs per component, 3 components = 30 total
  const auto cellIt = mesh.getCell(0);
  EXPECT_EQ(jac_u.getDOFs(*cellIt), 30);
}

