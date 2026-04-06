#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/H1/Derivative.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_H1_Derivative, GridFunction_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);
    gf = RealFunction(1.0);

    auto dx_gf = Dx(gf);
    auto dy_gf = Dy(gf);

    EXPECT_EQ(&dx_gf.getOperand(), &gf);
    EXPECT_EQ(&dy_gf.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Derivative, GridFunction_ConstantFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    gf = RealFunction(42.0);

    auto dx_gf = Dx(gf);
    auto dy_gf = Dy(gf);

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(dx_gf.getValue(p), 0.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(dy_gf.getValue(p), 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Derivative, GridFunction_LinearFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // f(x,y) = 3x + 4y
    RealFunction linear_func(
        [](const Geometry::Point& p) { return 3.0 * p.x() + 4.0 * p.y(); });
    gf.project(linear_func);

    auto dx_gf = Dx(gf);
    auto dy_gf = Dy(gf);

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);

    // df/dx = 3, df/dy = 4
    EXPECT_NEAR(dx_gf.getValue(p), 3.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(dy_gf.getValue(p), 4.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_H1_Derivative, GridFunction_QuadraticFunction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 8, 8 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // H1<2> can represent quadratic functions exactly
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // f(x,y) = x^2 + 2xy
    RealFunction quad_func(
        [](const Geometry::Point& p) { return p.x() * p.x() + 2.0 * p.x() * p.y(); });
    gf.project(quad_func);

    auto dx_gf = Dx(gf);
    auto dy_gf = Dy(gf);

    // Evaluate at the centroid of an interior element
    auto it = mesh.getPolytope(mesh.getDimension(), 30);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0 / 3.0, 1.0 / 3.0}};
    Point p(polytope, rc);

    // df/dx = 2x + 2y, df/dy = 2x
    const Real x = p.x();
    const Real y = p.y();
    EXPECT_NEAR(dx_gf.getValue(p), 2.0 * x + 2.0 * y, 1e-4);
    EXPECT_NEAR(dy_gf.getValue(p), 2.0 * x, 1e-4);
  }

  TEST(Rodin_Variational_H1_Derivative, GridFunction_Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);
    gf = RealFunction(1.0);

    auto dx_gf = Dx(gf);
    auto copied = dx_gf;

    EXPECT_EQ(&copied.getOperand(), &gf);
  }

  TEST(Rodin_Variational_H1_Derivative, GridFunction_GetOrder)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);
    gf = RealFunction(1.0);

    auto dx_gf = Dx(gf);

    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    auto order = dx_gf.getOrder(polytope);

    // H1<2> has degree 2, derivative should be degree 1
    ASSERT_TRUE(order.has_value());
    EXPECT_EQ(*order, 1);
  }

  TEST(Rodin_Variational_H1_Derivative, DxDyConsistentWithGrad)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // f(x,y) = x^2 + y^2
    RealFunction func(
        [](const Geometry::Point& p) { return p.x() * p.x() + p.y() * p.y(); });
    gf.project(func);

    auto dx_gf = Dx(gf);
    auto dy_gf = Dy(gf);
    auto grad_gf = Grad(gf);

    // Check at several elements that Dx and Dy match Grad components
    for (Index idx = 0; idx < 4; ++idx)
    {
      auto it = mesh.getPolytope(mesh.getDimension(), idx);
      const auto& polytope = *it;
      const Math::Vector<Real> rc{{0.25, 0.25}};
      Point p(polytope, rc);

      auto grad_val = grad_gf.getValue(p);
      EXPECT_NEAR(dx_gf.getValue(p), grad_val(0), RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(dy_gf.getValue(p), grad_val(1), RODIN_FUZZY_CONSTANT);
    }
  }
}
