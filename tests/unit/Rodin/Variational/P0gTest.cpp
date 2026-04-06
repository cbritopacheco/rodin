/*
 * Unit tests for P0g (global constant) finite element space.
 *
 * Tests cover:
 *   - Scalar P0g space construction, size, and DOFs
 *   - Vector P0g space construction, size, and DOFs
 *   - P0g element properties (constant basis, order 0)
 *   - GridFunction projection onto P0g
 *   - Grad of P0g GridFunction is identically zero
 *   - Div of P0g vector GridFunction is identically zero
 *   - Jacobian of P0g vector GridFunction is identically zero
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Assembly.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  // ---- Scalar P0g space ----

  TEST(Rodin_Variational_P0g, ScalarConstruction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P0g Vh(mesh);

    // Scalar P0g has exactly 1 global DOF
    EXPECT_EQ(Vh.getSize(), 1u);
    EXPECT_EQ(Vh.getVectorDimension(), 1u);
  }

  TEST(Rodin_Variational_P0g, ScalarGridFunctionProjection)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P0g Vh(mesh);

    // Project a constant function
    GridFunction u(Vh);
    u = RealFunction(3.14);

    // The single DOF should hold the projected value
    EXPECT_EQ(u.getData().size(), 1);
    EXPECT_NEAR(u.getData()(0), 3.14, 1e-10);
  }

  TEST(Rodin_Variational_P0g, ScalarDOFsMapping)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P0g Vh(mesh);

    // Every cell maps to the same global DOF 0
    const auto& dofs = Vh.getDOFs(mesh.getDimension(), 0);
    EXPECT_EQ(dofs.size(), 1);
    EXPECT_EQ(dofs[0], 0);

    // Different cell, same DOF
    const auto& dofs2 = Vh.getDOFs(mesh.getDimension(), 1);
    EXPECT_EQ(dofs2.size(), 1);
    EXPECT_EQ(dofs2[0], 0);
  }

  // ---- Vector P0g space ----

  TEST(Rodin_Variational_P0g, VectorConstruction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P0g<Math::Vector<Real>, Geometry::Mesh<Context::Local>> Vh(mesh, 2);

    // Vector P0g in 2D has vdim=2 DOFs
    EXPECT_EQ(Vh.getSize(), 2u);
    EXPECT_EQ(Vh.getVectorDimension(), 2u);
  }

  TEST(Rodin_Variational_P0g, VectorDOFsMapping)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P0g<Math::Vector<Real>, Geometry::Mesh<Context::Local>> Vh(mesh, 2);

    // Every cell maps to DOFs {0, 1}
    const auto& dofs = Vh.getDOFs(mesh.getDimension(), 0);
    EXPECT_EQ(dofs.size(), 2);
    EXPECT_EQ(dofs[0], 0);
    EXPECT_EQ(dofs[1], 1);
  }

  // ---- P0g element properties ----

  TEST(Rodin_Variational_P0g, ScalarElementProperties)
  {
    P0gElement<Real> elem(Polytope::Type::Triangle);

    // Scalar P0g element has 1 basis function
    EXPECT_EQ(elem.getCount(), 1u);

    // Polynomial order is 0
    EXPECT_EQ(elem.getOrder(), 0u);

    // Basis function evaluates to 1 everywhere
    Math::Vector<Real> pt(2);
    pt << 0.25, 0.25;
    EXPECT_NEAR(elem.getBasis(0)(pt), 1.0, 1e-14);
  }

  // ---- Grad of P0g scalar GridFunction is zero ----

  TEST(Rodin_Variational_P0g, GradOfScalarGridFunctionIsZero)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);
    P0g Vh(mesh);

    GridFunction u(Vh);
    u = RealFunction(5.0);

    // Integrate |Grad(u)|^2 over the domain — should be zero for constant u
    P1 Wh(mesh);
    GridFunction grad_norm(Wh);
    grad_norm = Frobenius(Grad(u));

    auto integral_val = Integral(grad_norm).compute();
    EXPECT_NEAR(integral_val, 0.0, 1e-14);
  }

  TEST(Rodin_Variational_P0g, GradOrder)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P0g Vh(mesh);

    GridFunction u(Vh);
    u = RealFunction(1.0);

    auto grad_u = Grad(u);

    // Grad of P0g has order 0 (zero polynomial)
    const auto& polytope = *mesh.getCell();
    auto order = grad_u.getOrder(polytope);
    ASSERT_TRUE(order.has_value());
    EXPECT_EQ(*order, 0u);
  }

  // ---- Div of P0g vector GridFunction is zero ----

  TEST(Rodin_Variational_P0g, DivOfVectorGridFunctionIsZero)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);
    P0g<Math::Vector<Real>, Geometry::Mesh<Context::Local>> Vh(mesh, 2);

    GridFunction u(Vh);
    // Set to constant vector [1, 2]
    u.getData().resize(2);
    u.getData()(0) = 1.0;
    u.getData()(1) = 2.0;

    // Integrate |Div(u)|^2 — should be zero for constant vector
    P1 Wh(mesh);
    GridFunction div_norm(Wh);
    div_norm = Pow(Div(u), 2);

    auto integral_val = Integral(div_norm).compute();
    EXPECT_NEAR(integral_val, 0.0, 1e-14);
  }

  // ---- Jacobian of P0g vector GridFunction is zero ----

  TEST(Rodin_Variational_P0g, JacobianOfVectorGridFunctionIsZero)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);
    P0g<Math::Vector<Real>, Geometry::Mesh<Context::Local>> Vh(mesh, 2);

    GridFunction u(Vh);
    u.getData().resize(2);
    u.getData()(0) = 1.0;
    u.getData()(1) = 2.0;

    // Integrate |Jacobian(u)|^2 — should be zero for constant vector
    P1 Wh(mesh);
    GridFunction jac_norm(Wh);
    jac_norm = Frobenius(Jacobian(u));

    auto integral_val = Integral(jac_norm).compute();
    EXPECT_NEAR(integral_val, 0.0, 1e-14);
  }
}
