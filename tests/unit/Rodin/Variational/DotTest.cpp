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
  TEST(Rodin_Variational_Dot, VectorFunctions_Construction)
  {
    VectorFunction vf1{1.0, 2.0};
    VectorFunction vf2{3.0, 4.0};

    auto dot_product = Dot(vf1, vf2);
  }

  TEST(Rodin_Variational_Dot, VectorFunctions_Value)
  {
    VectorFunction vf1{3.0, 4.0};
    VectorFunction vf2{2.0, 1.0};

    auto dot_product = Dot(vf1, vf2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    Real value = dot_product.getValue(p);
    // (3,4) · (2,1) = 3*2 + 4*1 = 6 + 4 = 10
    EXPECT_NEAR(value, 10.0, 1e-10);
  }

  TEST(Rodin_Variational_Dot, ShapeFunctions_Construction)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);

    auto dot_uv = Dot(u, v);
  }

  TEST(Rodin_Variational_Dot, Gradients_Construction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    auto dot_grads = Dot(Grad(u), Grad(v));
  }

  TEST(Rodin_Variational_Dot, Copy)
  {
    VectorFunction vf1{1.0, 2.0};
    VectorFunction vf2{3.0, 4.0};

    auto dot_product = Dot(vf1, vf2);
    auto copied = dot_product.copy();

    EXPECT_NE(copied, nullptr);

    delete copied;
  }

  TEST(Rodin_Variational_Dot, WithVectorFunctionAndShapeFunction)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TestFunction v(fes);

    VectorFunction force{1.0, 2.0};
    auto dot_force_v = Dot(force, v);
  }

  TEST(Rodin_Variational_Dot, InLinearForm)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TestFunction v(fes);
    LinearForm lf(v);

    VectorFunction force{5.0, 10.0};
    lf = Integral(Dot(force, v));
    lf.assemble();

    const auto& vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);
  }

  TEST(Rodin_Variational_Dot, InBilinearForm)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);

    // Vector mass matrix: ∫ u · v dx
    bf = Integral(Dot(u, v));
    bf.assemble();

    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
  }

  TEST(Rodin_Variational_Dot, ZeroVectors)
  {
    VectorFunction vf1{0.0, 0.0};
    VectorFunction vf2{1.0, 2.0};

    auto dot_product = Dot(vf1, vf2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    Real value = dot_product.getValue(p);
    EXPECT_NEAR(value, 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Dot, OrthogonalVectors)
  {
    VectorFunction vf1{1.0, 0.0};
    VectorFunction vf2{0.0, 1.0};

    auto dot_product = Dot(vf1, vf2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);

    Real value = dot_product.getValue(p);
    EXPECT_NEAR(value, 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Dot, ParallelVectors)
  {
    VectorFunction vf1{2.0, 3.0};
    VectorFunction vf2{4.0, 6.0};  // 2 * vf1

    auto dot_product = Dot(vf1, vf2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.75}};
    Point p(polytope, rc);

    Real value = dot_product.getValue(p);
    // (2,3) · (4,6) = 2*4 + 3*6 = 8 + 18 = 26
    EXPECT_NEAR(value, 26.0, 1e-10);
  }

  TEST(Rodin_Variational_Dot, 3D_Vectors)
  {
    VectorFunction vf1{1.0, 2.0, 3.0};
    VectorFunction vf2{4.0, 5.0, 6.0};

    auto dot_product = Dot(vf1, vf2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);

    Real value = dot_product.getValue(p);
    // (1,2,3) · (4,5,6) = 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    EXPECT_NEAR(value, 32.0, 1e-10);
  }

  TEST(Rodin_Variational_Dot, CopyConstructor)
  {
    VectorFunction vf1{7.5, -2.3};
    VectorFunction vf2{1.2, 4.8};

    auto dot_product = Dot(vf1, vf2);
    auto dot_copy(dot_product);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    Real value_orig = dot_product.getValue(p);
    Real value_copy = dot_copy.getValue(p);

    EXPECT_NEAR(value_orig, value_copy, 1e-10);
  }
}
