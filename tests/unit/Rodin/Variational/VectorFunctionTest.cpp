#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_VectorFunction, ConstantVector_2D_Construction)
  {
    VectorFunction vf{1.0, 2.0};
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    auto value = vf.getValue(p);
    EXPECT_NEAR(value(0), 1.0, 1e-10);
    EXPECT_NEAR(value(1), 2.0, 1e-10);
    EXPECT_EQ(vf.getDimension(), 2);
  }

  TEST(Rodin_Variational_VectorFunction, ConstantVector_3D_Construction)
  {
    VectorFunction vf{3.14, -2.71, 1.41};
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    auto value = vf.getValue(p);
    EXPECT_NEAR(value(0), 3.14, 1e-10);
    EXPECT_NEAR(value(1), -2.71, 1e-10);
    EXPECT_NEAR(value(2), 1.41, 1e-10);
    EXPECT_EQ(vf.getDimension(), 3);
  }

  TEST(Rodin_Variational_VectorFunction, MixedTypes_Construction)
  {
    VectorFunction vf{1.0, 2.5, 3.0};
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    auto value = vf.getValue(p);
    EXPECT_NEAR(value(0), 1.0, 1e-10);
    EXPECT_NEAR(value(1), 2.5, 1e-10);
    EXPECT_NEAR(value(2), 3.0, 1e-10);
    EXPECT_EQ(vf.getDimension(), 3);
  }

  TEST(Rodin_Variational_VectorFunction, ComponentAccess_2D)
  {
    VectorFunction vf{10.0, 20.0};
    auto x_comp = vf.x();
    auto y_comp = vf.y();
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    EXPECT_NEAR(x_comp.getValue(p), 10.0, 1e-10);
    EXPECT_NEAR(y_comp.getValue(p), 20.0, 1e-10);
  }

  TEST(Rodin_Variational_VectorFunction, ComponentAccess_3D)
  {
    VectorFunction vf{100.0, 200.0, 300.0};
    auto x_comp = vf.x();
    auto y_comp = vf.y();
    auto z_comp = vf.z();
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);
    EXPECT_NEAR(x_comp.getValue(p), 100.0, 1e-10);
    EXPECT_NEAR(y_comp.getValue(p), 200.0, 1e-10);
    EXPECT_NEAR(z_comp.getValue(p), 300.0, 1e-10);
  }

  TEST(Rodin_Variational_VectorFunction, IndexAccess)
  {
    VectorFunction vf{5.0, 15.0, 25.0};
    auto comp0 = vf(0);
    auto comp1 = vf(1);
    auto comp2 = vf(2);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    EXPECT_NEAR(comp0.getValue(p), 5.0, 1e-10);
    EXPECT_NEAR(comp1.getValue(p), 15.0, 1e-10);
    EXPECT_NEAR(comp2.getValue(p), 25.0, 1e-10);
  }

  TEST(Rodin_Variational_VectorFunction, CopyConstructor)
  {
    VectorFunction vf{7.5, -12.3};
    VectorFunction vf_copy(vf);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);
    auto value_orig = vf.getValue(p);
    auto value_copy = vf_copy.getValue(p);
    EXPECT_NEAR(value_orig(0), value_copy(0), 1e-10);
    EXPECT_NEAR(value_orig(1), value_copy(1), 1e-10);
    EXPECT_EQ(vf_copy.getDimension(), vf.getDimension());
  }

  TEST(Rodin_Variational_VectorFunction, MoveConstructor)
  {
    VectorFunction vf{42.0, 84.0};
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);
    auto original_value = vf.getValue(p);
    VectorFunction vf_moved(std::move(vf));
    auto moved_value = vf_moved.getValue(p);
    EXPECT_NEAR(original_value(0), moved_value(0), 1e-10);
    EXPECT_NEAR(original_value(1), moved_value(1), 1e-10);
    EXPECT_EQ(vf_moved.getDimension(), 2);
  }

  TEST(Rodin_Variational_VectorFunction, Copy)
  {
    VectorFunction vf{99.9, -77.7, 55.5};
    auto copied = vf.copy();
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.75}};
    Point p(polytope, rc);
    EXPECT_NE(copied, nullptr);
    auto value_orig = vf.getValue(p);
    auto value_copy = copied->getValue(p);
    EXPECT_NEAR(value_orig(0), value_copy(0), 1e-10);
    EXPECT_NEAR(value_orig(1), value_copy(1), 1e-10);
    EXPECT_NEAR(value_orig(2), value_copy(2), 1e-10);
    EXPECT_EQ(copied->getDimension(), vf.getDimension());
    delete copied;
  }

  TEST(Rodin_Variational_VectorFunction, TraceToBoundary)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D - 1, D);
    VectorFunction vf{1.0, 2.0};
    vf.traceOf(RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
    // Create a boundary point for testing
    auto it = mesh.getPolytope(D - 1, 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(D - 1, 0);
    const Math::Vector<Real> rc{{0.5}};
    Point p(polytope, rc);
    auto value = vf.getValue(p);
    EXPECT_NEAR(value(0), 1.0, 1e-10);
    EXPECT_NEAR(value(1), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_VectorFunction, ZeroVector)
  {
    VectorFunction vf{0.0, 0.0, 0.0};
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);
    auto value = vf.getValue(p);
    EXPECT_NEAR(value.norm(), 0.0, 1e-10);
  }
}
