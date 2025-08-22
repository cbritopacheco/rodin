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
  TEST(Rodin_Variational_Component, VectorFunction_ComponentAccess)
  {
    VectorFunction vf{10.0, 20.0, 30.0};

    auto comp0 = Component(vf, 0);
    auto comp1 = Component(vf, 1);
    auto comp2 = Component(vf, 2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(comp0.getValue(p), 10.0, 1e-10);
    EXPECT_NEAR(comp1.getValue(p), 20.0, 1e-10);
    EXPECT_NEAR(comp2.getValue(p), 30.0, 1e-10);
  }

  TEST(Rodin_Variational_Component, VectorShapeFunction_ComponentAccess)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);

    auto u_x = u.x();  // This internally uses Component
    auto u_y = u.y();  // This internally uses Component
  }

  TEST(Rodin_Variational_Component, Copy)
  {
    VectorFunction vf{5.5, -3.3};
    auto comp = Component(vf, 1);
    auto copied = comp.copy();

    EXPECT_NE(copied, nullptr);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(copied->getValue(p), comp.getValue(p), 1e-10);

    delete copied;
  }

  TEST(Rodin_Variational_Component, CopyConstructor)
  {
    VectorFunction vf{7.0, 8.0};
    auto comp = Component(vf, 0);
    auto comp_copy(comp);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(comp_copy.getValue(p), comp.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Component, MoveConstructor)
  {
    VectorFunction vf{42.0, 84.0};
    auto comp = Component(vf, 1);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);

    Real original_value = comp.getValue(p);

    auto comp_moved(std::move(comp));
    EXPECT_NEAR(comp_moved.getValue(p), original_value, 1e-10);
  }

  TEST(Rodin_Variational_Component, GetOperand)
  {
    VectorFunction vf{11.0, 22.0};
    auto comp = Component(vf, 0);

    const auto& operand = comp.getOperand();

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    auto operand_value = operand.getValue(p);
    EXPECT_NEAR(operand_value(0), 11.0, 1e-10);
    EXPECT_NEAR(operand_value(1), 22.0, 1e-10);
  }

  TEST(Rodin_Variational_Component, UsageInLinearForm)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TestFunction v(fes);
    LinearForm lf(v);

    // Apply a force only in the x-direction
    RealFunction force_x(5.0);
    lf = Integral(force_x, v.x());
    lf.assemble();

    const auto& vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);
  }

  TEST(Rodin_Variational_Component, UsageInBilinearForm)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm bf(u, v);

    // Mass matrix for x-component only
    bf = Integral(u.x(), v.x());
    bf.assemble();

    const auto& op = bf.getOperator();
    EXPECT_GT(op.rows(), 0);
    EXPECT_GT(op.cols(), 0);
  }

  TEST(Rodin_Variational_Component, ZeroComponent)
  {
    VectorFunction vf{0.0, 15.0, 0.0};
    auto comp0 = Component(vf, 0);
    auto comp2 = Component(vf, 2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.75}};
    Point p(polytope, rc);

    EXPECT_NEAR(comp0.getValue(p), 0.0, 1e-10);
    EXPECT_NEAR(comp2.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Component, NegativeComponents)
  {
    VectorFunction vf{-1.5, -2.5, -3.5};
    auto comp0 = Component(vf, 0);
    auto comp1 = Component(vf, 1);
    auto comp2 = Component(vf, 2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);

    EXPECT_NEAR(comp0.getValue(p), -1.5, 1e-10);
    EXPECT_NEAR(comp1.getValue(p), -2.5, 1e-10);
    EXPECT_NEAR(comp2.getValue(p), -3.5, 1e-10);
  }

  TEST(Rodin_Variational_Component, LargeVectorDimension)
  {
    VectorFunction vf{1.0, 2.0, 3.0, 4.0, 5.0};

    std::vector<Real> expected_values = {1.0, 2.0, 3.0, 4.0, 5.0};

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    for (size_t i = 0; i < expected_values.size(); i++)
    {
      auto comp = Component(vf, i);
      EXPECT_NEAR(comp.getValue(p), expected_values[i], 1e-10);
    }
  }
}
