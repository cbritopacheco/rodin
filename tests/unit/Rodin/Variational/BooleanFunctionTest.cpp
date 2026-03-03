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
  TEST(Rodin_Variational_BooleanFunction, TrueConstant_Construction)
  {
    BooleanFunction bf(true);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_TRUE(bf.getValue(p));
  }

  TEST(Rodin_Variational_BooleanFunction, FalseConstant_Construction)
  {
    BooleanFunction bf(false);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);

    EXPECT_FALSE(bf.getValue(p));
  }

  TEST(Rodin_Variational_BooleanFunction, CopyConstructor_True)
  {
    BooleanFunction bf(true);
    BooleanFunction bf_copy(bf);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);

    EXPECT_EQ(bf_copy.getValue(p), bf.getValue(p));
    EXPECT_TRUE(bf_copy.getValue(p));
  }

  TEST(Rodin_Variational_BooleanFunction, CopyConstructor_False)
  {
    BooleanFunction bf(false);
    BooleanFunction bf_copy(bf);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);

    EXPECT_EQ(bf_copy.getValue(p), bf.getValue(p));
    EXPECT_FALSE(bf_copy.getValue(p));
  }

  TEST(Rodin_Variational_BooleanFunction, MoveConstructor_True)
  {
    BooleanFunction bf(true);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    Boolean original_value = bf.getValue(p);

    BooleanFunction bf_moved(std::move(bf));
    EXPECT_EQ(bf_moved.getValue(p), original_value);
    EXPECT_TRUE(bf_moved.getValue(p));
  }

  TEST(Rodin_Variational_BooleanFunction, MoveConstructor_False)
  {
    BooleanFunction bf(false);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.75}};
    Point p(polytope, rc);

    Boolean original_value = bf.getValue(p);

    BooleanFunction bf_moved(std::move(bf));
    EXPECT_EQ(bf_moved.getValue(p), original_value);
    EXPECT_FALSE(bf_moved.getValue(p));
  }

  TEST(Rodin_Variational_BooleanFunction, Copy_True)
  {
    BooleanFunction bf(true);
    auto copied = bf.copy();

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NE(copied, nullptr);
    EXPECT_EQ(copied->getValue(p), bf.getValue(p));
    EXPECT_TRUE(copied->getValue(p));

    delete copied;
  }

  TEST(Rodin_Variational_BooleanFunction, Copy_False)
  {
    BooleanFunction bf(false);
    auto copied = bf.copy();

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);

    EXPECT_NE(copied, nullptr);
    EXPECT_EQ(copied->getValue(p), bf.getValue(p));
    EXPECT_FALSE(copied->getValue(p));

    delete copied;
  }

  TEST(Rodin_Variational_BooleanFunction, ConstantValue_MultiplePoints)
  {
    BooleanFunction bf_true(true);
    BooleanFunction bf_false(false);

    // Create multiple points to test that the constant value is returned everywhere
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;

    std::vector<Math::Vector<Real>> test_coords =
    {
      Math::Vector<Real>{{0.0, 0.0}},
      Math::Vector<Real>{{1.0, 0.0}},
      Math::Vector<Real>{{0.0, 1.0}},
      Math::Vector<Real>{{0.5, 0.5}},
      Math::Vector<Real>{{0.33, 0.67}}
    };

    for (const auto& rc : test_coords)
    {
      Point p(polytope, rc);
      EXPECT_TRUE(bf_true.getValue(p));
      EXPECT_FALSE(bf_false.getValue(p));
    }
  }

  TEST(Rodin_Variational_BooleanFunction, PolymorphicUsage)
  {
    // Test that we can use BooleanFunction polymorphically as BooleanFunctionBase
    std::unique_ptr<BooleanFunctionBase<BooleanFunction<Boolean>>> bool_func = 
      std::make_unique<BooleanFunction<Boolean>>(true);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_TRUE(bool_func->getValue(p));
  }

  TEST(Rodin_Variational_BooleanFunction, TraceToBoundary)
  {
    const Attribute interior_attr = 1;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    for (auto it = mesh.getCell(); !it.end(); ++it)
      mesh.setAttribute(it.key(), interior_attr);
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D - 1, D);

    BooleanFunction bf(true);
    auto traced_bf = bf.traceOf(interior_attr);

    // Create a boundary point for testing
    auto it = mesh.getPolytope(D - 1, 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5}};
    Point p(polytope, rc);

    EXPECT_TRUE(traced_bf.getValue(p));
  }

  TEST(Rodin_Variational_BooleanFunction, CTAD_True)
  {
    // Test Class Template Argument Deduction
    auto bf = BooleanFunction(true);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);

    EXPECT_TRUE(bf.getValue(p));
  }

  TEST(Rodin_Variational_BooleanFunction, CTAD_False)
  {
    // Test Class Template Argument Deduction
    auto bf = BooleanFunction(false);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);

    EXPECT_FALSE(bf.getValue(p));
  }
}
