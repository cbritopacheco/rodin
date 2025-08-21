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
  TEST(Rodin_Variational_ScalarFunction, RealFunction_IsScalarFunction)
  {
    RealFunction rf(3.14);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, trans, rc);
    
    // RealFunction inherits from ScalarFunctionBase
    EXPECT_NEAR(rf.getValue(p), 3.14, 1e-10);
  }

  TEST(Rodin_Variational_ScalarFunction, RealFunction_CopyOperation)
  {
    RealFunction rf(2.718);
    auto copied = rf.copy();
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, trans, rc);
    
    EXPECT_NE(copied, nullptr);
    EXPECT_NEAR(copied->getValue(p), rf.getValue(p), 1e-10);
    
    delete copied;
  }

  TEST(Rodin_Variational_ScalarFunction, RealFunction_PolymorphicUsage)
  {
    // Test that we can use RealFunction polymorphically as ScalarFunctionBase
    std::unique_ptr<ScalarFunctionBase<Real, RealFunction<Real>>> scalar_func = 
      std::make_unique<RealFunction<Real>>(42.0);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, trans, rc);
    
    EXPECT_NEAR(scalar_func->getValue(p), 42.0, 1e-10);
  }

  TEST(Rodin_Variational_ScalarFunction, RealFunction_MoveSemantics)
  {
    RealFunction rf1(99.99);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{0.25, 0.75}};
    Point p(polytope, trans, rc);
    
    Real original_value = rf1.getValue(p);
    
    RealFunction rf2(std::move(rf1));
    EXPECT_NEAR(rf2.getValue(p), original_value, 1e-10);
  }

  TEST(Rodin_Variational_ScalarFunction, RealFunction_ConstantValue_MultiplePoints)
  {
    RealFunction rf(123.456);
    
    // Create multiple points to test that the constant value is returned everywhere
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mesh.getDimension(), 0);
    
    std::vector<Math::Vector<Real>> test_coords = {
      {{0.0, 0.0}},
      {{1.0, 0.0}},
      {{0.0, 1.0}},
      {{0.5, 0.5}},
      {{0.33, 0.67}}
    };
    
    for (const auto& rc : test_coords)
    {
      Point p(polytope, trans, rc);
      EXPECT_NEAR(rf.getValue(p), 123.456, 1e-10);
    }
  }

  TEST(Rodin_Variational_ScalarFunction, RealFunction_IntegerConstant)
  {
    RealFunction rf(789);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, trans, rc);
    
    EXPECT_NEAR(rf.getValue(p), 789.0, 1e-10);
  }

  TEST(Rodin_Variational_ScalarFunction, RealFunction_ZeroValue)
  {
    RealFunction rf(0.0);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, trans, rc);
    
    EXPECT_NEAR(rf.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_ScalarFunction, RealFunction_NegativeValue)
  {
    RealFunction rf(-456.789);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, trans, rc);
    
    EXPECT_NEAR(rf.getValue(p), -456.789, 1e-10);
  }
}