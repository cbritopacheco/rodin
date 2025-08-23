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
  TEST(Rodin_Variational_RealFunction, ConstantReal_Construction)
  {
    RealFunction f(3.14);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    EXPECT_NEAR(f.getValue(p), 3.14, 1e-10);
  }

  TEST(Rodin_Variational_RealFunction, ConstantInteger_Construction)
  {
    RealFunction f(42);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    EXPECT_NEAR(f.getValue(p), 42.0, 1e-10);
  }

  TEST(Rodin_Variational_RealFunction, CopyConstructor_Real)
  {
    RealFunction f(2.718);
    RealFunction f_copy(f);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    EXPECT_NEAR(f_copy.getValue(p), f.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_RealFunction, MoveConstructor_Real)
  {
    RealFunction f(1.414);
    Real original_value =
      f.getValue(
          Point(
            *LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 }).getPolytope(2, 0),
            Math::Vector<Real>{{0.0, 0.0}}
      )
    );
    RealFunction f_moved(std::move(f));
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    EXPECT_NEAR(f_moved.getValue(p), original_value, 1e-10);
  }

  TEST(Rodin_Variational_RealFunction, Copy_Real)
  {
    RealFunction f(123.456);
    auto copied = f.copy();
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    EXPECT_NE(copied, nullptr);
    EXPECT_NEAR(copied->getValue(p), f.getValue(p), 1e-10);
    delete copied;
  }

  TEST(Rodin_Variational_RealFunction, Copy_Integer)
  {
    RealFunction f(789);
    auto copied = f.copy();
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    EXPECT_NE(copied, nullptr);
    EXPECT_NEAR(copied->getValue(p), 789.0, 1e-10);
    delete copied;
  }

  TEST(Rodin_Variational_RealFunction, ZeroValue)
  {
    RealFunction f(0.0);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 1.0}};
    Point p(polytope, rc);
    EXPECT_NEAR(f.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_RealFunction, NegativeValue)
  {
    RealFunction f(-99.99);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);
    EXPECT_NEAR(f.getValue(p), -99.99, 1e-10);
  }

  TEST(Rodin_Variational_RealFunction, ConstantFunction_TracedToBoundary)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D - 1, D);
    RealFunction f(5.0);
    auto traced_f = f.traceOf(RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
    // Create a boundary point for testing
    auto it = mesh.getPolytope(D - 1, 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5}};
    Point p(polytope, rc);
    EXPECT_NEAR(traced_f.getValue(p), 5.0, 1e-10);
  }
}
