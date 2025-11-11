#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational/P1.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_RealP1Element, SanityTest_1D_Reference_Segment)
  {
    RealP1Element k(Polytope::Type::Segment);

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0}}), 0, RODIN_FUZZY_CONSTANT);
    }

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
    }

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{1}}), 1, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_RealP1Element, FuzzyTest_1D_Reference_Segment)
  {
    constexpr size_t n = 25;

    RandomFloat gen(0.0, 1.0);
    RealP1Element k(Polytope::Type::Segment);

    {
      for (size_t i = 0; i < n; i++)
      {
        const auto& s = gen();
        Math::Vector<Real> p = (1 - s) * Math::Vector<Real>{{1}};
        EXPECT_NEAR(k.getBasis(0)(p), 1 - p.x(), RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR(k.getBasis(1)(p), p.x(), RODIN_FUZZY_CONSTANT);
      }
    }
  }

  TEST(Rodin_Variational_RealP1Element, SanityTest_2D_Reference_Triangle)
  {
    RealP1Element k(Polytope::Type::Triangle);

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{1, 0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{1, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 1}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0, 1}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0, 1}}), 1, RODIN_FUZZY_CONSTANT);
    }

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.5, 0}}), 0.5, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0.5, 0}}), 0.5, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0.5, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.5, 0.5}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0.5, 0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0.5, 0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
    }

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0, 0.5}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0, 0.5}}), 0.5, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_RealP1Element, FuzzyTest_2D_Reference_Triangle)
  {
    constexpr size_t n = 25;

    RandomFloat gen(0.0, 1.0);
    RealP1Element k(Polytope::Type::Triangle);

    {
      for (size_t i = 0; i < n; i++)
      {
        const auto& s = gen();
        Math::Vector<Real> p = s * Math::Vector<Real>{{0, 0}} + (1 - s) * Math::Vector<Real>{{1, 0}};
        EXPECT_NEAR(k.getBasis(0)(p), -p.x() - p.y() + 1, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR(k.getBasis(1)(p), p.x(), RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR(k.getBasis(2)(p), p.y(), RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P1 element properties: DOF count, order, nodes
  TEST(Rodin_Variational_RealP1Element, ElementProperties_Segment)
  {
    RealP1Element k(Polytope::Type::Segment);
    
    // Check DOF count (should be number of vertices = 2)
    EXPECT_EQ(k.getCount(), 2);
    
    // Check order (linear = order 1)
    EXPECT_EQ(k.getOrder(), 1);
    
    // Check nodes (vertices of reference segment)
    const auto& node0 = k.getNode(0);
    const auto& node1 = k.getNode(1);
    EXPECT_NEAR(node0.x(), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node1.x(), 1, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Variational_RealP1Element, ElementProperties_Triangle)
  {
    RealP1Element k(Polytope::Type::Triangle);
    
    // Check DOF count (should be number of vertices = 3)
    EXPECT_EQ(k.getCount(), 3);
    
    // Check order (linear = order 1)
    EXPECT_EQ(k.getOrder(), 1);
    
    // Check nodes (vertices of reference triangle)
    const auto& node0 = k.getNode(0);
    const auto& node1 = k.getNode(1);
    const auto& node2 = k.getNode(2);
    EXPECT_NEAR(node0.x(), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node0.y(), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node1.x(), 1, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node1.y(), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node2.x(), 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node2.y(), 1, RODIN_FUZZY_CONSTANT);
  }

  // Test P1 element derivatives on Segment
  TEST(Rodin_Variational_RealP1Element, DerivativeTest_1D_Reference_Segment)
  {
    RealP1Element k(Polytope::Type::Segment);
    
    // First basis function derivative
    {
      auto deriv = k.getBasis(0).getDerivative<1>(0);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0}}), -1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0.5}}), -1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{1}}), -1, RODIN_FUZZY_CONSTANT);
    }
    
    // Second basis function derivative
    {
      auto deriv = k.getBasis(1).getDerivative<1>(0);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0.5}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{1}}), 1, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P1 element derivatives on Triangle
  TEST(Rodin_Variational_RealP1Element, DerivativeTest_2D_Reference_Triangle)
  {
    RealP1Element k(Polytope::Type::Triangle);
    
    Math::Vector<Real> p{{0.25, 0.25}};
    
    // First basis function: phi_0 = 1 - x - y
    {
      auto deriv_x = k.getBasis(0).getDerivative<1>(0);
      auto deriv_y = k.getBasis(0).getDerivative<1>(1);
      EXPECT_NEAR(deriv_x(p), -1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv_y(p), -1, RODIN_FUZZY_CONSTANT);
    }
    
    // Second basis function: phi_1 = x
    {
      auto deriv_x = k.getBasis(1).getDerivative<1>(0);
      auto deriv_y = k.getBasis(1).getDerivative<1>(1);
      EXPECT_NEAR(deriv_x(p), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv_y(p), 0, RODIN_FUZZY_CONSTANT);
    }
    
    // Third basis function: phi_2 = y
    {
      auto deriv_x = k.getBasis(2).getDerivative<1>(0);
      auto deriv_y = k.getBasis(2).getDerivative<1>(1);
      EXPECT_NEAR(deriv_x(p), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv_y(p), 1, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test gradient function on Triangle
  TEST(Rodin_Variational_RealP1Element, GradientTest_2D_Reference_Triangle)
  {
    RealP1Element k(Polytope::Type::Triangle);
    
    Math::Vector<Real> p{{0.25, 0.25}};
    
    // First basis function gradient: grad(phi_0) = [-1, -1]
    {
      auto grad = k.getBasis(0).getGradient();
      const auto& grad_val = grad(p);
      EXPECT_EQ(grad_val.size(), 2);
      EXPECT_NEAR(grad_val(0), -1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_val(1), -1, RODIN_FUZZY_CONSTANT);
    }
    
    // Second basis function gradient: grad(phi_1) = [1, 0]
    {
      auto grad = k.getBasis(1).getGradient();
      const auto& grad_val = grad(p);
      EXPECT_EQ(grad_val.size(), 2);
      EXPECT_NEAR(grad_val(0), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_val(1), 0, RODIN_FUZZY_CONSTANT);
    }
    
    // Third basis function gradient: grad(phi_2) = [0, 1]
    {
      auto grad = k.getBasis(2).getGradient();
      const auto& grad_val = grad(p);
      EXPECT_EQ(grad_val.size(), 2);
      EXPECT_NEAR(grad_val(0), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(grad_val(1), 1, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test partition of unity property: sum of basis functions = 1
  TEST(Rodin_Variational_RealP1Element, PartitionOfUnity_Segment)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealP1Element k(Polytope::Type::Segment);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      Math::Vector<Real> p{{s}};
      Real sum = k.getBasis(0)(p) + k.getBasis(1)(p);
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_RealP1Element, PartitionOfUnity_Triangle)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealP1Element k(Polytope::Type::Triangle);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      const auto& t = gen();
      if (s + t <= 1.0)
      {
        Math::Vector<Real> p{{s, t}};
        Real sum = k.getBasis(0)(p) + k.getBasis(1)(p) + k.getBasis(2)(p);
        EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P1 element on Tetrahedron
  TEST(Rodin_Variational_RealP1Element, SanityTest_3D_Reference_Tetrahedron)
  {
    RealP1Element k(Polytope::Type::Tetrahedron);
    
    // Check DOF count (should be number of vertices = 4)
    EXPECT_EQ(k.getCount(), 4);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 1);
    
    // Test at vertices
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0, 0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(3)(Math::Vector<Real>{{0, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{1, 0, 0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{1, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(3)(Math::Vector<Real>{{1, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 1, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0, 1, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0, 1, 0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(3)(Math::Vector<Real>{{0, 1, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0, 1}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0, 0, 1}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0, 0, 1}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(3)(Math::Vector<Real>{{0, 0, 1}}), 1, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test partition of unity on Tetrahedron
  TEST(Rodin_Variational_RealP1Element, PartitionOfUnity_Tetrahedron)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealP1Element k(Polytope::Type::Tetrahedron);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      const auto& t = gen();
      const auto& u = gen();
      if (s + t + u <= 1.0)
      {
        Math::Vector<Real> p{{s, t, u}};
        Real sum = k.getBasis(0)(p) + k.getBasis(1)(p) + k.getBasis(2)(p) + k.getBasis(3)(p);
        EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P1 element on Quadrilateral
  TEST(Rodin_Variational_RealP1Element, SanityTest_2D_Reference_Quadrilateral)
  {
    RealP1Element k(Polytope::Type::Quadrilateral);
    
    // Check DOF count (should be number of vertices = 4)
    EXPECT_EQ(k.getCount(), 4);
    
    // Check order (bilinear = order 2 for quadrilateral)
    EXPECT_EQ(k.getOrder(), 2);
    
    // Test at vertices
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(3)(Math::Vector<Real>{{0, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{1, 0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{1, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(3)(Math::Vector<Real>{{1, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test partition of unity on Quadrilateral
  TEST(Rodin_Variational_RealP1Element, PartitionOfUnity_Quadrilateral)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealP1Element k(Polytope::Type::Quadrilateral);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      const auto& t = gen();
      Math::Vector<Real> p{{s, t}};
      Real sum = k.getBasis(0)(p) + k.getBasis(1)(p) + k.getBasis(2)(p) + k.getBasis(3)(p);
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }
}
