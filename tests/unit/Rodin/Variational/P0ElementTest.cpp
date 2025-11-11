#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational/P0/P0Element.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  // Test P0 element on Point geometry
  TEST(Rodin_Variational_RealP0Element, SanityTest_0D_Reference_Point)
  {
    RealP0Element k(Polytope::Type::Point);
    
    // P0 element should always return 1 for basis function evaluation
    EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0}}), 1, RODIN_FUZZY_CONSTANT);
    
    // Check DOF count
    EXPECT_EQ(k.getCount(), 1);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 0);
  }

  // Test P0 element on Segment geometry
  TEST(Rodin_Variational_RealP0Element, SanityTest_1D_Reference_Segment)
  {
    RealP0Element k(Polytope::Type::Segment);
    
    // P0 element basis function should always return 1
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.5}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    // Check that node is at barycenter
    const auto& node = k.getNode(0);
    EXPECT_NEAR(node.x(), 0.5, RODIN_FUZZY_CONSTANT);
    
    // Check DOF count
    EXPECT_EQ(k.getCount(), 1);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 0);
  }

  // Test P0 element derivatives on Segment (should all be zero)
  TEST(Rodin_Variational_RealP0Element, DerivativeTest_1D_Reference_Segment)
  {
    RealP0Element k(Polytope::Type::Segment);
    
    // Derivatives of constant functions should be zero
    {
      auto deriv = k.getBasis(0).getDerivative<1>(0);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0.5}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{1}}), 0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P0 element on Triangle geometry
  TEST(Rodin_Variational_RealP0Element, SanityTest_2D_Reference_Triangle)
  {
    RealP0Element k(Polytope::Type::Triangle);
    
    // P0 element basis function should always return 1
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 1}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.5, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.5, 0.5}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0.5}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    // Check that node is at barycenter (1/3, 1/3)
    const auto& node = k.getNode(0);
    EXPECT_NEAR(node.x(), Real(1) / Real(3), RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node.y(), Real(1) / Real(3), RODIN_FUZZY_CONSTANT);
    
    // Check DOF count
    EXPECT_EQ(k.getCount(), 1);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 0);
  }

  // Test P0 element derivatives on Triangle (should all be zero)
  TEST(Rodin_Variational_RealP0Element, DerivativeTest_2D_Reference_Triangle)
  {
    RealP0Element k(Polytope::Type::Triangle);
    
    // Derivatives of constant functions should be zero
    {
      auto deriv_x = k.getBasis(0).getDerivative<1>(0);
      auto deriv_y = k.getBasis(0).getDerivative<1>(1);
      
      Math::Vector<Real> p{{0.25, 0.25}};
      EXPECT_NEAR(deriv_x(p), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv_y(p), 0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P0 element on Quadrilateral geometry
  TEST(Rodin_Variational_RealP0Element, SanityTest_2D_Reference_Quadrilateral)
  {
    RealP0Element k(Polytope::Type::Quadrilateral);
    
    // P0 element basis function should always return 1
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1, 1}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 1}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.5, 0.5}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    // Check that node is at barycenter (0.5, 0.5)
    const auto& node = k.getNode(0);
    EXPECT_NEAR(node.x(), 0.5, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node.y(), 0.5, RODIN_FUZZY_CONSTANT);
    
    // Check DOF count
    EXPECT_EQ(k.getCount(), 1);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 0);
  }

  // Test P0 element on Tetrahedron geometry
  TEST(Rodin_Variational_RealP0Element, SanityTest_3D_Reference_Tetrahedron)
  {
    RealP0Element k(Polytope::Type::Tetrahedron);
    
    // P0 element basis function should always return 1
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1, 0, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 1, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0, 1}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    // Check that node is at barycenter (0.25, 0.25, 0.25)
    const auto& node = k.getNode(0);
    EXPECT_NEAR(node.x(), 0.25, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node.y(), 0.25, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node.z(), 0.25, RODIN_FUZZY_CONSTANT);
    
    // Check DOF count
    EXPECT_EQ(k.getCount(), 1);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 0);
  }

  // Test P0 element derivatives on Tetrahedron (should all be zero)
  TEST(Rodin_Variational_RealP0Element, DerivativeTest_3D_Reference_Tetrahedron)
  {
    RealP0Element k(Polytope::Type::Tetrahedron);
    
    // Derivatives of constant functions should be zero
    {
      auto deriv_x = k.getBasis(0).getDerivative<1>(0);
      auto deriv_y = k.getBasis(0).getDerivative<1>(1);
      auto deriv_z = k.getBasis(0).getDerivative<1>(2);
      
      Math::Vector<Real> p{{0.25, 0.25, 0.25}};
      EXPECT_NEAR(deriv_x(p), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv_y(p), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv_z(p), 0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Fuzzy test for P0 element on different geometries
  TEST(Rodin_Variational_RealP0Element, FuzzyTest_AllGeometries)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    
    // Test on Segment
    {
      RealP0Element k(Polytope::Type::Segment);
      for (size_t i = 0; i < n; i++)
      {
        const auto& s = gen();
        Math::Vector<Real> p{{s}};
        EXPECT_NEAR(k.getBasis(0)(p), 1, RODIN_FUZZY_CONSTANT);
      }
    }
    
    // Test on Triangle
    {
      RealP0Element k(Polytope::Type::Triangle);
      for (size_t i = 0; i < n; i++)
      {
        const auto& s = gen();
        const auto& t = gen();
        if (s + t <= 1.0)
        {
          Math::Vector<Real> p{{s, t}};
          EXPECT_NEAR(k.getBasis(0)(p), 1, RODIN_FUZZY_CONSTANT);
        }
      }
    }
  }

  // Test P0 element on Wedge geometry
  TEST(Rodin_Variational_RealP0Element, SanityTest_3D_Reference_Wedge)
  {
    RealP0Element k(Polytope::Type::Wedge);
    
    // P0 element basis function should always return 1
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1, 0, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 1, 0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0, 1}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    // Check that node is at barycenter (1/3, 1/3, 0.5)
    const auto& node = k.getNode(0);
    EXPECT_NEAR(node.x(), Real(1) / Real(3), RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node.y(), Real(1) / Real(3), RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(node.z(), 0.5, RODIN_FUZZY_CONSTANT);
    
    // Check DOF count
    EXPECT_EQ(k.getCount(), 1);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 0);
  }
}
