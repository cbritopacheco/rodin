#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

// Disable strict aliasing and array-bounds warnings from Eigen for tests
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif

#include "Rodin/Variational/Pk/PkElement.h"
#include "Rodin/Variational/P0/P0Element.h"
#include "Rodin/Variational/P1/P1Element.h"

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  // Test P2 element (K=2) on Point geometry
  TEST(Rodin_Variational_RealPkElement, SanityTest_P2_0D_Reference_Point)
  {
    RealPkElement<2> k(Polytope::Type::Point);
    
    // P2 element on Point should have 1 DOF
    EXPECT_EQ(k.getCount(), 1);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 2);
    
    // Check basis function value
    EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0}}), 1, RODIN_FUZZY_CONSTANT);
  }

  // Test P2 element (K=2) on Segment geometry
  TEST(Rodin_Variational_RealPkElement, SanityTest_P2_1D_Reference_Segment)
  {
    RealPkElement<2> k(Polytope::Type::Segment);
    
    // P2 element on Segment should have 3 DOFs (k+1 = 2+1 = 3)
    EXPECT_EQ(k.getCount(), 3);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 2);
    
    // Test basis functions at nodes (Lagrange property)
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0.0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0.0}}), 0, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.5}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0.5}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0.5}}), 0, RODIN_FUZZY_CONSTANT);
    }
    
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1.0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{1.0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{1.0}}), 1, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 element partition of unity on Segment
  TEST(Rodin_Variational_RealPkElement, PartitionOfUnity_P2_Segment)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealPkElement<2> k(Polytope::Type::Segment);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      Math::Vector<Real> p{{s}};
      Real sum = k.getBasis(0)(p) + k.getBasis(1)(p) + k.getBasis(2)(p);
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 element derivatives on Segment
  TEST(Rodin_Variational_RealPkElement, DerivativeTest_P2_1D_Reference_Segment)
  {
    RealPkElement<2> k(Polytope::Type::Segment);
    
    // First basis function: phi_0(x) = (1-x)(1-2x) = 2x^2 - 3x + 1
    // phi_0'(x) = 4x - 3
    {
      auto deriv = k.getBasis(0).getDerivative<1>(0);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0.0}}), -3, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0.5}}), -1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{1.0}}), 1, RODIN_FUZZY_CONSTANT);
    }
    
    // Second basis function: phi_1(x) = 4x(1-x) = -4x^2 + 4x
    // phi_1'(x) = -8x + 4
    {
      auto deriv = k.getBasis(1).getDerivative<1>(0);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0.0}}), 4, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0.5}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{1.0}}), -4, RODIN_FUZZY_CONSTANT);
    }
    
    // Third basis function: phi_2(x) = x(2x-1) = 2x^2 - x
    // phi_2'(x) = 4x - 1
    {
      auto deriv = k.getBasis(2).getDerivative<1>(0);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0.0}}), -1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0.5}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{1.0}}), 3, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 element on Triangle geometry
  TEST(Rodin_Variational_RealPkElement, SanityTest_P2_2D_Reference_Triangle)
  {
    RealPkElement<2> k(Polytope::Type::Triangle);
    
    // P2 element on Triangle should have 6 DOFs: (k+1)(k+2)/2 = 3*4/2 = 6
    EXPECT_EQ(k.getCount(), 6);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 2);
  }

  // Test P2 element on Quadrilateral geometry
  TEST(Rodin_Variational_RealPkElement, SanityTest_P2_2D_Reference_Quadrilateral)
  {
    RealPkElement<2> k(Polytope::Type::Quadrilateral);
    
    // P2 element on Quadrilateral should have 9 DOFs: (k+1)^2 = 3^2 = 9
    EXPECT_EQ(k.getCount(), 9);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 2);
    
    // Test tensor product property at corner
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.0, 0.0}}), 1, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 element partition of unity on Quadrilateral
  TEST(Rodin_Variational_RealPkElement, PartitionOfUnity_P2_Quadrilateral)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealPkElement<2> k(Polytope::Type::Quadrilateral);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      const auto& t = gen();
      Math::Vector<Real> p{{s, t}};
      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test that P1 element matches PkElement<1>
  TEST(Rodin_Variational_RealPkElement, P1_Consistency_Segment)
  {
    RealPkElement<1> pk(Polytope::Type::Segment);
    RealP1Element p1(Polytope::Type::Segment);
    
    // Same number of DOFs
    EXPECT_EQ(pk.getCount(), p1.getCount());
    
    // Same order
    EXPECT_EQ(pk.getOrder(), p1.getOrder());
    
    // Same basis function values at several points
    constexpr size_t n = 10;
    RandomFloat gen(0.0, 1.0);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      Math::Vector<Real> p{{x}};
      for (size_t j = 0; j < pk.getCount(); j++)
      {
        EXPECT_NEAR(pk.getBasis(j)(p), p1.getBasis(j)(p), RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test that P0 element matches PkElement<0>
  TEST(Rodin_Variational_RealPkElement, P0_Consistency_Segment)
  {
    RealPkElement<0> pk(Polytope::Type::Segment);
    RealP0Element p0(Polytope::Type::Segment);
    
    // Same number of DOFs
    EXPECT_EQ(pk.getCount(), p0.getCount());
    
    // Same order
    EXPECT_EQ(pk.getOrder(), p0.getOrder());
    
    // Same basis function value (constant = 1)
    constexpr size_t n = 10;
    RandomFloat gen(0.0, 1.0);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      Math::Vector<Real> p{{x}};
      EXPECT_NEAR(pk.getBasis(0)(p), p0.getBasis(0)(p), RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P3 element on Segment
  TEST(Rodin_Variational_RealPkElement, SanityTest_P3_1D_Reference_Segment)
  {
    RealPkElement<3> k(Polytope::Type::Segment);
    
    // P3 element on Segment should have 4 DOFs (k+1 = 3+1 = 4)
    EXPECT_EQ(k.getCount(), 4);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 3);
    
    // Test Lagrange property at nodes
    for (size_t i = 0; i < 4; i++)
    {
      for (size_t j = 0; j < 4; j++)
      {
        const auto& node = k.getNode(j);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(k.getBasis(i)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P3 element partition of unity on Segment
  TEST(Rodin_Variational_RealPkElement, PartitionOfUnity_P3_Segment)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealPkElement<3> k(Polytope::Type::Segment);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      Math::Vector<Real> p{{s}};
      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 element on Tetrahedron geometry
  TEST(Rodin_Variational_RealPkElement, SanityTest_P2_3D_Reference_Tetrahedron)
  {
    RealPkElement<2> k(Polytope::Type::Tetrahedron);
    
    // P2 element on Tetrahedron should have 10 DOFs: (k+1)(k+2)(k+3)/6 = 3*4*5/6 = 10
    EXPECT_EQ(k.getCount(), 10);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 2);
  }

  // Test P2 element on Wedge geometry
  TEST(Rodin_Variational_RealPkElement, SanityTest_P2_3D_Reference_Wedge)
  {
    RealPkElement<2> k(Polytope::Type::Wedge);
    
    // P2 element on Wedge: (k+1) * (k+1)(k+2)/2 = 3 * 3*4/2 = 3 * 6 = 18
    EXPECT_EQ(k.getCount(), 18);
    
    // Check order
    EXPECT_EQ(k.getOrder(), 2);
  }
}
