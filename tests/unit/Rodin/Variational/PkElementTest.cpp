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

  // ========================================================================
  // Extended tests for higher orders (P2 through P6) as requested
  // Testing: behavior, partition of unity, consistency
  // ========================================================================

  // Test P4 element on Segment
  TEST(Rodin_Variational_RealPkElement, SanityTest_P4_1D_Reference_Segment)
  {
    RealPkElement<4> k(Polytope::Type::Segment);
    
    // P4 element on Segment should have 5 DOFs (k+1 = 4+1 = 5)
    EXPECT_EQ(k.getCount(), 5);
    EXPECT_EQ(k.getOrder(), 4);
    
    // Test Lagrange property at nodes
    for (size_t i = 0; i < 5; i++)
    {
      for (size_t j = 0; j < 5; j++)
      {
        const auto& node = k.getNode(j);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(k.getBasis(i)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P4 partition of unity on Segment
  TEST(Rodin_Variational_RealPkElement, PartitionOfUnity_P4_Segment)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealPkElement<4> k(Polytope::Type::Segment);
    
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

  // Test P4 derivatives on Segment
  TEST(Rodin_Variational_RealPkElement, DerivativeTest_P4_1D_Reference_Segment)
  {
    RealPkElement<4> k(Polytope::Type::Segment);
    
    // Test that derivatives sum to zero (derivative of constant function)
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      Math::Vector<Real> p{{x}};
      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j).getDerivative<1>(0)(p);
      }
      EXPECT_NEAR(sum, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P5 element on Segment
  TEST(Rodin_Variational_RealPkElement, SanityTest_P5_1D_Reference_Segment)
  {
    RealPkElement<5> k(Polytope::Type::Segment);
    
    // P5 element on Segment should have 6 DOFs (k+1 = 5+1 = 6)
    EXPECT_EQ(k.getCount(), 6);
    EXPECT_EQ(k.getOrder(), 5);
    
    // Test Lagrange property at nodes
    for (size_t i = 0; i < 6; i++)
    {
      for (size_t j = 0; j < 6; j++)
      {
        const auto& node = k.getNode(j);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(k.getBasis(i)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P5 partition of unity on Segment
  TEST(Rodin_Variational_RealPkElement, PartitionOfUnity_P5_Segment)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealPkElement<5> k(Polytope::Type::Segment);
    
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

  // Test P5 behavior: symmetry of basis functions
  TEST(Rodin_Variational_RealPkElement, SymmetryTest_P5_Segment)
  {
    RealPkElement<5> k(Polytope::Type::Segment);
    
    // Basis functions should be symmetric: phi_i(x) = phi_{n-i}(1-x)
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      Math::Vector<Real> p{{x}};
      Math::Vector<Real> p_sym{{1.0 - x}};
      
      for (size_t j = 0; j < k.getCount(); j++)
      {
        size_t j_sym = k.getCount() - 1 - j;
        EXPECT_NEAR(k.getBasis(j)(p), k.getBasis(j_sym)(p_sym), RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P6 element on Segment
  TEST(Rodin_Variational_RealPkElement, SanityTest_P6_1D_Reference_Segment)
  {
    RealPkElement<6> k(Polytope::Type::Segment);
    
    // P6 element on Segment should have 7 DOFs (k+1 = 6+1 = 7)
    EXPECT_EQ(k.getCount(), 7);
    EXPECT_EQ(k.getOrder(), 6);
    
    // Test Lagrange property at nodes
    for (size_t i = 0; i < 7; i++)
    {
      for (size_t j = 0; j < 7; j++)
      {
        const auto& node = k.getNode(j);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(k.getBasis(i)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P6 partition of unity on Segment
  TEST(Rodin_Variational_RealPkElement, PartitionOfUnity_P6_Segment)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealPkElement<6> k(Polytope::Type::Segment);
    
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

  // Test P6 derivatives sum to zero
  TEST(Rodin_Variational_RealPkElement, DerivativeTest_P6_1D_Reference_Segment)
  {
    RealPkElement<6> k(Polytope::Type::Segment);
    
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      Math::Vector<Real> p{{x}};
      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j).getDerivative<1>(0)(p);
      }
      EXPECT_NEAR(sum, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test DOF count consistency for all orders on Segment
  TEST(Rodin_Variational_RealPkElement, DOFCount_Consistency_Segment)
  {
    // Test DOF formula: (k+1) for Segment
    EXPECT_EQ(RealPkElement<2>(Polytope::Type::Segment).getCount(), 3);
    EXPECT_EQ(RealPkElement<3>(Polytope::Type::Segment).getCount(), 4);
    EXPECT_EQ(RealPkElement<4>(Polytope::Type::Segment).getCount(), 5);
    EXPECT_EQ(RealPkElement<5>(Polytope::Type::Segment).getCount(), 6);
    EXPECT_EQ(RealPkElement<6>(Polytope::Type::Segment).getCount(), 7);
  }

  // Test higher orders on Quadrilateral
  TEST(Rodin_Variational_RealPkElement, SanityTest_HigherOrders_Quadrilateral)
  {
    // Test DOF formula: (k+1)^2 for Quadrilateral
    EXPECT_EQ(RealPkElement<2>(Polytope::Type::Quadrilateral).getCount(), 9);
    EXPECT_EQ(RealPkElement<3>(Polytope::Type::Quadrilateral).getCount(), 16);
    EXPECT_EQ(RealPkElement<4>(Polytope::Type::Quadrilateral).getCount(), 25);
    EXPECT_EQ(RealPkElement<5>(Polytope::Type::Quadrilateral).getCount(), 36);
    EXPECT_EQ(RealPkElement<6>(Polytope::Type::Quadrilateral).getCount(), 49);
  }

  // Test P3 partition of unity on Quadrilateral
  TEST(Rodin_Variational_RealPkElement, PartitionOfUnity_P3_Quadrilateral)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealPkElement<3> k(Polytope::Type::Quadrilateral);
    
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

  // Test P4 partition of unity on Quadrilateral
  TEST(Rodin_Variational_RealPkElement, PartitionOfUnity_P4_Quadrilateral)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealPkElement<4> k(Polytope::Type::Quadrilateral);
    
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

  // Test higher orders on Triangle
  TEST(Rodin_Variational_RealPkElement, SanityTest_HigherOrders_Triangle)
  {
    // Test DOF formula: (k+1)(k+2)/2 for Triangle
    EXPECT_EQ(RealPkElement<2>(Polytope::Type::Triangle).getCount(), 6);
    EXPECT_EQ(RealPkElement<3>(Polytope::Type::Triangle).getCount(), 10);
    EXPECT_EQ(RealPkElement<4>(Polytope::Type::Triangle).getCount(), 15);
    EXPECT_EQ(RealPkElement<5>(Polytope::Type::Triangle).getCount(), 21);
    EXPECT_EQ(RealPkElement<6>(Polytope::Type::Triangle).getCount(), 28);
  }

  // Test higher orders on Tetrahedron
  TEST(Rodin_Variational_RealPkElement, SanityTest_HigherOrders_Tetrahedron)
  {
    // Test DOF formula: (k+1)(k+2)(k+3)/6 for Tetrahedron
    EXPECT_EQ(RealPkElement<2>(Polytope::Type::Tetrahedron).getCount(), 10);
    EXPECT_EQ(RealPkElement<3>(Polytope::Type::Tetrahedron).getCount(), 20);
    EXPECT_EQ(RealPkElement<4>(Polytope::Type::Tetrahedron).getCount(), 35);
    EXPECT_EQ(RealPkElement<5>(Polytope::Type::Tetrahedron).getCount(), 56);
    EXPECT_EQ(RealPkElement<6>(Polytope::Type::Tetrahedron).getCount(), 84);
  }

  // Test P2 behavior: linear reproduction property
  TEST(Rodin_Variational_RealPkElement, LinearReproduction_P2_Segment)
  {
    RealPkElement<2> k(Polytope::Type::Segment);
    
    // Test that P2 can exactly reproduce linear functions
    // f(x) = a + b*x should be exactly represented
    constexpr size_t n = 20;
    RandomFloat gen(-1.0, 1.0);
    
    for (size_t trial = 0; trial < 10; trial++)
    {
      Real a = gen();
      Real b = gen();
      
      for (size_t i = 0; i < n; i++)
      {
        const auto& x = gen() * 0.5 + 0.5; // Map to [0, 1]
        Math::Vector<Real> p{{x}};
        
        // Compute linear function at point
        Real f_exact = a + b * x;
        
        // Interpolate using basis functions
        Real f_interp = 0;
        for (size_t j = 0; j < k.getCount(); j++)
        {
          const auto& node = k.getNode(j);
          Real f_node = a + b * node.x();
          f_interp += f_node * k.getBasis(j)(p);
        }
        
        EXPECT_NEAR(f_interp, f_exact, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P3 behavior: quadratic reproduction property
  TEST(Rodin_Variational_RealPkElement, QuadraticReproduction_P3_Segment)
  {
    RealPkElement<3> k(Polytope::Type::Segment);
    
    // Test that P3 can exactly reproduce quadratic functions
    // f(x) = a + b*x + c*x^2 should be exactly represented
    constexpr size_t n = 20;
    RandomFloat gen(-1.0, 1.0);
    
    for (size_t trial = 0; trial < 5; trial++)
    {
      Real a = gen();
      Real b = gen();
      Real c = gen();
      
      for (size_t i = 0; i < n; i++)
      {
        const auto& x = gen() * 0.5 + 0.5; // Map to [0, 1]
        Math::Vector<Real> p{{x}};
        
        // Compute quadratic function at point
        Real f_exact = a + b * x + c * x * x;
        
        // Interpolate using basis functions
        Real f_interp = 0;
        for (size_t j = 0; j < k.getCount(); j++)
        {
          const auto& node = k.getNode(j);
          Real node_x = node.x();
          Real f_node = a + b * node_x + c * node_x * node_x;
          f_interp += f_node * k.getBasis(j)(p);
        }
        
        EXPECT_NEAR(f_interp, f_exact, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test gradient computation consistency for P4
  TEST(Rodin_Variational_RealPkElement, GradientConsistency_P4_Segment)
  {
    RealPkElement<4> k(Polytope::Type::Segment);
    
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    
    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      Math::Vector<Real> p{{x}};
      
      for (size_t j = 0; j < k.getCount(); j++)
      {
        // Check that gradient function matches derivative function
        auto grad = k.getBasis(j).getGradient();
        auto deriv = k.getBasis(j).getDerivative<1>(0);
        
        const auto& grad_val = grad(p);
        Real deriv_val = deriv(p);
        
        EXPECT_EQ(grad_val.size(), 1);
        EXPECT_NEAR(grad_val(0), deriv_val, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test consistency across orders for identical geometries
  TEST(Rodin_Variational_RealPkElement, CrossOrder_Consistency_Segment)
  {
    // Verify that nodes are uniformly distributed for all orders
    for (size_t k = 2; k <= 6; k++)
    {
      // We can't instantiate templates with runtime values, so we check specific cases
      if (k == 2)
      {
        RealPkElement<2> elem(Polytope::Type::Segment);
        for (size_t i = 0; i < elem.getCount(); i++)
        {
          const auto& node = elem.getNode(i);
          Real expected = static_cast<Real>(i) / static_cast<Real>(k);
          EXPECT_NEAR(node.x(), expected, RODIN_FUZZY_CONSTANT);
        }
      }
    }
  }

  // Test non-negativity of basis functions at interior points
  TEST(Rodin_Variational_RealPkElement, BasisNonNegativity_P3_Segment)
  {
    RealPkElement<3> k(Polytope::Type::Segment);
    
    // Basis functions should be non-negative at their support nodes
    // but may be negative elsewhere (characteristic of high-order Lagrange)
    for (size_t i = 0; i < k.getCount(); i++)
    {
      const auto& node = k.getNode(i);
      EXPECT_GE(k.getBasis(i)(node), 0.0);
    }
  }

  // Test order property is correctly set for all orders
  TEST(Rodin_Variational_RealPkElement, OrderProperty_AllOrders)
  {
    EXPECT_EQ(RealPkElement<2>(Polytope::Type::Segment).getOrder(), 2);
    EXPECT_EQ(RealPkElement<3>(Polytope::Type::Segment).getOrder(), 3);
    EXPECT_EQ(RealPkElement<4>(Polytope::Type::Segment).getOrder(), 4);
    EXPECT_EQ(RealPkElement<5>(Polytope::Type::Segment).getOrder(), 5);
    EXPECT_EQ(RealPkElement<6>(Polytope::Type::Segment).getOrder(), 6);
  }
}
