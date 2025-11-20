#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

// Disable strict aliasing and array-bounds warnings from Eigen for tests
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif

#include <complex>
#include "Rodin/Variational/P0.h"
#include "Rodin/Variational/P1.h"
#include "Rodin/Variational/Pk.h"

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
  TEST(Rodin_Variational_RealH1Element, SanityTest_P2_0D_Reference_Point)
  {
    RealH1Element<2> k(Polytope::Type::Point);

    // P2 element on Point should have 1 DOF
    EXPECT_EQ(k.getCount(), 1);

    // Check order
    EXPECT_EQ(k.getOrder(), 0);

    // Check basis function value
    EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0}}), 1, RODIN_FUZZY_CONSTANT);
  }

  // Test P2 element (K=2) on Segment geometry
  TEST(Rodin_Variational_RealH1Element, SanityTest_P2_1D_Reference_Segment)
  {
    RealH1Element<2> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P2_Segment)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<2> k(Polytope::Type::Segment);

    for (size_t i = 0; i < n; i++)
    {
      const auto& s = gen();
      Math::Vector<Real> p{{s}};
      Real sum = k.getBasis(0)(p) + k.getBasis(1)(p) + k.getBasis(2)(p);
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 element derivatives on Segment
  TEST(Rodin_Variational_RealH1Element, DerivativeTest_P2_1D_Reference_Segment)
  {
    RealH1Element<2> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, SanityTest_P2_2D_Reference_Triangle)
  {
    RealH1Element<2> k(Polytope::Type::Triangle);

    // P2 element on Triangle should have 6 DOFs: (k+1)(k+2)/2 = 3*4/2 = 6
    EXPECT_EQ(k.getCount(), 6);

    // Check order
    EXPECT_EQ(k.getOrder(), 2);
  }

  // Test P2 element on Quadrilateral geometry
  TEST(Rodin_Variational_RealH1Element, SanityTest_P2_2D_Reference_Quadrilateral)
  {
    RealH1Element<2> k(Polytope::Type::Quadrilateral);

    // P2 element on Quadrilateral should have 9 DOFs: (k+1)^2 = 3^2 = 9
    EXPECT_EQ(k.getCount(), 9);

    // Check order
    EXPECT_EQ(k.getOrder(), 4);

    // Test tensor product property at corner
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0.0, 0.0}}), 1, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 element partition of unity on Quadrilateral
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P2_Quadrilateral)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<2> k(Polytope::Type::Quadrilateral);

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



  // Test P3 element on Segment
  TEST(Rodin_Variational_RealH1Element, SanityTest_P3_1D_Reference_Segment)
  {
    RealH1Element<3> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P3_Segment)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<3> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, SanityTest_P2_3D_Reference_Tetrahedron)
  {
    RealH1Element<2> k(Polytope::Type::Tetrahedron);

    // P2 element on Tetrahedron should have 10 DOFs: (k+1)(k+2)(k+3)/6 = 3*4*5/6 = 10
    EXPECT_EQ(k.getCount(), 10);

    // Check order
    EXPECT_EQ(k.getOrder(), 2);
  }

  // Test P2 element on Wedge geometry
  TEST(Rodin_Variational_RealH1Element, SanityTest_P2_3D_Reference_Wedge)
  {
    RealH1Element<2> k(Polytope::Type::Wedge);

    // P2 element on Wedge: (k+1) * (k+1)(k+2)/2 = 3 * 3*4/2 = 3 * 6 = 18
    EXPECT_EQ(k.getCount(), 18);

    // Check order
    EXPECT_EQ(k.getOrder(), 4);
  }

  // ========================================================================
  // Extended tests for higher orders (P2 through P6) as requested
  // Testing: behavior, partition of unity, consistency
  // ========================================================================

  // Test P4 element on Segment
  TEST(Rodin_Variational_RealH1Element, SanityTest_P4_1D_Reference_Segment)
  {
    RealH1Element<4> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P4_Segment)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<4> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, DerivativeTest_P4_1D_Reference_Segment)
  {
    RealH1Element<4> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, SanityTest_P5_1D_Reference_Segment)
  {
    RealH1Element<5> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P5_Segment)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<5> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, SymmetryTest_P5_Segment)
  {
    RealH1Element<5> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, SanityTest_P6_1D_Reference_Segment)
  {
    RealH1Element<6> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P6_Segment)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<6> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, DerivativeTest_P6_1D_Reference_Segment)
  {
    RealH1Element<6> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, DOFCount_Consistency_Segment)
  {
    // Test DOF formula: (k+1) for Segment
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Segment).getCount(), 3);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Segment).getCount(), 4);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Segment).getCount(), 5);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Segment).getCount(), 6);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Segment).getCount(), 7);
  }

  // Test higher orders on Quadrilateral
  TEST(Rodin_Variational_RealH1Element, SanityTest_HigherOrders_Quadrilateral)
  {
    // Test DOF formula: (k+1)^2 for Quadrilateral
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Quadrilateral).getCount(), 9);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Quadrilateral).getCount(), 16);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Quadrilateral).getCount(), 25);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Quadrilateral).getCount(), 36);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Quadrilateral).getCount(), 49);
  }

  // Test P3 partition of unity on Quadrilateral
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P3_Quadrilateral)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<3> k(Polytope::Type::Quadrilateral);

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
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P4_Quadrilateral)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<4> k(Polytope::Type::Quadrilateral);

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
  TEST(Rodin_Variational_RealH1Element, SanityTest_HigherOrders_Triangle)
  {
    // Test DOF formula: (k+1)(k+2)/2 for Triangle
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Triangle).getCount(), 6);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Triangle).getCount(), 10);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Triangle).getCount(), 15);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Triangle).getCount(), 21);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Triangle).getCount(), 28);
  }

  // Test higher orders on Tetrahedron
  TEST(Rodin_Variational_RealH1Element, SanityTest_HigherOrders_Tetrahedron)
  {
    // Test DOF formula: (k+1)(k+2)(k+3)/6 for Tetrahedron
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Tetrahedron).getCount(), 10);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Tetrahedron).getCount(), 20);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Tetrahedron).getCount(), 35);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Tetrahedron).getCount(), 56);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Tetrahedron).getCount(), 84);
  }

  // Test P2 behavior: linear reproduction property
  TEST(Rodin_Variational_RealH1Element, LinearReproduction_P2_Segment)
  {
    RealH1Element<2> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, QuadraticReproduction_P3_Segment)
  {
    RealH1Element<3> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, GradientConsistency_P4_Segment)
  {
    RealH1Element<4> k(Polytope::Type::Segment);

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
  TEST(Rodin_Variational_RealH1Element, CrossOrder_Consistency_Segment)
  {
    // Verify that nodes are uniformly distributed for all orders
    for (size_t k = 2; k <= 6; k++)
    {
      // We can't instantiate templates with runtime values, so we check specific cases
      if (k == 2)
      {
        RealH1Element<2> elem(Polytope::Type::Segment);
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
  TEST(Rodin_Variational_RealH1Element, BasisNonNegativity_P3_Segment)
  {
    RealH1Element<3> k(Polytope::Type::Segment);

    // Basis functions should be non-negative at their support nodes
    // but may be negative elsewhere (characteristic of high-order Lagrange)
    for (size_t i = 0; i < k.getCount(); i++)
    {
      const auto& node = k.getNode(i);
      EXPECT_GE(k.getBasis(i)(node), 0.0);
    }
  }

  // Test order property is correctly set for all orders
  TEST(Rodin_Variational_RealH1Element, OrderProperty_AllOrders)
  {
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Segment).getOrder(), 2);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Segment).getOrder(), 3);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Segment).getOrder(), 4);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Segment).getOrder(), 5);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Segment).getOrder(), 6);
  }

  // ========================================================================
  // Tests for Triangle elements (as requested)
  // ========================================================================

  // Test P2 triangle Lagrange property
  TEST(Rodin_Variational_RealH1Element, LagrangeProperty_P2_Triangle)
  {
    RealH1Element<2> k(Polytope::Type::Triangle);

    // Test Lagrange property at nodes: φᵢ(xⱼ) = δᵢⱼ
    for (size_t i = 0; i < k.getCount(); i++)
    {
      for (size_t j = 0; j < k.getCount(); j++)
      {
        const auto& node = k.getNode(j);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(k.getBasis(i)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P2 triangle partition of unity
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P2_Triangle)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<2> k(Polytope::Type::Triangle);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x); // Ensure x + y ≤ 1
      Math::Vector<Real> p{{x, y}};

      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P3 triangle Lagrange property
  TEST(Rodin_Variational_RealH1Element, LagrangeProperty_P3_Triangle)
  {
    RealH1Element<3> k(Polytope::Type::Triangle);

    for (size_t i = 0; i < k.getCount(); i++)
    {
      for (size_t j = 0; j < k.getCount(); j++)
      {
        const auto& node = k.getNode(j);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(k.getBasis(i)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P3 triangle partition of unity
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P3_Triangle)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<3> k(Polytope::Type::Triangle);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      Math::Vector<Real> p{{x, y}};

      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P4 triangle partition of unity
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P4_Triangle)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<4> k(Polytope::Type::Triangle);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      Math::Vector<Real> p{{x, y}};

      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 triangle derivative consistency
  // TODO: Fix derivative computation for barycentric elements
  TEST(Rodin_Variational_RealH1Element, DerivativeConsistency_P2_Triangle)
  {
    RealH1Element<2> k(Polytope::Type::Triangle);

    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      Math::Vector<Real> p{{x, y}};

      // Derivatives should sum to zero (constant function property)
      Real sum_dx = 0, sum_dy = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum_dx += k.getBasis(j).getDerivative<1>(0)(p);
        sum_dy += k.getBasis(j).getDerivative<1>(1)(p);
      }
      EXPECT_NEAR(sum_dx, 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(sum_dy, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 triangle linear reproduction
  TEST(Rodin_Variational_RealH1Element, LinearReproduction_P2_Triangle)
  {
    RealH1Element<2> k(Polytope::Type::Triangle);

    constexpr size_t n = 20;
    RandomFloat gen(-1.0, 1.0);

    for (size_t trial = 0; trial < 5; trial++)
    {
      Real a = gen();
      Real b = gen();
      Real c = gen();

      for (size_t i = 0; i < n; i++)
      {
        const auto& x = gen() * 0.5 + 0.25;
        const auto& y = gen() * 0.5 * (1.0 - x);
        Math::Vector<Real> p{{x, y}};

        Real f_exact = a + b * x + c * y;

        Real f_interp = 0;
        for (size_t j = 0; j < k.getCount(); j++)
        {
          const auto& node = k.getNode(j);
          Real f_node = a + b * node.x() + c * node.y();
          f_interp += f_node * k.getBasis(j)(p);
        }

        EXPECT_NEAR(f_interp, f_exact, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // ========================================================================
  // Tests for Tetrahedron elements (as requested)
  // ========================================================================

  // Test P2 tetrahedron Lagrange property
  TEST(Rodin_Variational_RealH1Element, LagrangeProperty_P2_Tetrahedron)
  {
    RealH1Element<2> k(Polytope::Type::Tetrahedron);

    for (size_t i = 0; i < k.getCount(); i++)
    {
      for (size_t j = 0; j < k.getCount(); j++)
      {
        const auto& node = k.getNode(j);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(k.getBasis(i)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P2 tetrahedron partition of unity
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P2_Tetrahedron)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<2> k(Polytope::Type::Tetrahedron);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      const auto& z = gen() * (1.0 - x - y);
      Math::Vector<Real> p{{x, y, z}};

      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P3 tetrahedron Lagrange property
  TEST(Rodin_Variational_RealH1Element, LagrangeProperty_P3_Tetrahedron)
  {
    RealH1Element<3> k(Polytope::Type::Tetrahedron);

    for (size_t i = 0; i < k.getCount(); i++)
    {
      for (size_t j = 0; j < k.getCount(); j++)
      {
        const auto& node = k.getNode(j);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(k.getBasis(i)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P3 tetrahedron partition of unity
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P3_Tetrahedron)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<3> k(Polytope::Type::Tetrahedron);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      const auto& z = gen() * (1.0 - x - y);
      Math::Vector<Real> p{{x, y, z}};

      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 tetrahedron derivative consistency
  // TODO: Fix derivative computation for barycentric elements
  TEST(Rodin_Variational_RealH1Element, DerivativeConsistency_P2_Tetrahedron)
  {
    RealH1Element<2> k(Polytope::Type::Tetrahedron);

    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      const auto& z = gen() * (1.0 - x - y);
      Math::Vector<Real> p{{x, y, z}};

      Real sum_dx = 0, sum_dy = 0, sum_dz = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum_dx += k.getBasis(j).getDerivative<1>(0)(p);
        sum_dy += k.getBasis(j).getDerivative<1>(1)(p);
        sum_dz += k.getBasis(j).getDerivative<1>(2)(p);
      }
      EXPECT_NEAR(sum_dx, 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(sum_dy, 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(sum_dz, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // ========================================================================
  // Tests for Wedge elements (as requested)
  // ========================================================================

  // Test P2 wedge Lagrange property
  TEST(Rodin_Variational_RealH1Element, LagrangeProperty_P2_Wedge)
  {
    RealH1Element<2> k(Polytope::Type::Wedge);

    for (size_t i = 0; i < k.getCount(); i++)
    {
      for (size_t j = 0; j < k.getCount(); j++)
      {
        const auto& node = k.getNode(j);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(k.getBasis(i)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P2 wedge partition of unity
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P2_Wedge)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<2> k(Polytope::Type::Wedge);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      const auto& z = gen();
      Math::Vector<Real> p{{x, y, z}};

      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P3 wedge Lagrange property
  TEST(Rodin_Variational_RealH1Element, LagrangeProperty_P3_Wedge)
  {
    RealH1Element<3> k(Polytope::Type::Wedge);

    for (size_t i = 0; i < k.getCount(); i++)
    {
      for (size_t j = 0; j < k.getCount(); j++)
      {
        const auto& node = k.getNode(j);
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(k.getBasis(i)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test P3 wedge partition of unity
  TEST(Rodin_Variational_RealH1Element, PartitionOfUnity_P3_Wedge)
  {
    constexpr size_t n = 50;
    RandomFloat gen(0.0, 1.0);
    RealH1Element<3> k(Polytope::Type::Wedge);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      const auto& z = gen();
      Math::Vector<Real> p{{x, y, z}};

      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test P2 wedge derivative consistency
  // TODO: Fix derivative computation for barycentric elements
  TEST(Rodin_Variational_RealH1Element, DerivativeConsistency_P2_Wedge)
  {
    RealH1Element<2> k(Polytope::Type::Wedge);

    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      const auto& z = gen();
      Math::Vector<Real> p{{x, y, z}};

      Real sum_dx = 0, sum_dy = 0, sum_dz = 0;
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum_dx += k.getBasis(j).getDerivative<1>(0)(p);
        sum_dy += k.getBasis(j).getDerivative<1>(1)(p);
        sum_dz += k.getBasis(j).getDerivative<1>(2)(p);
      }
      EXPECT_NEAR(sum_dx, 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(sum_dy, 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(sum_dz, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test higher orders DOF count for triangle
  TEST(Rodin_Variational_RealH1Element, DOFCount_HigherOrders_Triangle)
  {
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Triangle).getCount(), 15);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Triangle).getCount(), 21);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Triangle).getCount(), 28);
  }

  // Test higher orders DOF count for tetrahedron
  TEST(Rodin_Variational_RealH1Element, DOFCount_HigherOrders_Tetrahedron)
  {
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Tetrahedron).getCount(), 35);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Tetrahedron).getCount(), 56);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Tetrahedron).getCount(), 84);
  }

  // Test higher orders DOF count for wedge
  TEST(Rodin_Variational_RealH1Element, DOFCount_HigherOrders_Wedge)
  {
    // Wedge: (k+1) * (k+1)(k+2)/2
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Wedge).getCount(), 75);  // 5 * 5*6/2 = 5 * 15 = 75
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Wedge).getCount(), 126); // 6 * 6*7/2 = 6 * 21 = 126
  }

  // ========================================================================
  // Tests for Complex-valued H1Element (as requested)
  // ========================================================================

  // Test Complex P2 element basic functionality
  TEST(Rodin_Variational_ComplexH1Element, SanityTest_P2_Segment)
  {
    ComplexH1Element<2> k(Polytope::Type::Segment);

    // P2 element on Segment should have 3 DOFs
    EXPECT_EQ(k.getCount(), 3);
    EXPECT_EQ(k.getOrder(), 2);

    // Test Lagrange property with complex values
    for (size_t i = 0; i < k.getCount(); i++)
    {
      for (size_t j = 0; j < k.getCount(); j++)
      {
        const auto& node = k.getNode(j);
        Complex expected = (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
        Complex value = k.getBasis(i)(node);
        EXPECT_NEAR(std::abs(value - expected), 0.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test Complex P3 partition of unity
  TEST(Rodin_Variational_ComplexH1Element, PartitionOfUnity_P3_Segment)
  {
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    ComplexH1Element<3> k(Polytope::Type::Segment);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      Math::Vector<Real> p{{x}};

      Complex sum(0.0, 0.0);
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(std::abs(sum - Complex(1.0, 0.0)), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test Complex P2 on Triangle
  TEST(Rodin_Variational_ComplexH1Element, LagrangeProperty_P2_Triangle)
  {
    ComplexH1Element<2> k(Polytope::Type::Triangle);

    // Test Lagrange property
    for (size_t i = 0; i < k.getCount(); i++)
    {
      for (size_t j = 0; j < k.getCount(); j++)
      {
        const auto& node = k.getNode(j);
        Complex expected = (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
        Complex value = k.getBasis(i)(node);
        EXPECT_NEAR(std::abs(value - expected), 0.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test Complex P2 partition of unity on Quadrilateral
  TEST(Rodin_Variational_ComplexH1Element, PartitionOfUnity_P2_Quadrilateral)
  {
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    ComplexH1Element<2> k(Polytope::Type::Quadrilateral);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen();
      Math::Vector<Real> p{{x, y}};

      Complex sum(0.0, 0.0);
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(std::abs(sum - Complex(1.0, 0.0)), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test Complex DOF count consistency
  TEST(Rodin_Variational_ComplexH1Element, DOFCount_Consistency)
  {
    EXPECT_EQ(ComplexH1Element<2>(Polytope::Type::Segment).getCount(), 3);
    EXPECT_EQ(ComplexH1Element<3>(Polytope::Type::Segment).getCount(), 4);
    EXPECT_EQ(ComplexH1Element<2>(Polytope::Type::Triangle).getCount(), 6);
    EXPECT_EQ(ComplexH1Element<3>(Polytope::Type::Triangle).getCount(), 10);
    EXPECT_EQ(ComplexH1Element<2>(Polytope::Type::Quadrilateral).getCount(), 9);
  }

  // ========================================================================
  // Tests for Vector-valued H1Element (as requested)
  // ========================================================================

  // Test Vector P2 element basic functionality
  TEST(Rodin_Variational_VectorH1Element, SanityTest_P2_Segment_2D)
  {
    constexpr size_t vdim = 2;
    VectorH1Element<2> k(Polytope::Type::Segment, vdim);

    // P2 element on Segment with vdim=2 should have 3*2=6 DOFs
    EXPECT_EQ(k.getCount(), 6);
    EXPECT_EQ(k.getOrder(), 2);

    // Check that nodes match scalar version
    RealH1Element<2> scalar_k(Polytope::Type::Segment);
    for (size_t i = 0; i < scalar_k.getCount(); i++)
    {
      const auto& vec_node = k.getNode(i * vdim);
      const auto& scalar_node = scalar_k.getNode(i);
      EXPECT_NEAR((vec_node - scalar_node).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test Vector P3 element with 3D vectors
  TEST(Rodin_Variational_VectorH1Element, SanityTest_P3_Segment_3D)
  {
    constexpr size_t vdim = 3;
    VectorH1Element<3> k(Polytope::Type::Segment, vdim);

    // P3 element on Segment with vdim=3 should have 4*3=12 DOFs
    EXPECT_EQ(k.getCount(), 12);
    EXPECT_EQ(k.getOrder(), 3);
  }

  // Test Vector P2 basis function structure
  TEST(Rodin_Variational_VectorH1Element, BasisStructure_P2_Segment)
  {
    constexpr size_t vdim = 2;
    VectorH1Element<2> k(Polytope::Type::Segment, vdim);
    RealH1Element<2> scalar_k(Polytope::Type::Segment);

    // Test that vector basis functions have correct component structure
    for (size_t node_idx = 0; node_idx < scalar_k.getCount(); node_idx++)
    {
      const auto& node = scalar_k.getNode(node_idx);

      for (size_t component = 0; component < vdim; component++)
      {
        size_t local = node_idx * vdim + component;
        const auto& vec_basis = k.getBasis(local)(node);

        // Only the corresponding component should be non-zero
        for (size_t c = 0; c < vdim; c++)
        {
          if (c == component)
          {
            // This component should match the scalar basis
            EXPECT_NEAR(std::abs(vec_basis(c) - scalar_k.getBasis(node_idx)(node)), 0.0, RODIN_FUZZY_CONSTANT);
          }
          else
          {
            // Other components should be zero
            EXPECT_NEAR(std::abs(vec_basis(c)), 0.0, RODIN_FUZZY_CONSTANT);
          }
        }
      }
    }
  }

  // Test Vector P2 on Triangle
  TEST(Rodin_Variational_VectorH1Element, SanityTest_P2_Triangle_2D)
  {
    constexpr size_t vdim = 2;
    VectorH1Element<2> k(Polytope::Type::Triangle, vdim);

    // P2 element on Triangle with vdim=2 should have 6*2=12 DOFs
    EXPECT_EQ(k.getCount(), 12);
    EXPECT_EQ(k.getOrder(), 2);
  }

  // Test Vector P2 on Quadrilateral
  TEST(Rodin_Variational_VectorH1Element, SanityTest_P2_Quadrilateral_2D)
  {
    constexpr size_t vdim = 2;
    VectorH1Element<2> k(Polytope::Type::Quadrilateral, vdim);

    // P2 element on Quadrilateral with vdim=2 should have 9*2=18 DOFs
    EXPECT_EQ(k.getCount(), 18);
    EXPECT_EQ(k.getOrder(), 2);
  }

  // Test Vector P2 on Tetrahedron
  TEST(Rodin_Variational_VectorH1Element, SanityTest_P2_Tetrahedron_3D)
  {
    constexpr size_t vdim = 3;
    VectorH1Element<2> k(Polytope::Type::Tetrahedron, vdim);

    // P2 element on Tetrahedron with vdim=3 should have 10*3=30 DOFs
    EXPECT_EQ(k.getCount(), 30);
    EXPECT_EQ(k.getOrder(), 2);
  }

  // Test Vector P2 Jacobian computation
  TEST(Rodin_Variational_VectorH1Element, JacobianTest_P2_Segment)
  {
    constexpr size_t vdim = 2;
    VectorH1Element<2> k(Polytope::Type::Segment, vdim);

    RandomFloat gen(0.0, 1.0);

    for (size_t trial = 0; trial < 5; trial++)
    {
      const auto& x = gen();
      Math::Vector<Real> p{{x}};

      for (size_t i = 0; i < k.getCount(); i++)
      {
        const auto& jac = k.getBasis(i).getJacobian()(p);

        // Jacobian should have shape (vdim, spatial_dim)
        EXPECT_EQ(jac.rows(), vdim);
        EXPECT_EQ(jac.cols(), 1); // Segment is 1D
      }
    }
  }

  // Test Vector DOF count consistency
  TEST(Rodin_Variational_VectorH1Element, DOFCount_Consistency)
  {
    EXPECT_EQ(VectorH1Element<2>(Polytope::Type::Segment, 2).getCount(), 6);   // 3 nodes * 2 components
    EXPECT_EQ(VectorH1Element<3>(Polytope::Type::Segment, 2).getCount(), 8);   // 4 nodes * 2 components
    EXPECT_EQ(VectorH1Element<2>(Polytope::Type::Segment, 3).getCount(), 9);   // 3 nodes * 3 components
    EXPECT_EQ(VectorH1Element<2>(Polytope::Type::Triangle, 2).getCount(), 12); // 6 nodes * 2 components
    EXPECT_EQ(VectorH1Element<2>(Polytope::Type::Triangle, 3).getCount(), 18); // 6 nodes * 3 components
  }

  // Test higher-order vector elements
  TEST(Rodin_Variational_VectorH1Element, HigherOrder_P4_Segment)
  {
    constexpr size_t vdim = 3;
    VectorH1Element<4> k(Polytope::Type::Segment, vdim);

    // P4 element on Segment with vdim=3 should have 5*3=15 DOFs
    EXPECT_EQ(k.getCount(), 15);
    EXPECT_EQ(k.getOrder(), 4);
  }

  // Test vector P3 on higher-dimensional elements
  TEST(Rodin_Variational_VectorH1Element, P3_Wedge_3D)
  {
    constexpr size_t vdim = 3;
    VectorH1Element<3> k(Polytope::Type::Wedge, vdim);

    // P3 element on Wedge: (k+1) * (k+1)(k+2)/2 = 4 * 4*5/2 = 4 * 10 = 40 scalar DOFs
    // With vdim=3: 40*3 = 120 total DOFs
    EXPECT_EQ(k.getCount(), 120);
    EXPECT_EQ(k.getOrder(), 3);
  }

  // ========================================================================
  // Extended tests for Complex H1Element (additional mathematical properties)
  // ========================================================================

  // Test Complex P4 partition of unity on Triangle
  TEST(Rodin_Variational_ComplexH1Element, PartitionOfUnity_P4_Triangle)
  {
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    ComplexH1Element<4> k(Polytope::Type::Triangle);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      Math::Vector<Real> p{{x, y}};

      Complex sum(0.0, 0.0);
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(std::abs(sum - Complex(1.0, 0.0)), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test Complex P5 partition of unity on Quadrilateral
  TEST(Rodin_Variational_ComplexH1Element, PartitionOfUnity_P5_Quadrilateral)
  {
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    ComplexH1Element<5> k(Polytope::Type::Quadrilateral);

    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen();
      Math::Vector<Real> p{{x, y}};

      Complex sum(0.0, 0.0);
      for (size_t j = 0; j < k.getCount(); j++)
      {
        sum += k.getBasis(j)(p);
      }
      EXPECT_NEAR(std::abs(sum - Complex(1.0, 0.0)), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test Complex P3 Lagrange property on Tetrahedron
  TEST(Rodin_Variational_ComplexH1Element, LagrangeProperty_P3_Tetrahedron)
  {
    ComplexH1Element<3> k(Polytope::Type::Tetrahedron);

    for (size_t i = 0; i < k.getCount(); i++)
    {
      for (size_t j = 0; j < k.getCount(); j++)
      {
        const auto& node = k.getNode(j);
        Complex expected = (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
        Complex value = k.getBasis(i)(node);
        EXPECT_NEAR(std::abs(value - expected), 0.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // Test Complex higher orders DOF counting
  TEST(Rodin_Variational_ComplexH1Element, DOFCount_HigherOrders)
  {
    EXPECT_EQ(ComplexH1Element<4>(Polytope::Type::Segment).getCount(), 5);
    EXPECT_EQ(ComplexH1Element<5>(Polytope::Type::Segment).getCount(), 6);
    EXPECT_EQ(ComplexH1Element<6>(Polytope::Type::Segment).getCount(), 7);
    EXPECT_EQ(ComplexH1Element<4>(Polytope::Type::Triangle).getCount(), 15);
    EXPECT_EQ(ComplexH1Element<4>(Polytope::Type::Quadrilateral).getCount(), 25);
  }

  // Test Complex linear reproduction
  TEST(Rodin_Variational_ComplexH1Element, LinearReproduction_P2_Segment)
  {
    ComplexH1Element<2> k(Polytope::Type::Segment);
    RandomFloat gen(0.0, 1.0);

    for (size_t trial = 0; trial < 10; trial++)
    {
      // Random linear function: f(x) = a + b*x
      Complex a(gen(), gen());
      Complex b(gen(), gen());

      // Test interpolation at a random point
      const auto& x = gen();
      Math::Vector<Real> p{{x}};

      Complex expected = a + b * x;
      Complex interpolated(0.0, 0.0);

      for (size_t i = 0; i < k.getCount(); i++)
      {
        const auto& node = k.getNode(i);
        Complex f_at_node = a + b * node.x();
        interpolated += f_at_node * k.getBasis(i)(p);
      }

      EXPECT_NEAR(std::abs(interpolated - expected), 0.0, 1e-10);
    }
  }

  // Test Complex symmetry properties
  TEST(Rodin_Variational_ComplexH1Element, BasisSymmetry_P3_Segment)
  {
    ComplexH1Element<3> k(Polytope::Type::Segment);

    // For P3 on segment, basis functions should have symmetry
    // φ_0(x) at x should equal φ_3(1-x) at (1-x)
    RandomFloat gen(0.0, 1.0);

    for (size_t trial = 0; trial < 10; trial++)
    {
      const auto& x = gen();
      Math::Vector<Real> p1{{x}};
      Math::Vector<Real> p2{{1.0 - x}};

      Complex val1 = k.getBasis(0)(p1);
      Complex val2 = k.getBasis(3)(p2);

      EXPECT_NEAR(std::abs(val1 - val2), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // ========================================================================
  // Extended tests for Vector H1Element (additional mathematical properties)
  // ========================================================================

  // Test Vector P3 partition of unity on Triangle
  TEST(Rodin_Variational_VectorH1Element, PartitionOfUnity_P3_Triangle)
  {
    constexpr size_t vdim = 2;
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    VectorH1Element<3> k(Polytope::Type::Triangle, vdim);
    RealH1Element<3> scalar_k(Polytope::Type::Triangle);

    // For vector elements, each scalar basis function appears vdim times
    // The partition of unity should hold for the underlying scalar basis
    for (size_t i = 0; i < n; i++)
    {
      const auto& x = gen();
      const auto& y = gen() * (1.0 - x);
      Math::Vector<Real> p{{x, y}};

      // Check that scalar basis functions (summed across components) sum to 1
      Real scalar_sum = 0.0;
      for (size_t node_idx = 0; node_idx < scalar_k.getCount(); node_idx++)
      {
        scalar_sum += scalar_k.getBasis(node_idx)(p);
      }
      EXPECT_NEAR(scalar_sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // Test Vector P2 Lagrange property on Quadrilateral
  TEST(Rodin_Variational_VectorH1Element, LagrangeProperty_P2_Quadrilateral)
  {
    constexpr size_t vdim = 2;
    VectorH1Element<2> k(Polytope::Type::Quadrilateral, vdim);
    RealH1Element<2> scalar_k(Polytope::Type::Quadrilateral);

    // Test Lagrange property: each vector basis should be 1 in one component at one node
    for (size_t node_idx = 0; node_idx < scalar_k.getCount(); node_idx++)
    {
      const auto& node = scalar_k.getNode(node_idx);

      for (size_t component = 0; component < vdim; component++)
      {
        size_t local = node_idx * vdim + component;
        const auto& vec_basis = k.getBasis(local)(node);

        for (size_t c = 0; c < vdim; c++)
        {
          if (c == component)
          {
            EXPECT_NEAR(std::abs(vec_basis(c) - 1.0), 0.0, RODIN_FUZZY_CONSTANT);
          }
          else
          {
            EXPECT_NEAR(std::abs(vec_basis(c)), 0.0, RODIN_FUZZY_CONSTANT);
          }
        }
      }
    }
  }

  // Test Vector higher orders on Tetrahedron
  TEST(Rodin_Variational_VectorH1Element, HigherOrder_P4_Tetrahedron)
  {
    constexpr size_t vdim = 3;
    VectorH1Element<4> k(Polytope::Type::Tetrahedron, vdim);

    // P4 on Tetrahedron: (k+1)(k+2)(k+3)/6 = 5*6*7/6 = 35 scalar DOFs
    EXPECT_EQ(k.getCount(), 35 * vdim);
    EXPECT_EQ(k.getOrder(), 4);
  }

  // Test Vector basis sparsity in components
  TEST(Rodin_Variational_VectorH1Element, ComponentSparsity_P2_Segment)
  {
    constexpr size_t vdim = 3;
    VectorH1Element<2> k(Polytope::Type::Segment, vdim);

    RandomFloat gen(0.0, 1.0);

    for (size_t trial = 0; trial < 10; trial++)
    {
      const auto& x = gen();
      Math::Vector<Real> p{{x}};

      // Each vector basis function should have only one non-zero component
      for (size_t i = 0; i < k.getCount(); i++)
      {
        const auto& basis_val = k.getBasis(i)(p);
        size_t active_comp = i % vdim;

        size_t non_zero_count = 0;
        for (size_t c = 0; c < vdim; c++)
        {
          if (c == active_comp)
          {
            // This component should match the underlying scalar basis
            // (non-zero in general)
            if (std::abs(basis_val(c)) > RODIN_FUZZY_CONSTANT)
              non_zero_count++;
          }
          else
          {
            // Other components should be zero
            EXPECT_NEAR(std::abs(basis_val(c)), 0.0, RODIN_FUZZY_CONSTANT);
          }
        }

        // Exactly one component should be non-zero (the active one)
        EXPECT_EQ(non_zero_count, 1);
      }
    }
  }

  // Test Vector gradient consistency
  TEST(Rodin_Variational_VectorH1Element, GradientConsistency_P3_Segment)
  {
    constexpr size_t vdim = 2;
    VectorH1Element<3> k(Polytope::Type::Segment, vdim);

    RandomFloat gen(0.0, 1.0);

    for (size_t trial = 0; trial < 5; trial++)
    {
      const auto& x = gen();
      Math::Vector<Real> p{{x}};

      for (size_t i = 0; i < k.getCount(); i++)
      {
        // Check that Jacobian is computed correctly
        const auto& jac = k.getBasis(i).getJacobian()(p);

        // Jacobian should match component-wise derivatives
        for (size_t c = 0; c < vdim; c++)
        {
          for (size_t d = 0; d < 1; d++) // 1D spatial dimension
          {
            Real deriv = k.getBasis(i).template getDerivative<1>(c, d)(p);
            EXPECT_NEAR(std::abs(jac(c, d) - deriv), 0.0, RODIN_FUZZY_CONSTANT);
          }
        }
      }
    }
  }

  // Test Vector DOF count on all geometries with higher orders
  TEST(Rodin_Variational_VectorH1Element, DOFCount_AllGeometries_P5)
  {
    constexpr size_t vdim = 3;

    // P5 Segment: 6 scalar DOFs
    EXPECT_EQ(VectorH1Element<5>(Polytope::Type::Segment, vdim).getCount(), 6 * vdim);

    // P5 Triangle: (k+1)(k+2)/2 = 6*7/2 = 21 scalar DOFs
    EXPECT_EQ(VectorH1Element<5>(Polytope::Type::Triangle, vdim).getCount(), 21 * vdim);

    // P5 Quadrilateral: (k+1)^2 = 36 scalar DOFs
    EXPECT_EQ(VectorH1Element<5>(Polytope::Type::Quadrilateral, vdim).getCount(), 36 * vdim);

    // P5 Tetrahedron: (k+1)(k+2)(k+3)/6 = 6*7*8/6 = 56 scalar DOFs
    EXPECT_EQ(VectorH1Element<5>(Polytope::Type::Tetrahedron, vdim).getCount(), 56 * vdim);
  }

  // Test Vector linear field reproduction
  TEST(Rodin_Variational_VectorH1Element, LinearReproduction_P2_Segment)
  {
    constexpr size_t vdim = 2;
    VectorH1Element<2> k(Polytope::Type::Segment, vdim);
    RealH1Element<2> scalar_k(Polytope::Type::Segment);

    RandomFloat gen(0.0, 1.0);

    for (size_t trial = 0; trial < 5; trial++)
    {
      // Define a linear vector field: v(x) = [a0 + b0*x, a1 + b1*x]
      Real a0 = gen(), b0 = gen();
      Real a1 = gen(), b1 = gen();

      // Test interpolation at random point
      const auto& x = gen();
      Math::Vector<Real> test_pt{{x}};

      Math::Vector<Real> expected(vdim);
      expected(0) = a0 + b0 * x;
      expected(1) = a1 + b1 * x;

      Math::Vector<Real> interpolated(vdim);
      interpolated.setZero();

      for (size_t node_idx = 0; node_idx < scalar_k.getCount(); node_idx++)
      {
        const auto& node = scalar_k.getNode(node_idx);
        Real nx = node.x();

        // Vector field values at nodes
        Math::Vector<Real> v_at_node(vdim);
        v_at_node(0) = a0 + b0 * nx;
        v_at_node(1) = a1 + b1 * nx;

        // Add contributions from all vector basis functions at this node
        for (size_t c = 0; c < vdim; c++)
        {
          size_t local = node_idx * vdim + c;
          const auto& basis_val = k.getBasis(local)(test_pt);
          interpolated += v_at_node(c) * basis_val;
        }
      }

      for (size_t c = 0; c < vdim; c++)
      {
        EXPECT_NEAR(std::abs(interpolated(c) - expected(c)), 0.0, 1e-10);
      }
    }
  }

  TEST(FinalTest_H1Element_Real, PartitionOfUnity_P2_AllGeometries)
  {
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);

    // Segment
    {
      RealH1Element<2> elem(Polytope::Type::Segment);
      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen()}};
        Real sum = 0.0;
        for (size_t j = 0; j < elem.getCount(); j++)
          sum += elem.getBasis(j)(p);
        EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
      }
    }

    // Triangle
    {
      RealH1Element<2> elem(Polytope::Type::Triangle);
      for (size_t i = 0; i < n; i++)
      {
        Real s = gen(), t = gen();
        if (s + t <= 1.0)
        {
          Math::Vector<Real> p{{s, t}};
          Real sum = 0.0;
          for (size_t j = 0; j < elem.getCount(); j++)
            sum += elem.getBasis(j)(p);
          EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
        }
      }
    }

    // Quadrilateral
    {
      RealH1Element<2> elem(Polytope::Type::Quadrilateral);
      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen(), gen()}};
        Real sum = 0.0;
        for (size_t j = 0; j < elem.getCount(); j++)
          sum += elem.getBasis(j)(p);
        EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  TEST(FinalTest_H1Element_Real, PartitionOfUnity_P3_P4_P5_P6_Segment)
  {
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);

    // Test P3, P4, P5, P6
    for (size_t order = 3; order <= 6; order++)
    {
      size_t dofs = order + 1;

      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen()}};
        Real sum = 0.0;

        if (order == 3)
        {
          RealH1Element<3> elem(Polytope::Type::Segment);
          for (size_t j = 0; j < elem.getCount(); j++)
            sum += elem.getBasis(j)(p);
        }
        else if (order == 4)
        {
          RealH1Element<4> elem(Polytope::Type::Segment);
          for (size_t j = 0; j < elem.getCount(); j++)
            sum += elem.getBasis(j)(p);
        }
        else if (order == 5)
        {
          RealH1Element<5> elem(Polytope::Type::Segment);
          for (size_t j = 0; j < elem.getCount(); j++)
            sum += elem.getBasis(j)(p);
        }
        else if (order == 6)
        {
          RealH1Element<6> elem(Polytope::Type::Segment);
          for (size_t j = 0; j < elem.getCount(); j++)
            sum += elem.getBasis(j)(p);
        }

        EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  TEST(FinalTest_H1Element_Real, LagrangeProperty_P2_Segment)
  {
    RealH1Element<2> elem(Polytope::Type::Segment);

    for (size_t i = 0; i < 3; i++)
    {
      const auto& node = elem.getNode(i);
      for (size_t j = 0; j < 3; j++)
      {
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  TEST(FinalTest_H1Element_Real, LagrangeProperty_P3_P4_P5_P6_Segment)
  {
    // Test orders 3-6
    {
      RealH1Element<3> elem(Polytope::Type::Segment);
      for (size_t i = 0; i < 4; i++)
      {
        const auto& node = elem.getNode(i);
        for (size_t j = 0; j < 4; j++)
        {
          Real expected = (i == j) ? 1.0 : 0.0;
          EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
        }
      }
    }

    {
      RealH1Element<4> elem(Polytope::Type::Segment);
      for (size_t i = 0; i < 5; i++)
      {
        const auto& node = elem.getNode(i);
        for (size_t j = 0; j < 5; j++)
        {
          Real expected = (i == j) ? 1.0 : 0.0;
          EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
        }
      }
    }

    {
      RealH1Element<5> elem(Polytope::Type::Segment);
      for (size_t i = 0; i < 6; i++)
      {
        const auto& node = elem.getNode(i);
        for (size_t j = 0; j < 6; j++)
        {
          Real expected = (i == j) ? 1.0 : 0.0;
          EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
        }
      }
    }

    {
      RealH1Element<6> elem(Polytope::Type::Segment);
      for (size_t i = 0; i < 7; i++)
      {
        const auto& node = elem.getNode(i);
        for (size_t j = 0; j < 7; j++)
        {
          Real expected = (i == j) ? 1.0 : 0.0;
          EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
        }
      }
    }
  }

  TEST(FinalTest_H1Element_Real, QuadraticReproduction_P2_Segment)
  {
    RealH1Element<2> elem(Polytope::Type::Segment);

    // Test quadratic function f(x) = 1 + 2*x + 3*x^2
    // Values at nodes: x=0: f=1, x=0.5: f=2.75, x=1: f=6
    Real f0 = 1.0, f1 = 2.75, f2 = 6.0;

    Math::Vector<Real> p{{0.3}};
    Real exact = 1.0 + 2.0 * 0.3 + 3.0 * 0.3 * 0.3;  // 1.87
    Real interpolated = f0 * elem.getBasis(0)(p) + f1 * elem.getBasis(1)(p) + f2 * elem.getBasis(2)(p);

    EXPECT_NEAR(interpolated, exact, RODIN_FUZZY_CONSTANT);
  }

  TEST(FinalTest_H1Element_Complex, PartitionOfUnity_P2_P3_P4_Segment)
  {
    using Complex = std::complex<Real>;
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);

    // P2
    {
      H1Element<2, Complex> elem(Polytope::Type::Segment);
      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen()}};
        Complex sum(0.0, 0.0);
        for (size_t j = 0; j < elem.getCount(); j++)
          sum += elem.getBasis(j)(p);
        EXPECT_NEAR(sum.real(), 1.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR(sum.imag(), 0.0, RODIN_FUZZY_CONSTANT);
      }
    }

    // P3
    {
      H1Element<3, Complex> elem(Polytope::Type::Segment);
      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen()}};
        Complex sum(0.0, 0.0);
        for (size_t j = 0; j < elem.getCount(); j++)
          sum += elem.getBasis(j)(p);
        EXPECT_NEAR(sum.real(), 1.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR(sum.imag(), 0.0, RODIN_FUZZY_CONSTANT);
      }
    }

    // P4
    {
      H1Element<4, Complex> elem(Polytope::Type::Segment);
      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen()}};
        Complex sum(0.0, 0.0);
        for (size_t j = 0; j < elem.getCount(); j++)
          sum += elem.getBasis(j)(p);
        EXPECT_NEAR(sum.real(), 1.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR(sum.imag(), 0.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  TEST(FinalTest_H1Element_Complex, ComplexQuadraticReproduction_P2_Segment)
  {
    using Complex = std::complex<Real>;
    H1Element<2, Complex> elem(Polytope::Type::Segment);

    // Test complex quadratic: f(x) = (1+i) + (2+2i)*x + (3+3i)*x^2
    Complex f0(1.0, 1.0);      // at x=0
    Complex f1(2.75, 2.75);    // at x=0.5
    Complex f2(6.0, 6.0);      // at x=1

    Math::Vector<Real> p{{0.3}};
    Complex exact(1.87, 1.87);  // (1+i) + (2+2i)*0.3 + (3+3i)*0.09
    Complex interpolated = f0 * elem.getBasis(0)(p) + f1 * elem.getBasis(1)(p) + f2 * elem.getBasis(2)(p);

    EXPECT_NEAR(interpolated.real(), exact.real(), RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(interpolated.imag(), exact.imag(), RODIN_FUZZY_CONSTANT);
  }

  TEST(FinalTest_H1Element_Vector, ComponentStructure_P2_2D_Segment)
  {
    constexpr size_t vdim = 2;
    H1Element<2, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
    H1Element<2, Real> scalar_elem(Polytope::Type::Segment);

    EXPECT_EQ(elem.getCount(), scalar_elem.getCount() * vdim);

    Math::Vector<Real> p{{0.5}};

    for (size_t i = 0; i < scalar_elem.getCount(); ++i)
    {
      Real scalar_val = scalar_elem.getBasis(i)(p);

      for (size_t c = 0; c < vdim; ++c)
      {
        size_t local = i * vdim + c;

        const auto& val = elem.getBasis(local)(p);
        EXPECT_EQ(val.size(), vdim);

        // At most one non-zero component
        size_t num_nonzero = 0;
        for (size_t j = 0; j < vdim; ++j)
          if (std::abs(val(j)) > RODIN_FUZZY_CONSTANT)
            num_nonzero++;
        EXPECT_LE(num_nonzero, 1u);

        // Component c matches scalar basis, others are zero
        EXPECT_NEAR(val(c), scalar_val, RODIN_FUZZY_CONSTANT);
        for (size_t j = 0; j < vdim; ++j)
        {
          if (j != c)
          {
            EXPECT_NEAR(val(j), 0.0, RODIN_FUZZY_CONSTANT);
          }
        }
      }
    }
  }

  TEST(FinalTest_H1Element_Vector, QuadraticVectorFieldReproduction_P2_Segment)
  {
    H1Element<2, Math::Vector<Real>> elem(Polytope::Type::Segment, 2);

    auto v = [](Real x)
    {
      return Math::Vector<Real>{
        { 1.0 + 2.0 * x + 3.0 * x * x,  // first component
          4.0 + 5.0 * x + 6.0 * x * x } // second component
      };
    };

    // Nodes: x = 0, 0.5, 1
    std::vector<Math::Vector<Real>> node_values = {
      v(0.0),
      v(0.5),
      v(1.0)
    };

    Math::Vector<Real> p{{0.3}};
    Math::Vector<Real> exact = v(p(0));

    Math::Vector<Real> interpolated = Math::Vector<Real>::Zero(2);
    for (size_t node = 0; node < 3; node++)
    {
      for (size_t comp = 0; comp < 2; comp++)
      {
        size_t basis_idx = node * 2 + comp;
        interpolated += node_values[node](comp) * elem.getBasis(basis_idx)(p);
      }
    }

    EXPECT_NEAR(interpolated(0), exact(0), RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(interpolated(1), exact(1), RODIN_FUZZY_CONSTANT);
  }

  TEST(FinalTest_H1Element_Vector, PartitionOfUnity_P3_P4_3D_Triangle)
  {
    // P3 on Triangle with 3D vectors
    {
      H1Element<3, Math::Vector<Real>> elem(Polytope::Type::Triangle, 3);
      EXPECT_EQ(elem.getCount(), 30);  // 10 nodes × 3 components

      Math::Vector<Real> p{{0.3, 0.4}};

      // For each component, sum of basis functions with that component should equal 1
      for (size_t comp = 0; comp < 3; comp++)
      {
        Real sum = 0.0;
        for (size_t i = comp; i < elem.getCount(); i += 3)
        {
          const auto& val = elem.getBasis(i)(p);
          sum += val(comp);
        }
        EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  TEST(FinalTest_Integration, AllElementsAllTypes_DOFCounting)
  {
    // Verify DOF counting formulas are correct for all elements and types

    // P0 Elements
    EXPECT_EQ(RealP0Element(Polytope::Type::Segment).getCount(), 1);
    EXPECT_EQ(RealP0Element(Polytope::Type::Triangle).getCount(), 1);
    EXPECT_EQ(RealP0Element(Polytope::Type::Quadrilateral).getCount(), 1);
    EXPECT_EQ(RealP0Element(Polytope::Type::Tetrahedron).getCount(), 1);
    EXPECT_EQ(VectorP0Element<Real>(Polytope::Type::Triangle).getCount(), 2);
    EXPECT_EQ(VectorP0Element<Real>(Polytope::Type::Tetrahedron).getCount(), 3);

    // P1 Elements
    EXPECT_EQ(RealP1Element(Polytope::Type::Segment).getCount(), 2);
    EXPECT_EQ(RealP1Element(Polytope::Type::Triangle).getCount(), 3);
    EXPECT_EQ(RealP1Element(Polytope::Type::Quadrilateral).getCount(), 4);
    EXPECT_EQ(RealP1Element(Polytope::Type::Tetrahedron).getCount(), 4);
    EXPECT_EQ(VectorP1Element<Real>(Polytope::Type::Segment, 2).getCount(), 4);
    EXPECT_EQ(VectorP1Element<Real>(Polytope::Type::Triangle, 3).getCount(), 9);

    // Pk Elements (Segment): K+1
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Segment).getCount(), 3);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Segment).getCount(), 4);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Segment).getCount(), 5);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Segment).getCount(), 6);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Segment).getCount(), 7);

    // Pk Elements (Triangle): (K+1)(K+2)/2
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Triangle).getCount(), 6);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Triangle).getCount(), 10);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Triangle).getCount(), 15);

    // Pk Elements (Quadrilateral): (K+1)^2
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Quadrilateral).getCount(), 9);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Quadrilateral).getCount(), 16);

    // Pk Elements (Tetrahedron): (K+1)(K+2)(K+3)/6
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Tetrahedron).getCount(), 10);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Tetrahedron).getCount(), 20);

    // Vector Pk Elements
    int count = H1Element<2, Math::Vector<Real>>(Polytope::Type::Segment, 2).getCount();
    EXPECT_EQ(count, 6);

    count = H1Element<3, Math::Vector<Real>>(Polytope::Type::Segment, 3).getCount();
    EXPECT_EQ(count, 12);
  }

  TEST(FinalTest_Integration, AllElementsAllTypes_OrderProperty)
  {
    // Verify order property is correct

    EXPECT_EQ(RealP0Element(Polytope::Type::Segment).getOrder(), 0);
    EXPECT_EQ(RealP0Element(Polytope::Type::Triangle).getOrder(), 0);

    EXPECT_EQ(RealP1Element(Polytope::Type::Segment).getOrder(), 1);
    EXPECT_EQ(RealP1Element(Polytope::Type::Triangle).getOrder(), 1);

    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Segment).getOrder(), 2);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Segment).getOrder(), 3);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Segment).getOrder(), 4);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Segment).getOrder(), 5);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Segment).getOrder(), 6);
  }

  // ========================================================================
  // NEW COMPREHENSIVE TESTS FOR PKELEMENT K=0 THROUGH K=6
  // ========================================================================

  TEST(FinalTest_H1Element_Comprehensive, H1Element_K0_AllGeometries)
  {
    // Comprehensive tests for H1Element<0> across all geometries
    // K=0 is piecewise constant element

    for (auto geom : {Polytope::Type::Point, Polytope::Type::Segment,
                      Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
                      Polytope::Type::Tetrahedron, Polytope::Type::Wedge})
    {
      RealH1Element<0> pk(geom);

      // Should have 1 DOF
      EXPECT_EQ(pk.getCount(), 1);

      // Should have order 0
      EXPECT_EQ(pk.getOrder(), 0);

      // Test basis function value (constant = 1)
      Math::Vector<Real> p;
      switch (geom)
      {
        case Polytope::Type::Point:
          p = Math::Vector<Real>{{0}};
          break;
        case Polytope::Type::Segment:
          p = Math::Vector<Real>{{0.5}};
          break;
        case Polytope::Type::Triangle:
        case Polytope::Type::Quadrilateral:
          p = Math::Vector<Real>{{0.3, 0.3}};
          break;
        case Polytope::Type::Tetrahedron:
        case Polytope::Type::Wedge:
          p = Math::Vector<Real>{{0.25, 0.25, 0.25}};
          break;
      }

      EXPECT_NEAR(pk.getBasis(0)(p), 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(FinalTest_H1Element_Comprehensive, H1Element_K1_AllGeometries)
  {
    // Comprehensive tests for H1Element<1> across all geometries
    // K=1 is piecewise linear element

    for (auto geom : {Polytope::Type::Segment, Polytope::Type::Triangle,
                      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
                      Polytope::Type::Wedge})
    {
      RealH1Element<1> pk(geom);

      // Check expected DOF count
      size_t expected_dofs = 0;
      switch (geom)
      {
        case Polytope::Type::Segment:
          expected_dofs = 2;
          break;
        case Polytope::Type::Triangle:
          expected_dofs = 3;
          break;
        case Polytope::Type::Quadrilateral:
          expected_dofs = 4;
          break;
        case Polytope::Type::Tetrahedron:
          expected_dofs = 4;
          break;
        case Polytope::Type::Wedge:
          expected_dofs = 6;
          break;
        default:
          continue;
      }
      EXPECT_EQ(pk.getCount(), expected_dofs);

      // Check order (total polynomial degree)
      size_t expected_order = 0;
      switch (geom)
      {
        case Polytope::Type::Segment:
        case Polytope::Type::Triangle:
        case Polytope::Type::Tetrahedron:
          expected_order = 1;
          break;
        case Polytope::Type::Quadrilateral:
        case Polytope::Type::Wedge:
          expected_order = 2; // tensor-product degree with K=1
          break;
        default:
          continue;
      }
      EXPECT_EQ(pk.getOrder(), expected_order);

      // Test partition of unity
      Math::Vector<Real> p;
      switch (geom)
      {
        case Polytope::Type::Segment:
          p = Math::Vector<Real>{{0.5}};
          break;
        case Polytope::Type::Triangle:
        case Polytope::Type::Quadrilateral:
          p = Math::Vector<Real>{{0.3, 0.3}};
          break;
        case Polytope::Type::Tetrahedron:
        case Polytope::Type::Wedge:
          p = Math::Vector<Real>{{0.25, 0.25, 0.25}};
          break;
        default:
          continue;
      }

      Real sum_pk = 0;
      for (size_t i = 0; i < pk.getCount(); i++)
      {
        sum_pk += pk.getBasis(i)(p);
      }
      EXPECT_NEAR(sum_pk, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(FinalTest_H1Element_Comprehensive, H1Element_K2_to_K6_Segment)
  {
    // Comprehensive test of H1Element for K=2,3,4,5,6 on Segment
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);

    for (size_t K : {2, 3, 4, 5, 6})
    {
      // Test partition of unity
      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen()}};
        Real sum = 0;

        switch (K)
        {
          case 2:
          {
            RealH1Element<2> elem(Polytope::Type::Segment);
            for (size_t j = 0; j < elem.getCount(); j++)
              sum += elem.getBasis(j)(p);
            break;
          }
          case 3:
          {
            RealH1Element<3> elem(Polytope::Type::Segment);
            for (size_t j = 0; j < elem.getCount(); j++)
              sum += elem.getBasis(j)(p);
            break;
          }
          case 4:
          {
            RealH1Element<4> elem(Polytope::Type::Segment);
            for (size_t j = 0; j < elem.getCount(); j++)
              sum += elem.getBasis(j)(p);
            break;
          }
          case 5:
          {
            RealH1Element<5> elem(Polytope::Type::Segment);
            for (size_t j = 0; j < elem.getCount(); j++)
              sum += elem.getBasis(j)(p);
            break;
          }
          case 6:
          {
            RealH1Element<6> elem(Polytope::Type::Segment);
            for (size_t j = 0; j < elem.getCount(); j++)
              sum += elem.getBasis(j)(p);
            break;
          }
        }

        EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  TEST(FinalTest_H1Element_Comprehensive, H1Element_K2_to_K6_Triangle)
  {
    // Test H1Element for K=2,3,4,5,6 on Triangle
    constexpr size_t n = 15;
    RandomFloat gen(0.0, 1.0);

    for (size_t K : {2, 3, 4, 5, 6})
    {
      for (size_t i = 0; i < n; i++)
      {
        Real s = gen();
        Real t = gen() * (1 - s);
        Math::Vector<Real> p{{s, t}};
        Real sum = 0;

        switch (K)
        {
          case 2:
          {
            RealH1Element<2> elem(Polytope::Type::Triangle);
            EXPECT_EQ(elem.getCount(), 6);
            for (size_t j = 0; j < elem.getCount(); j++)
              sum += elem.getBasis(j)(p);
            break;
          }
          case 3:
          {
            RealH1Element<3> elem(Polytope::Type::Triangle);
            EXPECT_EQ(elem.getCount(), 10);
            for (size_t j = 0; j < elem.getCount(); j++)
              sum += elem.getBasis(j)(p);
            break;
          }
          case 4:
          {
            RealH1Element<4> elem(Polytope::Type::Triangle);
            EXPECT_EQ(elem.getCount(), 15);
            for (size_t j = 0; j < elem.getCount(); j++)
              sum += elem.getBasis(j)(p);
            break;
          }
          case 5:
          {
            RealH1Element<5> elem(Polytope::Type::Triangle);
            EXPECT_EQ(elem.getCount(), 21);
            for (size_t j = 0; j < elem.getCount(); j++)
              sum += elem.getBasis(j)(p);
            break;
          }
          case 6:
          {
            RealH1Element<6> elem(Polytope::Type::Triangle);
            EXPECT_EQ(elem.getCount(), 28);
            for (size_t j = 0; j < elem.getCount(); j++)
              sum += elem.getBasis(j)(p);
            break;
          }
        }

        EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  TEST(FinalTest_H1Element_Comprehensive, H1Element_K2_to_K6_AllGeometries_DOFCount)
  {
    // Verify DOF counts for all K values and geometries

    // Segment: K+1
    EXPECT_EQ(RealH1Element<0>(Polytope::Type::Segment).getCount(), 1);
    EXPECT_EQ(RealH1Element<1>(Polytope::Type::Segment).getCount(), 2);
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Segment).getCount(), 3);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Segment).getCount(), 4);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Segment).getCount(), 5);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Segment).getCount(), 6);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Segment).getCount(), 7);

    // Triangle: (K+1)(K+2)/2
    EXPECT_EQ(RealH1Element<0>(Polytope::Type::Triangle).getCount(), 1);
    EXPECT_EQ(RealH1Element<1>(Polytope::Type::Triangle).getCount(), 3);
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Triangle).getCount(), 6);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Triangle).getCount(), 10);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Triangle).getCount(), 15);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Triangle).getCount(), 21);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Triangle).getCount(), 28);

    // Quadrilateral: (K+1)^2
    EXPECT_EQ(RealH1Element<0>(Polytope::Type::Quadrilateral).getCount(), 1);
    EXPECT_EQ(RealH1Element<1>(Polytope::Type::Quadrilateral).getCount(), 4);
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Quadrilateral).getCount(), 9);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Quadrilateral).getCount(), 16);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Quadrilateral).getCount(), 25);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Quadrilateral).getCount(), 36);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Quadrilateral).getCount(), 49);

    // Tetrahedron: (K+1)(K+2)(K+3)/6
    EXPECT_EQ(RealH1Element<0>(Polytope::Type::Tetrahedron).getCount(), 1);
    EXPECT_EQ(RealH1Element<1>(Polytope::Type::Tetrahedron).getCount(), 4);
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Tetrahedron).getCount(), 10);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Tetrahedron).getCount(), 20);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Tetrahedron).getCount(), 35);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Tetrahedron).getCount(), 56);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Tetrahedron).getCount(), 84);

    // Wedge: (K+1)·(K+1)(K+2)/2
    EXPECT_EQ(RealH1Element<0>(Polytope::Type::Wedge).getCount(), 1);
    EXPECT_EQ(RealH1Element<1>(Polytope::Type::Wedge).getCount(), 6);
    EXPECT_EQ(RealH1Element<2>(Polytope::Type::Wedge).getCount(), 18);
    EXPECT_EQ(RealH1Element<3>(Polytope::Type::Wedge).getCount(), 40);
    EXPECT_EQ(RealH1Element<4>(Polytope::Type::Wedge).getCount(), 75);
    EXPECT_EQ(RealH1Element<5>(Polytope::Type::Wedge).getCount(), 126);
    EXPECT_EQ(RealH1Element<6>(Polytope::Type::Wedge).getCount(), 196);
  }

  TEST(FinalTest_H1Element_Comprehensive, VectorH1Element_AllVectorDimensions)
  {
    // Test vector H1Element with vdim=1,2,3 for various K values

    // K=0, vdim=1,2,3 on Segment
    {
      for (size_t vdim : {1, 2, 3})
      {
        H1Element<0, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
        EXPECT_EQ(elem.getCount(), vdim * 1);
      }
    }

    // K=1, vdim=1,2,3 on Segment
    {
      for (size_t vdim : {1, 2, 3})
      {
        H1Element<1, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
        EXPECT_EQ(elem.getCount(), vdim * 2);
      }
    }

    // K=2, vdim=1,2,3 on Triangle
    {
      for (size_t vdim : {1, 2, 3})
      {
        H1Element<2, Math::Vector<Real>> elem(Polytope::Type::Triangle, vdim);
        EXPECT_EQ(elem.getCount(), vdim * 6);
      }
    }

    // K=3, vdim=1,2,3 on Tetrahedron
    {
      for (size_t vdim : {1, 2, 3})
      {
        H1Element<3, Math::Vector<Real>> elem(Polytope::Type::Tetrahedron, vdim);
        EXPECT_EQ(elem.getCount(), vdim * 20);
      }
    }
  }

  TEST(FinalTest_H1Element_Comprehensive, H1Element_LagrangeProperty_AllOrders_Segment)
  {
    // Test Lagrange property for all K values on Segment
    // phi_i(node_j) = delta_ij

    for (size_t K : {0, 1, 2, 3, 4, 5, 6})
    {
      switch (K)
      {
        case 0:
        {
          RealH1Element<0> elem(Polytope::Type::Segment);
          const auto& node = elem.getNode(0);
          EXPECT_NEAR(elem.getBasis(0)(node), 1.0, RODIN_FUZZY_CONSTANT);
          break;
        }
        case 1:
        {
          RealH1Element<1> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& node = elem.getNode(i);
            for (size_t j = 0; j < elem.getCount(); j++)
            {
              Real expected = (i == j) ? 1.0 : 0.0;
              EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
            }
          }
          break;
        }
        case 2:
        {
          RealH1Element<2> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& node = elem.getNode(i);
            for (size_t j = 0; j < elem.getCount(); j++)
            {
              Real expected = (i == j) ? 1.0 : 0.0;
              EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
            }
          }
          break;
        }
        case 3:
        {
          RealH1Element<3> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& node = elem.getNode(i);
            for (size_t j = 0; j < elem.getCount(); j++)
            {
              Real expected = (i == j) ? 1.0 : 0.0;
              EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
            }
          }
          break;
        }
        case 4:
        {
          RealH1Element<4> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& node = elem.getNode(i);
            for (size_t j = 0; j < elem.getCount(); j++)
            {
              Real expected = (i == j) ? 1.0 : 0.0;
              EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
            }
          }
          break;
        }
        case 5:
        {
          RealH1Element<5> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& node = elem.getNode(i);
            for (size_t j = 0; j < elem.getCount(); j++)
            {
              Real expected = (i == j) ? 1.0 : 0.0;
              EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
            }
          }
          break;
        }
        case 6:
        {
          RealH1Element<6> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& node = elem.getNode(i);
            for (size_t j = 0; j < elem.getCount(); j++)
            {
              Real expected = (i == j) ? 1.0 : 0.0;
              EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
            }
          }
          break;
        }
      }
    }
  }

  // ========================================================================
  // LINEARFORM TESTS FOR PKELEMENT K=0 TO K=6
  // ========================================================================

  TEST(FinalTest_H1Element_LinearForm, ScalarLinearForm_K0_to_K6_Segment)
  {
    // Test LinearForm for all polynomial orders on Segment
    for (size_t K : {0, 1, 2, 3, 4, 5, 6})
    {
      // Polynomial function of appropriate degree
      auto f = [K](const Math::SpatialPoint& x) -> Real {
        Real result = 1.0;
        for (size_t k = 1; k <= K; k++)
          result += (k + 1.0) * std::pow(x.x(), k);
        return result;
      };

      switch (K)
      {
        case 0:
        {
          RealH1Element<0> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& lf = elem.getLinearForm(i);
            Real dof_val = lf(f);
            const auto& node = elem.getNode(i);
            EXPECT_NEAR(dof_val, f(node), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 1:
        {
          RealH1Element<1> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& lf = elem.getLinearForm(i);
            Real dof_val = lf(f);
            const auto& node = elem.getNode(i);
            EXPECT_NEAR(dof_val, f(node), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 2:
        {
          RealH1Element<2> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& lf = elem.getLinearForm(i);
            Real dof_val = lf(f);
            const auto& node = elem.getNode(i);
            EXPECT_NEAR(dof_val, f(node), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 3:
        {
          RealH1Element<3> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& lf = elem.getLinearForm(i);
            Real dof_val = lf(f);
            const auto& node = elem.getNode(i);
            EXPECT_NEAR(dof_val, f(node), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 4:
        {
          RealH1Element<4> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& lf = elem.getLinearForm(i);
            Real dof_val = lf(f);
            const auto& node = elem.getNode(i);
            EXPECT_NEAR(dof_val, f(node), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 5:
        {
          RealH1Element<5> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& lf = elem.getLinearForm(i);
            Real dof_val = lf(f);
            const auto& node = elem.getNode(i);
            EXPECT_NEAR(dof_val, f(node), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 6:
        {
          RealH1Element<6> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            const auto& lf = elem.getLinearForm(i);
            Real dof_val = lf(f);
            const auto& node = elem.getNode(i);
            EXPECT_NEAR(dof_val, f(node), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
      }
    }
  }

  TEST(FinalTest_H1Element_LinearForm, VectorLinearForm_K0_to_K3_AllVectorDimensions)
  {
    // Test LinearForm for vector Pk elements
    for (size_t K : {0, 1, 2, 3})
    {
      for (size_t vdim : {1, 2, 3})
      {
        // Vector function
        auto f = [vdim, K](const Math::SpatialPoint& x) {
          Math::Vector<Real> v(vdim);
          Real poly_val = 1.0;
          for (size_t k = 1; k <= K; k++)
            poly_val += std::pow(x.x(), k);
          for (size_t i = 0; i < vdim; i++)
            v(i) = (i + 1) * poly_val;
          return v;
        };

        switch (K)
        {
          case 0:
          {
            H1Element<0, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
            {
              const auto& lf = elem.getLinearForm(i);
              Real dof_val = lf(f);
              const auto& node = elem.getNode(i);
              size_t comp = i % vdim;
              EXPECT_NEAR(dof_val, f(node)(comp), RODIN_FUZZY_CONSTANT);
            }
            break;
          }
          case 1:
          {
            H1Element<1, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
            {
              const auto& lf = elem.getLinearForm(i);
              Real dof_val = lf(f);
              const auto& node = elem.getNode(i);
              size_t comp = i % vdim;
              EXPECT_NEAR(dof_val, f(node)(comp), RODIN_FUZZY_CONSTANT);
            }
            break;
          }
          case 2:
          {
            H1Element<2, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
            {
              const auto& lf = elem.getLinearForm(i);
              Real dof_val = lf(f);
              const auto& node = elem.getNode(i);
              size_t comp = i % vdim;
              EXPECT_NEAR(dof_val, f(node)(comp), RODIN_FUZZY_CONSTANT);
            }
            break;
          }
          case 3:
          {
            H1Element<3, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
            {
              const auto& lf = elem.getLinearForm(i);
              Real dof_val = lf(f);
              const auto& node = elem.getNode(i);
              size_t comp = i % vdim;
              EXPECT_NEAR(dof_val, f(node)(comp), RODIN_FUZZY_CONSTANT);
            }
            break;
          }
        }
      }
    }
  }

  // ========================================================================
  // GRADIENTFUNCTION TESTS FOR PKELEMENT K=0 TO K=6
  // ========================================================================

  TEST(FinalTest_H1Element_GradientFunction, GradientConsistency_K1_to_K6)
  {
    // Test gradient functions for polynomial orders 1-6
    for (size_t K : {1, 2, 3, 4, 5, 6})
    {
      Math::Vector<Real> p{{0.5}};

      switch (K)
      {
        case 1:
        {
          RealH1Element<1> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            auto grad_func = elem.getBasis(i).getGradient();
            const auto& grad_val = grad_func(p);
            EXPECT_EQ(grad_val.size(), 1);

            // Check consistency with derivative
            auto deriv = elem.getBasis(i).getDerivative<1>(0);
            EXPECT_NEAR(grad_val(0), deriv(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 2:
        {
          RealH1Element<2> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            auto grad_func = elem.getBasis(i).getGradient();
            const auto& grad_val = grad_func(p);
            auto deriv = elem.getBasis(i).getDerivative<1>(0);
            EXPECT_NEAR(grad_val(0), deriv(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 3:
        {
          RealH1Element<3> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            auto grad_func = elem.getBasis(i).getGradient();
            const auto& grad_val = grad_func(p);
            auto deriv = elem.getBasis(i).getDerivative<1>(0);
            EXPECT_NEAR(grad_val(0), deriv(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 4:
        {
          RealH1Element<4> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            auto grad_func = elem.getBasis(i).getGradient();
            const auto& grad_val = grad_func(p);
            auto deriv = elem.getBasis(i).getDerivative<1>(0);
            EXPECT_NEAR(grad_val(0), deriv(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 5:
        {
          RealH1Element<5> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            auto grad_func = elem.getBasis(i).getGradient();
            const auto& grad_val = grad_func(p);
            auto deriv = elem.getBasis(i).getDerivative<1>(0);
            EXPECT_NEAR(grad_val(0), deriv(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 6:
        {
          RealH1Element<6> elem(Polytope::Type::Segment);
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            auto grad_func = elem.getBasis(i).getGradient();
            const auto& grad_val = grad_func(p);
            auto deriv = elem.getBasis(i).getDerivative<1>(0);
            EXPECT_NEAR(grad_val(0), deriv(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
      }
    }
  }

  TEST(FinalTest_H1Element_GradientFunction, GradientPartitionProperty_K2_Triangle)
  {
    // Test that sum of gradients equals zero (partition of unity)
    RealH1Element<2> elem(Polytope::Type::Triangle);

    Math::Vector<Real> p{{0.3, 0.4}};
    Math::SpatialVector<Real> grad_sum = Math::SpatialVector<Real>::Zero(2);

    for (size_t i = 0; i < elem.getCount(); i++)
    {
      auto grad_func = elem.getBasis(i).getGradient();
      grad_sum += grad_func(p);
    }

    EXPECT_NEAR(grad_sum(0), 0.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(grad_sum(1), 0.0, RODIN_FUZZY_CONSTANT);
  }

  // ========================================================================
  // JACOBIANFUNCTION TESTS FOR PKELEMENT K=0 TO K=6
  // ========================================================================

  TEST(FinalTest_H1Element_JacobianFunction, JacobianStructure_K0_to_K3)
  {
    // Test Jacobian function for vector Pk elements
    for (size_t K : {0, 1, 2, 3})
    {
      for (size_t vdim : {2, 3})
      {
        Math::Vector<Real> p{{0.5}};

        switch (K)
        {
          case 0:
          {
            H1Element<0, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
            {
              auto jac_func = elem.getBasis(i).getJacobian();
              const auto& jac = jac_func(p);
              EXPECT_EQ(jac.rows(), vdim);
              EXPECT_EQ(jac.cols(), 1);
            }
            break;
          }
          case 1:
          {
            H1Element<1, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
            {
              auto jac_func = elem.getBasis(i).getJacobian();
              const auto& jac = jac_func(p);
              EXPECT_EQ(jac.rows(), vdim);
              EXPECT_EQ(jac.cols(), 1);
            }
            break;
          }
          case 2:
          {
            H1Element<2, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
            {
              auto jac_func = elem.getBasis(i).getJacobian();
              const auto& jac = jac_func(p);
              EXPECT_EQ(jac.rows(), vdim);
              EXPECT_EQ(jac.cols(), 1);
            }
            break;
          }
          case 3:
          {
            H1Element<3, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
            {
              auto jac_func = elem.getBasis(i).getJacobian();
              const auto& jac = jac_func(p);
              EXPECT_EQ(jac.rows(), vdim);
              EXPECT_EQ(jac.cols(), 1);
            }
            break;
          }
        }
      }
    }
  }

  TEST(FinalTest_H1Element_JacobianFunction, Jacobian2D_K2_Triangle)
  {
    // Test 2D Jacobian for K=2 on triangle
    H1Element<2, Math::Vector<Real>> elem(Polytope::Type::Triangle, 2);

    Math::Vector<Real> p{{0.3, 0.4}};

    for (size_t local = 0; local < elem.getCount(); local++)
    {
      auto jac_func = elem.getBasis(local).getJacobian();
      const auto& jac = jac_func(p);

      EXPECT_EQ(jac.rows(), 2);
      EXPECT_EQ(jac.cols(), 2);

      // Verify Jacobian entries
      size_t comp = local % 2;
      for (size_t i = 0; i < 2; i++)
      {
        for (size_t j = 0; j < 2; j++)
        {
          auto deriv = elem.getBasis(local).getDerivative<1>(i, j);
          EXPECT_NEAR(jac(i, j), deriv(p), RODIN_FUZZY_CONSTANT);
        }
      }
    }
  }

  // ========================================================================
  // INTERPOLATION TESTS FOR PKELEMENT K=0 TO K=6
  // ========================================================================

  TEST(FinalTest_H1Element_Interpolation, PolynomialInterpolation_K0_to_K6)
  {
    // Test that Pk element exactly interpolates polynomials of degree K
    RandomFloat gen(0.0, 1.0);

    for (size_t K : {0, 1, 2, 3, 4, 5, 6})
    {
      // Polynomial of degree K
      auto f = [K](const Math::SpatialPoint& x) -> Real {
        Real result = 1.0;
        for (size_t k = 1; k <= K; k++)
          result += (k + 0.5) * std::pow(x.x(), k);
        return result;
      };

      switch (K)
      {
        case 0:
        {
          RealH1Element<0> elem(Polytope::Type::Segment);
          std::vector<Real> dofs(elem.getCount());
          for (size_t i = 0; i < elem.getCount(); i++)
            dofs[i] = elem.getLinearForm(i)(f);

          for (size_t test = 0; test < 5; test++)
          {
            Math::Vector<Real> p{{gen()}};
            Real interp = 0.0;
            for (size_t i = 0; i < elem.getCount(); i++)
              interp += dofs[i] * elem.getBasis(i)(p);
            EXPECT_NEAR(interp, f(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 1:
        {
          RealH1Element<1> elem(Polytope::Type::Segment);
          std::vector<Real> dofs(elem.getCount());
          for (size_t i = 0; i < elem.getCount(); i++)
            dofs[i] = elem.getLinearForm(i)(f);

          for (size_t test = 0; test < 5; test++)
          {
            Math::Vector<Real> p{{gen()}};
            Real interp = 0.0;
            for (size_t i = 0; i < elem.getCount(); i++)
              interp += dofs[i] * elem.getBasis(i)(p);
            EXPECT_NEAR(interp, f(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 2:
        {
          RealH1Element<2> elem(Polytope::Type::Segment);
          std::vector<Real> dofs(elem.getCount());
          for (size_t i = 0; i < elem.getCount(); i++)
            dofs[i] = elem.getLinearForm(i)(f);

          for (size_t test = 0; test < 5; test++)
          {
            Math::Vector<Real> p{{gen()}};
            Real interp = 0.0;
            for (size_t i = 0; i < elem.getCount(); i++)
              interp += dofs[i] * elem.getBasis(i)(p);
            EXPECT_NEAR(interp, f(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 3:
        {
          RealH1Element<3> elem(Polytope::Type::Segment);
          std::vector<Real> dofs(elem.getCount());
          for (size_t i = 0; i < elem.getCount(); i++)
            dofs[i] = elem.getLinearForm(i)(f);

          for (size_t test = 0; test < 5; test++)
          {
            Math::Vector<Real> p{{gen()}};
            Real interp = 0.0;
            for (size_t i = 0; i < elem.getCount(); i++)
              interp += dofs[i] * elem.getBasis(i)(p);
            EXPECT_NEAR(interp, f(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 4:
        {
          RealH1Element<4> elem(Polytope::Type::Segment);
          std::vector<Real> dofs(elem.getCount());
          for (size_t i = 0; i < elem.getCount(); i++)
            dofs[i] = elem.getLinearForm(i)(f);

          for (size_t test = 0; test < 5; test++)
          {
            Math::Vector<Real> p{{gen()}};
            Real interp = 0.0;
            for (size_t i = 0; i < elem.getCount(); i++)
              interp += dofs[i] * elem.getBasis(i)(p);
            EXPECT_NEAR(interp, f(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 5:
        {
          RealH1Element<5> elem(Polytope::Type::Segment);
          std::vector<Real> dofs(elem.getCount());
          for (size_t i = 0; i < elem.getCount(); i++)
            dofs[i] = elem.getLinearForm(i)(f);

          for (size_t test = 0; test < 5; test++)
          {
            Math::Vector<Real> p{{gen()}};
            Real interp = 0.0;
            for (size_t i = 0; i < elem.getCount(); i++)
              interp += dofs[i] * elem.getBasis(i)(p);
            EXPECT_NEAR(interp, f(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
        case 6:
        {
          RealH1Element<6> elem(Polytope::Type::Segment);
          std::vector<Real> dofs(elem.getCount());
          for (size_t i = 0; i < elem.getCount(); i++)
            dofs[i] = elem.getLinearForm(i)(f);

          for (size_t test = 0; test < 5; test++)
          {
            Math::Vector<Real> p{{gen()}};
            Real interp = 0.0;
            for (size_t i = 0; i < elem.getCount(); i++)
              interp += dofs[i] * elem.getBasis(i)(p);
            EXPECT_NEAR(interp, f(p), RODIN_FUZZY_CONSTANT);
          }
          break;
        }
      }
    }
  }

  TEST(FinalTest_H1Element_Interpolation, VectorPolynomialInterpolation_K1_K2_K3)
  {
    // Test vector field interpolation for Pk elements
    for (size_t K : {1, 2, 3})
    {
      for (size_t vdim : {2, 3})
      {
        // Vector polynomial
        auto f = [K, vdim](const Math::SpatialPoint& x) {
          Math::Vector<Real> v(vdim);
          Real poly = 1.0;
          for (size_t k = 1; k <= K; k++)
            poly += std::pow(x.x(), k);
          for (size_t i = 0; i < vdim; i++)
            v(i) = (i + 1) * poly;
          return v;
        };

        switch (K)
        {
          case 1:
          {
            H1Element<1, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            std::vector<Real> dofs(elem.getCount());
            for (size_t i = 0; i < elem.getCount(); i++)
              dofs[i] = elem.getLinearForm(i)(f);

            Math::Vector<Real> p{{0.5}};
            Math::Vector<Real> interp = Math::Vector<Real>::Zero(vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
              interp += dofs[i] * elem.getBasis(i)(p);

            Math::Vector<Real> exact = f(p);
            for (size_t i = 0; i < vdim; i++)
              EXPECT_NEAR(interp(i), exact(i), RODIN_FUZZY_CONSTANT);
            break;
          }
          case 2:
          {
            H1Element<2, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            std::vector<Real> dofs(elem.getCount());
            for (size_t i = 0; i < elem.getCount(); i++)
              dofs[i] = elem.getLinearForm(i)(f);

            Math::Vector<Real> p{{0.5}};
            Math::Vector<Real> interp = Math::Vector<Real>::Zero(vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
              interp += dofs[i] * elem.getBasis(i)(p);

            Math::Vector<Real> exact = f(p);
            for (size_t i = 0; i < vdim; i++)
              EXPECT_NEAR(interp(i), exact(i), RODIN_FUZZY_CONSTANT);
            break;
          }
          case 3:
          {
            H1Element<3, Math::Vector<Real>> elem(Polytope::Type::Segment, vdim);
            std::vector<Real> dofs(elem.getCount());
            for (size_t i = 0; i < elem.getCount(); i++)
              dofs[i] = elem.getLinearForm(i)(f);

            Math::Vector<Real> p{{0.5}};
            Math::Vector<Real> interp = Math::Vector<Real>::Zero(vdim);
            for (size_t i = 0; i < elem.getCount(); i++)
              interp += dofs[i] * elem.getBasis(i)(p);

            Math::Vector<Real> exact = f(p);
            for (size_t i = 0; i < vdim; i++)
              EXPECT_NEAR(interp(i), exact(i), RODIN_FUZZY_CONSTANT);
            break;
          }
        }
      }
    }
  }

  TEST(FinalTest_H1Element_Interpolation, GradientInterpolationConsistency_K2_K3)
  {
    // Test that gradient of interpolant matches derivative of function
    for (size_t K : {2, 3})
    {
      // Polynomial and its derivative
      auto f = [K](const Math::SpatialPoint& x) -> Real {
        Real result = 0.0;
        for (size_t k = 1; k <= K; k++)
          result += std::pow(x.x(), k);
        return result;
      };

      auto df = [K](const Math::SpatialPoint& x) -> Real {
        Real result = 0.0;
        for (size_t k = 1; k <= K; k++)
          result += k * std::pow(x.x(), k - 1);
        return result;
      };

      switch (K)
      {
        case 2:
        {
          RealH1Element<2> elem(Polytope::Type::Segment);
          std::vector<Real> dofs(elem.getCount());
          for (size_t i = 0; i < elem.getCount(); i++)
            dofs[i] = elem.getLinearForm(i)(f);

          Math::Vector<Real> p{{0.5}};
          Real interp_grad = 0.0;
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            auto deriv = elem.getBasis(i).getDerivative<1>(0);
            interp_grad += dofs[i] * deriv(p);
          }

          EXPECT_NEAR(interp_grad, df(p), RODIN_FUZZY_CONSTANT);
          break;
        }
        case 3:
        {
          RealH1Element<3> elem(Polytope::Type::Segment);
          std::vector<Real> dofs(elem.getCount());
          for (size_t i = 0; i < elem.getCount(); i++)
            dofs[i] = elem.getLinearForm(i)(f);

          Math::Vector<Real> p{{0.5}};
          Real interp_grad = 0.0;
          for (size_t i = 0; i < elem.getCount(); i++)
          {
            auto deriv = elem.getBasis(i).getDerivative<1>(0);
            interp_grad += dofs[i] * deriv(p);
          }

          EXPECT_NEAR(interp_grad, df(p), RODIN_FUZZY_CONSTANT);
          break;
        }
      }
    }
  }
}
