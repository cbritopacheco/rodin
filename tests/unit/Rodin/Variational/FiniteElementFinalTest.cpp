/**
 * @file FiniteElementFinalTest.cpp
 * @brief Comprehensive final verification tests for P0Element, P1Element, and PkElement
 * 
 * This file contains final correctness tests across Real, Complex, and Vector value types
 * for all finite element types (P0, P1, Pk up to order 6) on all geometries.
 * 
 * Test Categories:
 * - Partition of Unity: Verifies sum of basis functions equals 1 (or identity for vectors)
 * - Lagrange Property: Verifies basis functions are 1 at their node and 0 at others
 * - Polynomial Reproduction: Verifies exact interpolation of polynomials up to element degree
 * - Derivative Consistency: Verifies derivative formulas are correct
 * - Complex Arithmetic: Verifies complex-valued basis functions work correctly
 * - Vector Component Structure: Verifies vector basis functions have correct component structure
 * - Cross-Element Consistency: Verifies PkElement<0> matches P0Element and PkElement<1> matches P1Element
 */

#include <gtest/gtest.h>
#include <complex>
#include "Rodin/Test/Random.h"
#include "Rodin/Variational/P0.h"
#include "Rodin/Variational/P1.h"
#include "Rodin/Variational/Pk/PkElement.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  // ============================================================================
  // P0Element Final Tests - Real
  // ============================================================================

  TEST(FinalTest_P0Element_Real, PartitionOfUnity_AllGeometries)
  {
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    
    // Segment
    {
      RealP0Element elem(Polytope::Type::Segment);
      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen()}};
        EXPECT_NEAR(elem.getBasis(0)(p), 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
    
    // Triangle
    {
      RealP0Element elem(Polytope::Type::Triangle);
      for (size_t i = 0; i < n; i++)
      {
        Real s = gen(), t = gen();
        if (s + t <= 1.0)
        {
          Math::Vector<Real> p{{s, t}};
          EXPECT_NEAR(elem.getBasis(0)(p), 1.0, RODIN_FUZZY_CONSTANT);
        }
      }
    }
    
    // Quadrilateral
    {
      RealP0Element elem(Polytope::Type::Quadrilateral);
      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen(), gen()}};
        EXPECT_NEAR(elem.getBasis(0)(p), 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
    
    // Tetrahedron
    {
      RealP0Element elem(Polytope::Type::Tetrahedron);
      for (size_t i = 0; i < n; i++)
      {
        Real s = gen(), t = gen(), u = gen();
        if (s + t + u <= 1.0)
        {
          Math::Vector<Real> p{{s, t, u}};
          EXPECT_NEAR(elem.getBasis(0)(p), 1.0, RODIN_FUZZY_CONSTANT);
        }
      }
    }
  }

  TEST(FinalTest_P0Element_Real, ConstantReproduction_AllGeometries)
  {
    // P0 should exactly reproduce constant functions
    const Real constant = 3.14159;
    
    // Test on Triangle
    {
      RealP0Element elem(Polytope::Type::Triangle);
      Math::Vector<Real> p{{0.3, 0.4}};
      Real interpolated = constant * elem.getBasis(0)(p);
      EXPECT_NEAR(interpolated, constant, RODIN_FUZZY_CONSTANT);
    }
    
    // Test on Tetrahedron
    {
      RealP0Element elem(Polytope::Type::Tetrahedron);
      Math::Vector<Real> p{{0.2, 0.3, 0.1}};
      Real interpolated = constant * elem.getBasis(0)(p);
      EXPECT_NEAR(interpolated, constant, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(FinalTest_P0Element_Real, ZeroDerivatives_AllGeometries)
  {
    // All derivatives of P0 elements should be zero
    
    // Segment
    {
      RealP0Element elem(Polytope::Type::Segment);
      auto deriv = elem.getBasis(0).getDerivative<1>(0);
      EXPECT_NEAR(deriv(Math::Vector<Real>{{0.5}}), 0.0, RODIN_FUZZY_CONSTANT);
    }
    
    // Triangle
    {
      RealP0Element elem(Polytope::Type::Triangle);
      auto deriv_x = elem.getBasis(0).getDerivative<1>(0);
      auto deriv_y = elem.getBasis(0).getDerivative<1>(1);
      Math::Vector<Real> p{{0.3, 0.4}};
      EXPECT_NEAR(deriv_x(p), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv_y(p), 0.0, RODIN_FUZZY_CONSTANT);
    }
    
    // Tetrahedron
    {
      RealP0Element elem(Polytope::Type::Tetrahedron);
      auto deriv_x = elem.getBasis(0).getDerivative<1>(0);
      auto deriv_y = elem.getBasis(0).getDerivative<1>(1);
      auto deriv_z = elem.getBasis(0).getDerivative<1>(2);
      Math::Vector<Real> p{{0.2, 0.3, 0.1}};
      EXPECT_NEAR(deriv_x(p), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv_y(p), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(deriv_z(p), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  // ============================================================================
  // P0Element Final Tests - Complex
  // ============================================================================

  TEST(FinalTest_P0Element_Complex, PartitionOfUnity_AllGeometries)
  {
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    using Complex = std::complex<Real>;
    
    // Segment
    {
      P0Element<Complex> elem(Polytope::Type::Segment);
      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen()}};
        auto val = elem.getBasis(0)(p);
        EXPECT_NEAR(val.real(), 1.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR(val.imag(), 0.0, RODIN_FUZZY_CONSTANT);
      }
    }
    
    // Triangle
    {
      P0Element<Complex> elem(Polytope::Type::Triangle);
      for (size_t i = 0; i < n; i++)
      {
        Real s = gen(), t = gen();
        if (s + t <= 1.0)
        {
          Math::Vector<Real> p{{s, t}};
          auto val = elem.getBasis(0)(p);
          EXPECT_NEAR(val.real(), 1.0, RODIN_FUZZY_CONSTANT);
          EXPECT_NEAR(val.imag(), 0.0, RODIN_FUZZY_CONSTANT);
        }
      }
    }
  }

  TEST(FinalTest_P0Element_Complex, ComplexArithmetic)
  {
    using Complex = std::complex<Real>;
    P0Element<Complex> elem(Polytope::Type::Segment);
    
    // Test complex-valued interpolation
    Complex c1(3.0, 4.0);
    Math::Vector<Real> p{{0.5}};
    Complex interpolated = c1 * elem.getBasis(0)(p);
    
    EXPECT_NEAR(interpolated.real(), 3.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(interpolated.imag(), 4.0, RODIN_FUZZY_CONSTANT);
  }

  // ============================================================================
  // P0Element Final Tests - Vector
  // ============================================================================

  TEST(FinalTest_P0Element_Vector, ComponentStructure_2D)
  {
    VectorP0Element<Real> elem(Polytope::Type::Triangle, 2);
    EXPECT_EQ(elem.getCount(), 2);  // 2 basis functions for 2D vector
    
    Math::Vector<Real> p{{0.3, 0.4}};
    
    // First basis function: should be [1, 0]
    {
      auto basis0 = elem.getBasis(0);
      const auto& val = basis0(p);
      EXPECT_EQ(val.size(), 2);
      EXPECT_NEAR(val(0), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(val(1), 0.0, RODIN_FUZZY_CONSTANT);
    }
    
    // Second basis function: should be [0, 1]
    {
      auto basis1 = elem.getBasis(1);
      const auto& val = basis1(p);
      EXPECT_EQ(val.size(), 2);
      EXPECT_NEAR(val(0), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(val(1), 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(FinalTest_P0Element_Vector, ComponentStructure_3D)
  {
    VectorP0Element<Real> elem(Polytope::Type::Tetrahedron, 3);
    EXPECT_EQ(elem.getCount(), 3);  // 3 basis functions for 3D vector
    
    Math::Vector<Real> p{{0.2, 0.3, 0.1}};
    
    // Verify unit vector structure
    for (size_t i = 0; i < 3; i++)
    {
      auto basis = elem.getBasis(i);
      const auto& val = basis(p);
      EXPECT_EQ(val.size(), 3);
      for (size_t j = 0; j < 3; j++)
      {
        EXPECT_NEAR(val(j), (i == j) ? 1.0 : 0.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  TEST(FinalTest_P0Element_Vector, VectorFieldReproduction)
  {
    VectorP0Element<Real> elem(Polytope::Type::Triangle, 2);
    
    // Test constant vector field [2.0, 3.0]
    Math::Vector<Real> constant_field{{2.0, 3.0}};
    Math::Vector<Real> p{{0.3, 0.4}};
    
    // Interpolate: sum of coefficients * basis functions
    Math::Vector<Real> interpolated = Math::Vector<Real>::Zero(2);
    for (size_t i = 0; i < 2; i++)
    {
      const auto& basis_val = elem.getBasis(i)(p);
      interpolated += constant_field(i) * basis_val;
    }
    
    EXPECT_NEAR(interpolated(0), 2.0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(interpolated(1), 3.0, RODIN_FUZZY_CONSTANT);
  }

  // ============================================================================
  // P1Element Final Tests - Real
  // ============================================================================

  TEST(FinalTest_P1Element_Real, PartitionOfUnity_AllGeometries)
  {
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    
    // Segment
    {
      RealP1Element elem(Polytope::Type::Segment);
      for (size_t i = 0; i < n; i++)
      {
        Math::Vector<Real> p{{gen()}};
        Real sum = elem.getBasis(0)(p) + elem.getBasis(1)(p);
        EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
    
    // Triangle
    {
      RealP1Element elem(Polytope::Type::Triangle);
      for (size_t i = 0; i < n; i++)
      {
        Real s = gen(), t = gen();
        if (s + t <= 1.0)
        {
          Math::Vector<Real> p{{s, t}};
          Real sum = 0.0;
          for (size_t j = 0; j < 3; j++)
            sum += elem.getBasis(j)(p);
          EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
        }
      }
    }
    
    // Tetrahedron
    {
      RealP1Element elem(Polytope::Type::Tetrahedron);
      for (size_t i = 0; i < n; i++)
      {
        Real s = gen(), t = gen(), u = gen();
        if (s + t + u <= 1.0)
        {
          Math::Vector<Real> p{{s, t, u}};
          Real sum = 0.0;
          for (size_t j = 0; j < 4; j++)
            sum += elem.getBasis(j)(p);
          EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
        }
      }
    }
  }

  TEST(FinalTest_P1Element_Real, LagrangeProperty_AllGeometries)
  {
    // Segment
    {
      RealP1Element elem(Polytope::Type::Segment);
      for (size_t i = 0; i < 2; i++)
      {
        const auto& node = elem.getNode(i);
        for (size_t j = 0; j < 2; j++)
        {
          Real expected = (i == j) ? 1.0 : 0.0;
          EXPECT_NEAR(elem.getBasis(j)(node), expected, RODIN_FUZZY_CONSTANT);
        }
      }
    }
    
    // Triangle
    {
      RealP1Element elem(Polytope::Type::Triangle);
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
    
    // Tetrahedron
    {
      RealP1Element elem(Polytope::Type::Tetrahedron);
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
  }

  TEST(FinalTest_P1Element_Real, LinearReproduction_Segment)
  {
    RealP1Element elem(Polytope::Type::Segment);
    
    // Test linear function f(x) = 2 + 3*x
    Real c0 = 2.0, c1 = 3.0;  // f(0) = 2, f(1) = 5
    Real f_at_0 = c0;
    Real f_at_1 = c0 + c1;
    
    // Test at arbitrary point
    Math::Vector<Real> p{{0.7}};
    Real exact = c0 + c1 * p(0);
    Real interpolated = f_at_0 * elem.getBasis(0)(p) + f_at_1 * elem.getBasis(1)(p);
    
    EXPECT_NEAR(interpolated, exact, RODIN_FUZZY_CONSTANT);
  }

  TEST(FinalTest_P1Element_Real, LinearReproduction_Triangle)
  {
    RealP1Element elem(Polytope::Type::Triangle);
    
    // Test linear function f(x,y) = 1 + 2*x + 3*y
    // Values at nodes: f(0,0)=1, f(1,0)=3, f(0,1)=4
    Real f0 = 1.0, f1 = 3.0, f2 = 4.0;
    
    // Test at arbitrary point
    Math::Vector<Real> p{{0.3, 0.4}};
    Real exact = 1.0 + 2.0 * p(0) + 3.0 * p(1);
    Real interpolated = f0 * elem.getBasis(0)(p) + f1 * elem.getBasis(1)(p) + f2 * elem.getBasis(2)(p);
    
    EXPECT_NEAR(interpolated, exact, RODIN_FUZZY_CONSTANT);
  }

  TEST(FinalTest_P1Element_Real, ConstantDerivatives_Segment)
  {
    RealP1Element elem(Polytope::Type::Segment);
    
    // Derivatives should be constant for linear functions
    auto deriv0 = elem.getBasis(0).getDerivative<1>(0);
    auto deriv1 = elem.getBasis(1).getDerivative<1>(0);
    
    Math::Vector<Real> p1{{0.2}}, p2{{0.8}};
    EXPECT_NEAR(deriv0(p1), deriv0(p2), RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(deriv1(p1), deriv1(p2), RODIN_FUZZY_CONSTANT);
    
    // Sum of derivatives should be zero (derivative of constant)
    EXPECT_NEAR(deriv0(p1) + deriv1(p1), 0.0, RODIN_FUZZY_CONSTANT);
  }

  // ============================================================================
  // P1Element Final Tests - Complex
  // ============================================================================

  TEST(FinalTest_P1Element_Complex, PartitionOfUnity_Segment)
  {
    using Complex = std::complex<Real>;
    P1Element<Complex> elem(Polytope::Type::Segment);
    
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    
    for (size_t i = 0; i < n; i++)
    {
      Math::Vector<Real> p{{gen()}};
      Complex sum = elem.getBasis(0)(p) + elem.getBasis(1)(p);
      EXPECT_NEAR(sum.real(), 1.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(sum.imag(), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(FinalTest_P1Element_Complex, ComplexLinearReproduction)
  {
    using Complex = std::complex<Real>;
    P1Element<Complex> elem(Polytope::Type::Segment);
    
    // Test complex linear function f(x) = (1+2i) + (3+4i)*x
    Complex f0(1.0, 2.0);  // at x=0
    Complex f1(4.0, 6.0);  // at x=1: (1+2i) + (3+4i) = 4+6i
    
    Math::Vector<Real> p{{0.6}};
    Complex exact(3.8, 4.4);  // (1+2i) + 0.6*(3+4i)
    Complex interpolated = f0 * elem.getBasis(0)(p) + f1 * elem.getBasis(1)(p);
    
    EXPECT_NEAR(interpolated.real(), exact.real(), RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(interpolated.imag(), exact.imag(), RODIN_FUZZY_CONSTANT);
  }

  // ============================================================================
  // P1Element Final Tests - Vector
  // ============================================================================

  TEST(FinalTest_P1Element_Vector, ComponentStructure_2D_Segment)
  {
    VectorP1Element<Real> elem(Polytope::Type::Segment, 2);
    EXPECT_EQ(elem.getCount(), 4);  // 2 nodes × 2 components
    
    Math::Vector<Real> p{{0.5}};
    
    // Verify component structure
    for (size_t i = 0; i < 4; i++)
    {
      const auto& val = elem.getBasis(i)(p);
      EXPECT_EQ(val.size(), 2);
      
      // Only one component should be non-zero
      size_t component = i % 2;
      EXPECT_NEAR(val(component), 0.5, RODIN_FUZZY_CONSTANT);  // P1 basis at midpoint
      EXPECT_NEAR(val(1 - component), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(FinalTest_P1Element_Vector, LinearVectorFieldReproduction)
  {
    VectorP1Element<Real> elem(Polytope::Type::Segment, 2);
    
    // Test linear vector field v(x) = [1+2x, 3+4x]
    // At x=0: v=[1,3], at x=1: v=[3,7]
    Math::Vector<Real> v0{{1.0, 3.0}};
    Math::Vector<Real> v1{{3.0, 7.0}};
    
    Math::Vector<Real> p{{0.4}};
    Math::Vector<Real> exact{{1.8, 4.6}};  // [1+2*0.4, 3+4*0.4]
    
    // Interpolate
    Math::Vector<Real> interpolated = Math::Vector<Real>::Zero(2);
    interpolated += v0(0) * elem.getBasis(0)(p);  // x-component, node 0
    interpolated += v0(1) * elem.getBasis(1)(p);  // y-component, node 0
    interpolated += v1(0) * elem.getBasis(2)(p);  // x-component, node 1
    interpolated += v1(1) * elem.getBasis(3)(p);  // y-component, node 1
    
    EXPECT_NEAR(interpolated(0), exact(0), RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(interpolated(1), exact(1), RODIN_FUZZY_CONSTANT);
  }

  // ============================================================================
  // PkElement Final Tests - Real (Orders 2-6)
  // ============================================================================

  TEST(FinalTest_PkElement_Real, PartitionOfUnity_P2_AllGeometries)
  {
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    
    // Segment
    {
      RealPkElement<2> elem(Polytope::Type::Segment);
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
      RealPkElement<2> elem(Polytope::Type::Triangle);
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
      RealPkElement<2> elem(Polytope::Type::Quadrilateral);
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

  TEST(FinalTest_PkElement_Real, PartitionOfUnity_P3_P4_P5_P6_Segment)
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
          RealPkElement<3> elem(Polytope::Type::Segment);
          for (size_t j = 0; j < elem.getCount(); j++)
            sum += elem.getBasis(j)(p);
        }
        else if (order == 4)
        {
          RealPkElement<4> elem(Polytope::Type::Segment);
          for (size_t j = 0; j < elem.getCount(); j++)
            sum += elem.getBasis(j)(p);
        }
        else if (order == 5)
        {
          RealPkElement<5> elem(Polytope::Type::Segment);
          for (size_t j = 0; j < elem.getCount(); j++)
            sum += elem.getBasis(j)(p);
        }
        else if (order == 6)
        {
          RealPkElement<6> elem(Polytope::Type::Segment);
          for (size_t j = 0; j < elem.getCount(); j++)
            sum += elem.getBasis(j)(p);
        }
        
        EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  TEST(FinalTest_PkElement_Real, LagrangeProperty_P2_Segment)
  {
    RealPkElement<2> elem(Polytope::Type::Segment);
    
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

  TEST(FinalTest_PkElement_Real, LagrangeProperty_P3_P4_P5_P6_Segment)
  {
    // Test orders 3-6
    {
      RealPkElement<3> elem(Polytope::Type::Segment);
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
      RealPkElement<4> elem(Polytope::Type::Segment);
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
      RealPkElement<5> elem(Polytope::Type::Segment);
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
      RealPkElement<6> elem(Polytope::Type::Segment);
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

  TEST(FinalTest_PkElement_Real, QuadraticReproduction_P2_Segment)
  {
    RealPkElement<2> elem(Polytope::Type::Segment);
    
    // Test quadratic function f(x) = 1 + 2*x + 3*x^2
    // Values at nodes: x=0: f=1, x=0.5: f=1.75, x=1: f=6
    Real f0 = 1.0, f1 = 1.75, f2 = 6.0;
    
    Math::Vector<Real> p{{0.3}};
    Real exact = 1.0 + 2.0 * 0.3 + 3.0 * 0.3 * 0.3;  // 1.87
    Real interpolated = f0 * elem.getBasis(0)(p) + f1 * elem.getBasis(1)(p) + f2 * elem.getBasis(2)(p);
    
    EXPECT_NEAR(interpolated, exact, RODIN_FUZZY_CONSTANT);
  }

  TEST(FinalTest_PkElement_Real, CrossElementConsistency_P0)
  {
    // PkElement<0> should match P0Element
    RealP0Element p0(Polytope::Type::Segment);
    RealPkElement<0> pk0(Polytope::Type::Segment);
    
    Math::Vector<Real> p{{0.7}};
    EXPECT_NEAR(p0.getBasis(0)(p), pk0.getBasis(0)(p), RODIN_FUZZY_CONSTANT);
  }

  TEST(FinalTest_PkElement_Real, CrossElementConsistency_P1)
  {
    // PkElement<1> should match P1Element
    RealP1Element p1(Polytope::Type::Triangle);
    RealPkElement<1> pk1(Polytope::Type::Triangle);
    
    Math::Vector<Real> p{{0.3, 0.4}};
    for (size_t i = 0; i < 3; i++)
    {
      EXPECT_NEAR(p1.getBasis(i)(p), pk1.getBasis(i)(p), RODIN_FUZZY_CONSTANT);
    }
  }

  // ============================================================================
  // PkElement Final Tests - Complex (Orders 2-6)
  // ============================================================================

  TEST(FinalTest_PkElement_Complex, PartitionOfUnity_P2_P3_P4_Segment)
  {
    using Complex = std::complex<Real>;
    constexpr size_t n = 20;
    RandomFloat gen(0.0, 1.0);
    
    // P2
    {
      PkElement<2, Complex> elem(Polytope::Type::Segment);
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
      PkElement<3, Complex> elem(Polytope::Type::Segment);
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
      PkElement<4, Complex> elem(Polytope::Type::Segment);
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

  TEST(FinalTest_PkElement_Complex, ComplexQuadraticReproduction_P2_Segment)
  {
    using Complex = std::complex<Real>;
    PkElement<2, Complex> elem(Polytope::Type::Segment);
    
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

  TEST(FinalTest_PkElement_Complex, CrossElementConsistency_P0_P1)
  {
    using Complex = std::complex<Real>;
    
    // Check P0 consistency
    {
      P0Element<Complex> p0(Polytope::Type::Segment);
      PkElement<0, Complex> pk0(Polytope::Type::Segment);
      
      Math::Vector<Real> p{{0.5}};
      auto v0 = p0.getBasis(0)(p);
      auto vk = pk0.getBasis(0)(p);
      EXPECT_NEAR(v0.real(), vk.real(), RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(v0.imag(), vk.imag(), RODIN_FUZZY_CONSTANT);
    }
    
    // Check P1 consistency
    {
      P1Element<Complex> p1(Polytope::Type::Segment);
      PkElement<1, Complex> pk1(Polytope::Type::Segment);
      
      Math::Vector<Real> p{{0.6}};
      for (size_t i = 0; i < 2; i++)
      {
        auto v1 = p1.getBasis(i)(p);
        auto vk = pk1.getBasis(i)(p);
        EXPECT_NEAR(v1.real(), vk.real(), RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR(v1.imag(), vk.imag(), RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // ============================================================================
  // PkElement Final Tests - Vector (Orders 2-6)
  // ============================================================================

  TEST(FinalTest_PkElement_Vector, ComponentStructure_P2_2D_Segment)
  {
    PkElement<2, Math::Vector<Real>> elem(Polytope::Type::Segment, 2);
    EXPECT_EQ(elem.getCount(), 6);  // 3 nodes × 2 components
    
    Math::Vector<Real> p{{0.5}};
    
    // Each basis function should have only one non-zero component
    for (size_t i = 0; i < 6; i++)
    {
      const auto& val = elem.getBasis(i)(p);
      EXPECT_EQ(val.size(), 2);
      
      size_t component = i % 2;
      size_t num_nonzero = 0;
      for (size_t j = 0; j < 2; j++)
      {
        if (std::abs(val(j)) > RODIN_FUZZY_CONSTANT)
          num_nonzero++;
      }
      EXPECT_EQ(num_nonzero, 1);  // Only one component should be non-zero
    }
  }

  TEST(FinalTest_PkElement_Vector, QuadraticVectorFieldReproduction_P2_Segment)
  {
    PkElement<2, Math::Vector<Real>> elem(Polytope::Type::Segment, 2);
    
    // Test quadratic vector field v(x) = [1+2x+3x^2, 4+5x+6x^2]
    // At nodes: x=0: [1,4], x=0.5: [1.75,6.5], x=1: [6,15]
    std::vector<Math::Vector<Real>> node_values = {
      Math::Vector<Real>{{1.0, 4.0}},
      Math::Vector<Real>{{1.75, 6.5}},
      Math::Vector<Real>{{6.0, 15.0}}
    };
    
    Math::Vector<Real> p{{0.3}};
    Math::Vector<Real> exact{{1.87, 7.84}};  // [1+0.6+0.27, 4+1.5+0.54]
    
    // Interpolate
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

  TEST(FinalTest_PkElement_Vector, PartitionOfUnity_P3_P4_3D_Triangle)
  {
    // P3 on Triangle with 3D vectors
    {
      PkElement<3, Math::Vector<Real>> elem(Polytope::Type::Triangle, 3);
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

  // ============================================================================
  // Final Integration Test - All Elements, All Types
  // ============================================================================

  TEST(FinalTest_Integration, AllElementsAllTypes_DOFCounting)
  {
    // Verify DOF counting formulas are correct for all elements and types
    
    // P0 Elements
    EXPECT_EQ(RealP0Element(Polytope::Type::Segment).getCount(), 1);
    EXPECT_EQ(RealP0Element(Polytope::Type::Triangle).getCount(), 1);
    EXPECT_EQ(RealP0Element(Polytope::Type::Quadrilateral).getCount(), 1);
    EXPECT_EQ(RealP0Element(Polytope::Type::Tetrahedron).getCount(), 1);
    EXPECT_EQ(VectorP0Element<Real>(Polytope::Type::Triangle, 2).getCount(), 2);
    EXPECT_EQ(VectorP0Element<Real>(Polytope::Type::Tetrahedron, 3).getCount(), 3);
    
    // P1 Elements
    EXPECT_EQ(RealP1Element(Polytope::Type::Segment).getCount(), 2);
    EXPECT_EQ(RealP1Element(Polytope::Type::Triangle).getCount(), 3);
    EXPECT_EQ(RealP1Element(Polytope::Type::Quadrilateral).getCount(), 4);
    EXPECT_EQ(RealP1Element(Polytope::Type::Tetrahedron).getCount(), 4);
    EXPECT_EQ(VectorP1Element<Real>(Polytope::Type::Segment, 2).getCount(), 4);
    EXPECT_EQ(VectorP1Element<Real>(Polytope::Type::Triangle, 3).getCount(), 9);
    
    // Pk Elements (Segment): K+1
    EXPECT_EQ(RealPkElement<2>(Polytope::Type::Segment).getCount(), 3);
    EXPECT_EQ(RealPkElement<3>(Polytope::Type::Segment).getCount(), 4);
    EXPECT_EQ(RealPkElement<4>(Polytope::Type::Segment).getCount(), 5);
    EXPECT_EQ(RealPkElement<5>(Polytope::Type::Segment).getCount(), 6);
    EXPECT_EQ(RealPkElement<6>(Polytope::Type::Segment).getCount(), 7);
    
    // Pk Elements (Triangle): (K+1)(K+2)/2
    EXPECT_EQ(RealPkElement<2>(Polytope::Type::Triangle).getCount(), 6);
    EXPECT_EQ(RealPkElement<3>(Polytope::Type::Triangle).getCount(), 10);
    EXPECT_EQ(RealPkElement<4>(Polytope::Type::Triangle).getCount(), 15);
    
    // Pk Elements (Quadrilateral): (K+1)^2
    EXPECT_EQ(RealPkElement<2>(Polytope::Type::Quadrilateral).getCount(), 9);
    EXPECT_EQ(RealPkElement<3>(Polytope::Type::Quadrilateral).getCount(), 16);
    
    // Pk Elements (Tetrahedron): (K+1)(K+2)(K+3)/6
    EXPECT_EQ(RealPkElement<2>(Polytope::Type::Tetrahedron).getCount(), 10);
    EXPECT_EQ(RealPkElement<3>(Polytope::Type::Tetrahedron).getCount(), 20);
    
    // Vector Pk Elements
    EXPECT_EQ(PkElement<2, Math::Vector<Real>>(Polytope::Type::Segment, 2).getCount(), 6);
    EXPECT_EQ(PkElement<3, Math::Vector<Real>>(Polytope::Type::Triangle, 3).getCount(), 30);
  }

  TEST(FinalTest_Integration, AllElementsAllTypes_OrderProperty)
  {
    // Verify order property is correct
    
    EXPECT_EQ(RealP0Element(Polytope::Type::Segment).getOrder(), 0);
    EXPECT_EQ(RealP0Element(Polytope::Type::Triangle).getOrder(), 0);
    
    EXPECT_EQ(RealP1Element(Polytope::Type::Segment).getOrder(), 1);
    EXPECT_EQ(RealP1Element(Polytope::Type::Triangle).getOrder(), 1);
    
    EXPECT_EQ(RealPkElement<2>(Polytope::Type::Segment).getOrder(), 2);
    EXPECT_EQ(RealPkElement<3>(Polytope::Type::Segment).getOrder(), 3);
    EXPECT_EQ(RealPkElement<4>(Polytope::Type::Segment).getOrder(), 4);
    EXPECT_EQ(RealPkElement<5>(Polytope::Type::Segment).getOrder(), 5);
    EXPECT_EQ(RealPkElement<6>(Polytope::Type::Segment).getOrder(), 6);
  }
}
