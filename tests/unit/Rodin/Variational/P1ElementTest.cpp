#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include <complex>
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

    Complex f0(1.0, 2.0);  // x=0
    Complex f1(4.0, 6.0);  // x=1

    Math::Vector<Real> p{{0.6}};         // if reference segment is [0,1]
    Complex exact(2.8, 4.4);             // correct value
    Complex interpolated = f0 * elem.getBasis(0)(p) + f1 * elem.getBasis(1)(p);

    EXPECT_NEAR(interpolated.real(), exact.real(), RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(interpolated.imag(), exact.imag(), RODIN_FUZZY_CONSTANT);
  }

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

  // ========================================================================
  // NEW COMPREHENSIVE TESTS FOR P1ELEMENT
  // ========================================================================

  // Test P1 element on Wedge geometry (3D prismatic element)
  TEST(Rodin_Variational_RealP1Element, SanityTest_3D_Reference_Wedge)
  {
    RealP1Element k(Polytope::Type::Wedge);

    // Wedge has 6 vertices, so 6 DOFs
    EXPECT_EQ(k.getCount(), 6);

    // Check order
    EXPECT_EQ(k.getOrder(), 2);  // Wedge is tensor product, so order 2

    // Test Lagrange property at vertices
    // Bottom triangle vertices: (0,0,0), (1,0,0), (0,1,0)
    // Top triangle vertices: (0,0,1), (1,0,1), (0,1,1)
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0, 0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(3)(Math::Vector<Real>{{0, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(4)(Math::Vector<Real>{{0, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(5)(Math::Vector<Real>{{0, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{1, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{1, 0, 0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{1, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(3)(Math::Vector<Real>{{1, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(4)(Math::Vector<Real>{{1, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(5)(Math::Vector<Real>{{1, 0, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }

    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 1, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0, 1, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0, 1, 0}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(3)(Math::Vector<Real>{{0, 1, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(4)(Math::Vector<Real>{{0, 1, 0}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(5)(Math::Vector<Real>{{0, 1, 0}}), 0, RODIN_FUZZY_CONSTANT);
    }

    // Test top triangle vertices
    {
      EXPECT_NEAR(k.getBasis(0)(Math::Vector<Real>{{0, 0, 1}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(1)(Math::Vector<Real>{{0, 0, 1}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(2)(Math::Vector<Real>{{0, 0, 1}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(3)(Math::Vector<Real>{{0, 0, 1}}), 1, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(4)(Math::Vector<Real>{{0, 0, 1}}), 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(k.getBasis(5)(Math::Vector<Real>{{0, 0, 1}}), 0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_RealP1Element, PartitionOfUnity_Wedge)
  {
    constexpr size_t n = 25;
    RandomFloat gen(0.0, 1.0);
    RealP1Element k(Polytope::Type::Wedge);

    for (size_t i = 0; i < n; i++)
    {
      // Generate random point inside wedge
      Real s = gen();
      Real t = gen() * (1 - s);  // Ensure s + t <= 1
      Real r = gen();
      Math::Vector<Real> p{{s, t, r}};

      Real sum = 0;
      for (size_t j = 0; j < k.getCount(); j++)
        sum += k.getBasis(j)(p);

      EXPECT_NEAR(sum, 1.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(Rodin_Variational_RealP1Element, DerivativeTest_3D_Reference_Wedge)
  {
    RealP1Element k(Polytope::Type::Wedge);
    Math::Vector<Real> p{{0.25, 0.25, 0.5}};

    // Test that derivatives are accessible and well-defined
    // Note: Wedge is a tensor product element, so derivatives are not
    // necessarily constant in all directions
    for (size_t i = 0; i < k.getCount(); i++)
    {
      auto deriv_x = k.getBasis(i).getDerivative<1>(0);
      auto deriv_y = k.getBasis(i).getDerivative<1>(1);
      auto deriv_z = k.getBasis(i).getDerivative<1>(2);

      // Just verify derivatives can be computed without error
      Real dx = deriv_x(p);
      Real dy = deriv_y(p);
      Real dz = deriv_z(p);

      // Derivatives should be finite
      EXPECT_TRUE(std::isfinite(dx));
      EXPECT_TRUE(std::isfinite(dy));
      EXPECT_TRUE(std::isfinite(dz));
    }
  }

  // Test vector P1 element with explicit vector dimensions
  TEST(FinalTest_P1Element_Vector, VectorDimensions_1D_2D_3D_Segment)
  {
    // Test vdim=1 (scalar-like)
    {
      VectorP1Element<Real> elem(Polytope::Type::Segment, 1);
      EXPECT_EQ(elem.getCount(), 2);  // 1 component × 2 nodes

      Math::Vector<Real> p{{0.5}};
      const auto& val = elem.getBasis(0)(p);
      EXPECT_EQ(val.size(), 1);
      EXPECT_NEAR(val(0), 0.5, RODIN_FUZZY_CONSTANT);
    }

    // Test vdim=2
    {
      VectorP1Element<Real> elem(Polytope::Type::Segment, 2);
      EXPECT_EQ(elem.getCount(), 4);  // 2 components × 2 nodes

      Math::Vector<Real> p{{0.5}};
      const auto& val0 = elem.getBasis(0)(p);  // First component, first node
      EXPECT_EQ(val0.size(), 2);
      EXPECT_NEAR(val0(0), 0.5, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(val0(1), 0.0, RODIN_FUZZY_CONSTANT);

      const auto& val1 = elem.getBasis(1)(p);  // Second component, first node
      EXPECT_NEAR(val1(0), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(val1(1), 0.5, RODIN_FUZZY_CONSTANT);
    }

    // Test vdim=3
    {
      VectorP1Element<Real> elem(Polytope::Type::Segment, 3);
      EXPECT_EQ(elem.getCount(), 6);  // 3 components × 2 nodes

      Math::Vector<Real> p{{0.5}};
      for (size_t i = 0; i < 3; i++)
      {
        const auto& val = elem.getBasis(i)(p);
        EXPECT_EQ(val.size(), 3);
        for (size_t j = 0; j < 3; j++)
        {
          if (i == j)
            EXPECT_NEAR(val(j), 0.5, RODIN_FUZZY_CONSTANT);
          else
            EXPECT_NEAR(val(j), 0.0, RODIN_FUZZY_CONSTANT);
        }
      }
    }
  }

  TEST(FinalTest_P1Element_Vector, VectorDimensions_Triangle)
  {
    // Test vdim=1 on Triangle
    {
      VectorP1Element<Real> elem(Polytope::Type::Triangle, 1);
      EXPECT_EQ(elem.getCount(), 3);  // 1 component × 3 nodes
    }

    // Test vdim=2 on Triangle
    {
      VectorP1Element<Real> elem(Polytope::Type::Triangle, 2);
      EXPECT_EQ(elem.getCount(), 6);  // 2 components × 3 nodes
    }

    // Test vdim=3 on Triangle
    {
      VectorP1Element<Real> elem(Polytope::Type::Triangle, 3);
      EXPECT_EQ(elem.getCount(), 9);  // 3 components × 3 nodes
    }
  }

  TEST(FinalTest_P1Element_Vector, VectorDimensions_Tetrahedron)
  {
    // Test vdim=3 on Tetrahedron (typical 3D case)
    VectorP1Element<Real> elem(Polytope::Type::Tetrahedron, 3);
    EXPECT_EQ(elem.getCount(), 12);  // 3 components × 4 nodes

    Math::Vector<Real> p{{0.25, 0.25, 0.25}};
    Real sum_per_component[3] = {0, 0, 0};

    for (size_t i = 0; i < elem.getCount(); i++)
    {
      const auto& val = elem.getBasis(i)(p);
      EXPECT_EQ(val.size(), 3);
      for (size_t j = 0; j < 3; j++)
        sum_per_component[j] += val(j);
    }

    // Each component should sum to 1 (partition of unity)
    for (size_t j = 0; j < 3; j++)
      EXPECT_NEAR(sum_per_component[j], 1.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(FinalTest_P1Element_Vector, VectorDimensions_Wedge)
  {
    // Test vdim=3 on Wedge
    VectorP1Element<Real> elem(Polytope::Type::Wedge, 3);
    EXPECT_EQ(elem.getCount(), 18);  // 3 components × 6 nodes
  }
}
