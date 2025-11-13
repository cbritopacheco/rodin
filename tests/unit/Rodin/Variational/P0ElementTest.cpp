#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include <complex>
#include "Rodin/Variational/P0.h"

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
}
