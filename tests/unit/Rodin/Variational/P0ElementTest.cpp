#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include <complex>
#include <functional>
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

  // ========================================================================
  // NEW COMPREHENSIVE TESTS FOR P0ELEMENT VECTOR DIMENSIONS
  // ========================================================================

  TEST(FinalTest_P0Element_Vector, VectorDimensions_1D_2D_3D_AllGeometries)
  {
    // Test vdim=1, 2, 3 on Segment
    {
      for (size_t vdim : {1, 2, 3})
      {
        VectorP0Element<Real> elem(Polytope::Type::Segment, vdim);
        EXPECT_EQ(elem.getCount(), vdim);
        EXPECT_EQ(elem.getOrder(), 0);

        Math::Vector<Real> p{{0.5}};
        for (size_t i = 0; i < vdim; i++)
        {
          const auto& val = elem.getBasis(i)(p);
          EXPECT_EQ(val.size(), vdim);
          for (size_t j = 0; j < vdim; j++)
          {
            if (i == j)
              EXPECT_NEAR(val(j), 1.0, RODIN_FUZZY_CONSTANT);
            else
              EXPECT_NEAR(val(j), 0.0, RODIN_FUZZY_CONSTANT);
          }
        }
      }
    }

    // Test vdim=1, 2, 3 on Triangle
    {
      for (size_t vdim : {1, 2, 3})
      {
        VectorP0Element<Real> elem(Polytope::Type::Triangle, vdim);
        EXPECT_EQ(elem.getCount(), vdim);
      }
    }

    // Test vdim=1, 2, 3 on Quadrilateral
    {
      for (size_t vdim : {1, 2, 3})
      {
        VectorP0Element<Real> elem(Polytope::Type::Quadrilateral, vdim);
        EXPECT_EQ(elem.getCount(), vdim);
      }
    }

    // Test vdim=1, 2, 3 on Tetrahedron
    {
      for (size_t vdim : {1, 2, 3})
      {
        VectorP0Element<Real> elem(Polytope::Type::Tetrahedron, vdim);
        EXPECT_EQ(elem.getCount(), vdim);
      }
    }

    // Test vdim=1, 2, 3 on Wedge
    {
      for (size_t vdim : {1, 2, 3})
      {
        VectorP0Element<Real> elem(Polytope::Type::Wedge, vdim);
        EXPECT_EQ(elem.getCount(), vdim);
      }
    }
  }

  TEST(FinalTest_P0Element_Vector, PartitionOfUnity_AllVectorDimensions)
  {
    for (auto geom : {Polytope::Type::Segment, Polytope::Type::Triangle,
                      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
                      Polytope::Type::Wedge})
    {
      for (size_t vdim : {1, 2, 3})
      {
        VectorP0Element<Real> elem(geom, vdim);

        // Create appropriate test point based on geometry
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

        // Each component should have partition of unity
        std::vector<Real> sum_per_component(vdim, 0.0);
        for (size_t i = 0; i < elem.getCount(); i++)
        {
          const auto& val = elem.getBasis(i)(p);
          for (size_t j = 0; j < vdim; j++)
            sum_per_component[j] += val(j);
        }

        for (size_t j = 0; j < vdim; j++)
          EXPECT_NEAR(sum_per_component[j], 1.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }

  // ========================================================================
  // LINEARFORM TESTS FOR P0ELEMENT
  // ========================================================================

  TEST(FinalTest_P0Element_LinearForm, LinearForm_AllGeometries)
  {
    // Test LinearForm evaluation for P0 elements across all geometries
    for (auto geom : {Polytope::Type::Point, Polytope::Type::Segment,
                      Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
                      Polytope::Type::Tetrahedron, Polytope::Type::Wedge})
    {
      RealP0Element elem(geom);

      // P0 element has only one DOF - test the linear form
      const auto& lf = elem.getLinearForm(0);

      // Create a test function that returns a constant value
      const Real constant_value = 5.0;
      auto test_func = [constant_value](const Math::SpatialPoint&) { return constant_value; };

      // LinearForm should evaluate the function at the barycenter
      Real result = lf(test_func);
      EXPECT_NEAR(result, 5.0, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(FinalTest_P0Element_LinearForm, VectorLinearForm_AllVectorDimensions)
  {
    // Test LinearForm for vector P0 elements with different vector dimensions
    for (auto geom : {Polytope::Type::Segment, Polytope::Type::Triangle,
                      Polytope::Type::Tetrahedron})
    {
      for (size_t vdim : {1, 2, 3})
      {
        VectorP0Element<Real> elem(geom, vdim);

        // Test each component's linear form
        for (size_t i = 0; i < vdim; i++)
        {
          decltype(auto) lf = elem.getLinearForm(i);

          // Create a test vector function
          const size_t captured_vdim = vdim;
          const size_t captured_i = i;
          const Real test_value = 3.0 + static_cast<Real>(i);
          std::function<Math::Vector<Real>(const Math::SpatialPoint&)> test_func =
            [captured_vdim, captured_i, test_value](const Math::SpatialPoint&) -> Math::Vector<Real> {
              Math::Vector<Real> v;
              v.resize(captured_vdim);
              for (size_t j = 0; j < captured_vdim; j++)
                v(j) = (j == captured_i) ? test_value : 0.0;
              return v;
            };

          Real result = lf(test_func);
          EXPECT_NEAR(result, test_value, RODIN_FUZZY_CONSTANT);
        }
      }
    }
  }

  // ========================================================================
  // INTERPOLATION TESTS FOR P0ELEMENT
  // ========================================================================

  TEST(FinalTest_P0Element_Interpolation, ScalarConstantInterpolation)
  {
    // Test that P0 element correctly interpolates constant functions
    for (auto geom : {Polytope::Type::Segment, Polytope::Type::Triangle,
                      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
                      Polytope::Type::Wedge})
    {
      RealP0Element elem(geom);

      // Constant function f(x) = 7.5
      const Real constant_value = 7.5;
      auto f = [constant_value](const Math::SpatialPoint&) { return constant_value; };

      // Get DOF value using linear form
      Real dof_value = elem.getLinearForm(0)(f);

      // Interpolate at various points
      Math::Vector<Real> p;
      switch (geom)
      {
        case Polytope::Type::Segment:
          p = Math::Vector<Real>{{0.3}};
          break;
        case Polytope::Type::Triangle:
        case Polytope::Type::Quadrilateral:
          p = Math::Vector<Real>{{0.4, 0.3}};
          break;
        case Polytope::Type::Tetrahedron:
        case Polytope::Type::Wedge:
          p = Math::Vector<Real>{{0.2, 0.3, 0.4}};
          break;
        default:
          continue;
      }

      Real interpolated = dof_value * elem.getBasis(0)(p);
      EXPECT_NEAR(interpolated, 7.5, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(FinalTest_P0Element_Interpolation, VectorConstantInterpolation)
  {
    // Test vector P0 element interpolation for constant vector fields
    for (size_t vdim : {1, 2, 3})
    {
      VectorP0Element<Real> elem(Polytope::Type::Triangle, vdim);

      // Constant vector field
      const size_t captured_vdim = vdim;
      std::function<Math::Vector<Real>(const Math::SpatialPoint&)> f = 
        [captured_vdim](const Math::SpatialPoint&) -> Math::Vector<Real> {
          Math::Vector<Real> v;
          v.resize(captured_vdim);
          for (size_t i = 0; i < captured_vdim; i++)
            v(i) = 2.0 + 0.5 * static_cast<Real>(i);
          return v;
        };

      // Get DOF values
      std::vector<Real> dof_values(vdim);
      for (size_t i = 0; i < vdim; i++)
        dof_values[i] = elem.getLinearForm(i)(f);

      // Interpolate
      Math::Vector<Real> p{{0.3, 0.4}};
      Math::Vector<Real> interpolated = Math::Vector<Real>::Zero(vdim);
      for (size_t i = 0; i < vdim; i++)
        interpolated += dof_values[i] * elem.getBasis(i)(p);
      ASSERT_EQ(interpolated.size(), vdim);

      // Check each component
      for (size_t i = 0; i < vdim; i++)
        EXPECT_NEAR(interpolated(i), 2.0 + 0.5 * i, RODIN_FUZZY_CONSTANT);
    }
  }

  TEST(FinalTest_P0Element_Interpolation, InterpolationAccuracy_MultiplePoints)
  {
    // Test interpolation accuracy at multiple points for P0 element
    RealP0Element elem(Polytope::Type::Segment);

    // Constant function
    const Real constant_value = 3.14;
    auto f = [constant_value](const Math::SpatialPoint&) { return constant_value; };
    Real dof_value = elem.getLinearForm(0)(f);

    // Test at multiple points - should all give same value
    RandomFloat gen(0.0, 1.0);
    for (size_t i = 0; i < 10; i++)
    {
      Math::Vector<Real> p{{gen()}};
      Real interpolated = dof_value * elem.getBasis(0)(p);
      EXPECT_NEAR(interpolated, 3.14, RODIN_FUZZY_CONSTANT);
    }
  }
}
