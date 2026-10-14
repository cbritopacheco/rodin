/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Certifies cached H1 basis tables against direct and analytic values.
 */
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <vector>
#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Owns reference samples including boundaries without claiming exactness.
   *
   * Unit weights are irrelevant to basis tabulation. The inherited formula
   * identity distinguishes this sample collection from every other collection.
   */
  class H1SampleFormula final : public QF::QuadratureFormulaBase
  {
    public:
    /**
     * @brief Owns a collection of reference samples.
     * @param points Interior and boundary points to tabulate.
     */
    explicit H1SampleFormula(std::vector<Math::SpatialPoint> points)
      : m_points(std::move(points))
    {}

    /**
     * @brief Copies samples with a new formula identity.
     * @param other Sample collection to copy.
     */
    H1SampleFormula(const H1SampleFormula& other) = default;

    /**
     * @brief Gets the sample count.
     * @returns Number of owned reference samples.
     */
    size_t getSize() const override
    {
      return m_points.size();
    }

    /**
     * @brief Gets an unused unit weight.
     * @param index Sample index, unused because every weight is one.
     * @returns One, without any quadrature-exactness claim.
     */
    Real getWeight([[maybe_unused]] size_t index) const override
    {
      return 1;
    }

    /**
     * @brief Gets a reference sample.
     * @param index Sample index.
     * @returns The owned reference coordinates.
     */
    const Math::SpatialPoint& getPoint(size_t index) const override
    {
      return m_points.at(index);
    }

    /**
     * @brief Copies the owned samples with independent identity.
     * @returns An owned copy.
     */
    H1SampleFormula* copy() const noexcept override
    {
      return new H1SampleFormula(*this);
    }

    private:
    std::vector<Math::SpatialPoint> m_points;
  };

  /**
   * @brief Checks scalar H1 degrees one through six over both scalar fields.
   *
   * Cached entries must equal every direct basis value and first derivative
   * exactly, including vertices, faces and the historical pyramid-apex
   * convention. A separate polynomial reconstruction uses the public DOF
   * functionals at interior quadrature samples. It assumes neither a nodal
   * coefficient layout nor differentiability of rational bases at the apex.
   *
   * @tparam Element Scalar H1 element under test.
   */
  template <class Element>
  class H1TabulationTest : public ::testing::Test
  {};

  using H1TabulationCases =
    ::testing::Types<Variational::H1Element<1, Real>, Variational::H1Element<2, Real>,
      Variational::H1Element<3, Real>, Variational::H1Element<4, Real>,
      Variational::H1Element<5, Real>, Variational::H1Element<6, Real>,
      Variational::H1Element<1, Complex>, Variational::H1Element<2, Complex>,
      Variational::H1Element<3, Complex>, Variational::H1Element<4, Complex>,
      Variational::H1Element<5, Complex>, Variational::H1Element<6, Complex>>;
  TYPED_TEST_SUITE(H1TabulationTest, H1TabulationCases);

  TYPED_TEST(H1TabulationTest, CachedTableMatchesDirectEverySample)
  {
    using Element = TypeParam;
    using G = Geometry::Polytope::Type;
    for (const auto geometry : {G::Point, G::Segment, G::Triangle, G::Quadrilateral,
           G::Tetrahedron, G::Pyramid, G::Hexahedron, G::Wedge})
    {
      SCOPED_TRACE(int(geometry));
      Element element(geometry);
      const Geometry::Polytope::Traits traits(geometry);
      const size_t count = element.getCount(), dim = traits.getDimension();
      std::vector<Math::SpatialPoint> samples{traits.getCentroid()};
      for (size_t a = 0; a < traits.getVertexCount(); ++a)
        samples.push_back(traits.getVertex(a));
      const QF::PolytopeQuadratureFormula interior(8, geometry);
      for (size_t qp = 0; qp < interior.getSize(); ++qp)
        samples.push_back(interior.getPoint(qp));
      const H1SampleFormula formula(std::move(samples));
      const auto& table = element.getTabulation(formula);
      ASSERT_EQ(table.nqp, formula.getSize());
      ASSERT_EQ(table.ndof, count);
      ASSERT_EQ(table.dim, dim);
      ASSERT_EQ(table.phi.size(), formula.getSize() * count);
      ASSERT_EQ(table.dphi.size(), formula.getSize() * count * dim);
      EXPECT_EQ(&table, &element.getTabulation(formula));
      for (size_t qp = 0; qp < formula.getSize(); ++qp)
      {
        SCOPED_TRACE(qp);
        const auto& r = formula.getPoint(qp);
        for (size_t a = 0; a < count; ++a)
        {
          EXPECT_EQ(table.getBasis(qp, a), element.getBasis(a)(r));
          EXPECT_EQ(table.getGradient(qp, a).size(), dim);
          for (size_t i = 0; i < dim; ++i)
          {
            EXPECT_EQ(table.template getDerivative<1>(qp, a, i),
              element.getBasis(a).template getDerivative<1>(i)(r));
          }
        }
      }
    }
  }

  TYPED_TEST(H1TabulationTest, IndependentPolynomialReproduction)
  {
    using Element = TypeParam;
    using Scalar = typename Element::ScalarType;
    using G = Geometry::Polytope::Type;
    // Dimensionless roundoff budget for the independent polynomial route,
    // not a convergence-rate or performance threshold.
    constexpr Real ReproductionTolerance = 1e-8;
    const Scalar amplitude = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(Real(5) / 4, -Real(1) / 2);
      else
        return Scalar(Real(5) / 4);
    }();
    for (const auto geometry : {G::Point, G::Segment, G::Triangle, G::Quadrilateral,
           G::Tetrahedron, G::Pyramid, G::Hexahedron, G::Wedge})
    {
      SCOPED_TRACE(int(geometry));
      Element element(geometry);
      const size_t count = element.getCount();
      const size_t dim = Geometry::Polytope::Traits(geometry).getDimension();
      const size_t degree = geometry == G::Point ? 1 : element.getFactorOrder();
      const auto polynomial = [&](const Math::SpatialPoint& r) {
        Real result = 1;
        for (size_t j = 0; j < dim; ++j)
          result += Real(j + 1) * r(j) + std::pow(r(j), degree);
        return amplitude * result;
      };
      std::vector<Scalar> coefficients(count);
      for (size_t a = 0; a < count; ++a)
        coefficients[a] = element.getLinearForm(a)(polynomial);
      for (const size_t order : {2u, 4u, 8u})
      {
        const QF::PolytopeQuadratureFormula formula(order, geometry);
        const auto& table = element.getTabulation(formula);
        for (size_t qp = 0; qp < formula.getSize(); ++qp)
        {
          const auto& r = formula.getPoint(qp);
          Scalar value = 0;
          for (size_t a = 0; a < count; ++a)
            value += coefficients[a] * table.getBasis(qp, a);
          const Scalar expected = polynomial(r);
          const Real budget =
            ReproductionTolerance * std::max(Real(1), std::abs(expected));
          EXPECT_LE(std::abs(value - expected), budget);
          EXPECT_GT(std::abs(value + coefficients.front() - expected), budget);
          for (size_t i = 0; i < dim; ++i)
          {
            Scalar derivative = 0;
            for (size_t a = 0; a < count; ++a)
              derivative += coefficients[a] * table.template getDerivative<1>(qp, a, i);
            const Scalar expectedDerivative =
              amplitude * (Real(i + 1) + Real(degree) * std::pow(r(i), degree - 1));
            const Real derivativeBudget =
              ReproductionTolerance * std::max(Real(1), std::abs(expectedDerivative));
            EXPECT_LE(std::abs(derivative - expectedDerivative), derivativeBudget);
            EXPECT_GT(std::abs(derivative + coefficients.front() - expectedDerivative),
              derivativeBudget);
          }
        }
      }
    }
  }
}
