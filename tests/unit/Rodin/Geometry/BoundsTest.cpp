/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <functional>
#include <utility>
#include <vector>
#include <gtest/gtest.h>
#include <Rodin/Alert/Exception.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational/H1/H1Element.h>
#include <Rodin/Variational/Integrator.h>
#include <Rodin/Variational/P0/P0Element.h>
#include <Rodin/Variational/P0g/P0gElement.h>
#include <Rodin/Variational/P1/P1Element.h>
#include <Rodin/QF/QuadratureFormula.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /**
   * @brief Exercises support enclosures on arbitrary basis coefficients.
   * @tparam K H1 family index.
   * @param geometry Reference cell, including the rational pyramid.
   * @param translation Physical translation testing centered conversion.
   */
  template <size_t K>
  void checkBounds(Polytope::Type geometry, Real translation)
  {
    const H1Element<K, Real> element(geometry);
    PointCloud coefficients(3, element.getCount());
    for (size_t j = 0; j < element.getCount(); ++j)
    {
      for (size_t i = 0; i < 3; ++i)
        coefficients(i, j) = translation + std::sin(Real((i + 1) * (j + 2)));
    }
    const ParametricTransformation transformation(coefficients, element);
    const auto enclosure = transformation.getBounds();
    ASSERT_TRUE(enclosure);
    EXPECT_EQ(transformation.isFirstOrder(), K == 1 || geometry == Polytope::Type::Point);
    EXPECT_EQ(transformation.getDegree(), element.getDegree());
    const PolytopeTransformation& base = transformation;
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    for (size_t direction = 0; direction < dimension; ++direction)
      EXPECT_EQ(base.getDegree({direction}), element.getDegree({direction}));
    EXPECT_FALSE(base.getDegree({3}));
    if (geometry == Polytope::Type::Pyramid)
      EXPECT_FALSE(transformation.getDegree());
    const auto& formula = QF::QuadratureFormula::get(12, geometry);
    for (const auto& direction : {Math::SpatialVector<Real>{1, 0, 0},
           Math::SpatialVector<Real>{-1, 0, 0}, Math::SpatialVector<Real>{0, 1, 0},
           Math::SpatialVector<Real>{0, 0, -1}, Math::SpatialVector<Real>{2, -3, 4}})
    {
      const Real support = enclosure->getSupport(direction);
      ASSERT_TRUE(std::isfinite(support));
      for (size_t q = 0; q < formula.getSize(); ++q)
      {
        Math::SpatialPoint mapped;
        transformation.transform(mapped, formula.getPoint(q));
        EXPECT_LE(mapped.dot(direction), support);
      }
      const Polytope::Traits traits(geometry);
      for (size_t v = 0; v < traits.getVertexCount(); ++v)
      {
        Math::SpatialPoint mapped;
        transformation.transform(mapped, traits.getVertex(v));
        EXPECT_LE(mapped.dot(direction), support);
      }
    }
    EXPECT_EQ(enclosure->getSupport(Math::SpatialVector<Real>{0, 0, 0}), 0);
    EXPECT_THROW(enclosure->getSupport(Math::SpatialVector<Real>{1}), Alert::Exception);
  }
}

TEST(Rodin_Geometry_Bounds, PolynomialAndRationalImages)
{
  for (const auto geometry :
    {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge, Polytope::Type::Pyramid})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    for (const Real translation : {Real(0), Real(1e9)})
    {
      SCOPED_TRACE(translation);
      checkBounds<1>(geometry, translation);
      checkBounds<2>(geometry, translation);
      checkBounds<3>(geometry, translation);
      checkBounds<4>(geometry, translation);
    }
  }
}

TEST(Rodin_Variational_Degree, FamilyIndexAndPolynomialContainment)
{
  for (const auto geometry :
    {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge, Polytope::Type::Pyramid})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const auto check = [&]<size_t K>() {
      const H1Element<K, Real> scalar(geometry);
      const H1Element<K, Complex> complex(geometry);
      const H1Element<K, Math::SpatialVector<Real>> vector(geometry, 3);
      const H1Element<K, Math::SpatialMatrix<Complex>> matrix(geometry, 2, 3);
      EXPECT_EQ(scalar.getOrder(), K);
      EXPECT_EQ(complex.getOrder(), K);
      EXPECT_EQ(vector.getOrder(), K);
      EXPECT_EQ(matrix.getOrder(), K);
      Optional<size_t> degree;
      if (geometry == Polytope::Type::Point)
        degree = 0;
      else if (geometry == Polytope::Type::Hexahedron)
        degree = 3 * K;
      else if (geometry == Polytope::Type::Quadrilateral ||
        geometry == Polytope::Type::Wedge)
        degree = 2 * K;
      else if (geometry != Polytope::Type::Pyramid)
        degree = K;
      EXPECT_EQ(scalar.getDegree(), degree);
      EXPECT_EQ(complex.getDegree(), degree);
      EXPECT_EQ(vector.getDegree(), degree);
      EXPECT_EQ(matrix.getDegree(), degree);
      const FiniteElementBase<H1Element<K, Real>>& base = scalar;
      EXPECT_EQ(base.getDegree(), degree);
      EXPECT_EQ(base.getDegree({}), 0);
      EXPECT_FALSE(base.getDegree({3}));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      for (size_t direction = 0; direction < dimension; ++direction)
      {
        const Optional<size_t> expected =
          geometry == Polytope::Type::Pyramid && direction == 2 ? Optional<size_t>{}
                                                                : Optional<size_t>{K};
        EXPECT_EQ(scalar.getDegree({direction}), expected);
        EXPECT_EQ(vector.getDegree({direction}), expected);
        EXPECT_EQ(matrix.getDegree({direction}), expected);
        EXPECT_EQ(base.getDegree({direction}), expected);
        EXPECT_EQ(scalar.getDegree({direction, direction}), expected);
      }
      if (geometry == Polytope::Type::Quadrilateral ||
        geometry == Polytope::Type::Hexahedron || geometry == Polytope::Type::Pyramid)
        EXPECT_EQ(base.getDegree({0, 1}), 2 * K);
      else if (geometry == Polytope::Type::Wedge)
      {
        EXPECT_EQ(base.getDegree({0, 1}), K);
        EXPECT_EQ(base.getDegree({0, 2}), 2 * K);
      }
      else if (geometry == Polytope::Type::Triangle ||
        geometry == Polytope::Type::Tetrahedron)
        EXPECT_EQ(base.getDegree({0, 1}), K);
    };
    check.template operator()<1>();
    check.template operator()<2>();
    check.template operator()<3>();
    const P0Element<Real> discontinuous(geometry);
    const P0gElement<Real> globalConstant(geometry);
    EXPECT_EQ(discontinuous.getOrder(), 0);
    EXPECT_EQ(globalConstant.getOrder(), 0);
    EXPECT_EQ(discontinuous.getDegree({}), 0);
    EXPECT_EQ(globalConstant.getDegree({}), 0);
    if (Polytope::Traits(geometry).getDimension() > 0)
    {
      EXPECT_EQ(discontinuous.getDegree({0}), 0);
      EXPECT_EQ(globalConstant.getDegree({0}), 0);
    }
    EXPECT_FALSE(discontinuous.getDegree({3}));
    EXPECT_FALSE(globalConstant.getDegree({3}));
    const P1Element<Real> linear(geometry);
    EXPECT_EQ(linear.getOrder(), 1);
    EXPECT_EQ(linear.getDegree(), (H1Element<1, Real>(geometry).getDegree()));
    for (size_t direction = 0; direction < Polytope::Traits(geometry).getDimension();
         ++direction)
    {
      EXPECT_EQ(linear.getDegree({direction}),
        (H1Element<1, Real>(geometry).getDegree({direction})));
    }
  }
}

namespace
{
  /**
   * @brief Records the formula used while forwarding a polynomial map.
   * @tparam K Geometry family order.
   */
  template <size_t K>
  class CountingTransformation final : public PolytopeTransformation
  {
    public:
      /**
       * @brief Constructs a wrapper retaining an external sample counter.
       * @param map Parametric map used for every geometric operation.
       * @param samples Receives the size of the selected quadrature formula.
       */
      CountingTransformation(
        ParametricTransformation<H1Element<K, Real>> map, size_t& samples)
        : PolytopeTransformation(map.getReferenceDimension(), map.getPhysicalDimension()),
          m_map(std::move(map)),
          m_samples(samples)
      {}

      Optional<size_t> getDegree() const override
      {
        return m_map.getDegree();
      }

      Optional<size_t> getDegree(std::initializer_list<size_t> directions) const override
      {
        return m_map.getDegree(directions);
      }

      bool isFirstOrder() const override
      {
        return m_map.isFirstOrder();
      }

      void transform(
        Math::SpatialPoint& physical, const Math::SpatialPoint& reference) const override
      {
        m_map.transform(physical, reference);
      }

      void jacobian(Math::SpatialMatrix<Real>& matrix,
        const Math::SpatialPoint& reference) const override
      {
        m_map.jacobian(matrix, reference);
      }

      void jacobian(std::vector<Math::SpatialMatrix<Real>>& matrices,
        const QF::QuadratureFormula& formula) const override
      {
        m_samples.get() = formula.getSize();
        m_map.jacobian(matrices, formula);
      }

      CountingTransformation* copy() const noexcept override
      {
        return new CountingTransformation(*this);
      }

    private:
      ParametricTransformation<H1Element<K, Real>> m_map;
      std::reference_wrapper<size_t> m_samples;
  };
}

TEST(Rodin_Geometry_Degree, DirectionalMeasureSelectsTightProductRules)
{
  for (const auto geometry : {Polytope::Type::Quadrilateral, Polytope::Type::Hexahedron,
         Polytope::Type::Wedge, Polytope::Type::Pyramid})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    const auto cell = *mesh.getCell();
    const H1Element<2, Real> element(geometry);
    PointCloud coefficients(dimension, element.getCount());
    for (size_t a = 0; a < element.getCount(); ++a)
    {
      const auto point = element.getNode(a);
      Math::SpatialPoint physical;
      cell.getTransformation().transform(physical, point);
      physical[0] += Real(0.05) * point[0] * point[0];
      for (size_t j = 0; j < dimension; ++j)
        coefficients(j, a) = physical[j];
    }
    size_t samples = 0;
    mesh.setPolytopeTransformation({dimension, cell.getIndex()},
      new CountingTransformation<2>(
        ParametricTransformation(std::move(coefficients), element), samples));
    const Real measure = cell.getMeasure();
    const size_t expectedOrder =
      geometry == Polytope::Type::Pyramid ? 10 : dimension * 2 - 1;
    EXPECT_EQ(samples, QF::QuadratureFormula::get(expectedOrder, geometry).getSize());
    if (geometry == Polytope::Type::Hexahedron)
      EXPECT_EQ(samples, 27);
    EXPECT_NEAR(measure, cell.getMeasure(20), 1e-12);
  }
}

TEST(Rodin_Variational_Degree, PyramidHorizontalBoundsMatchReferenceBasis)
{
  const auto check = []<size_t K>() {
    const H1Element<K, Real> element(Polytope::Type::Pyramid);
    const Real z = Real(0.37);
    const Real width = 1 - z;
    for (size_t basis = 0; basis < element.getCount(); ++basis)
    {
      for (size_t axis = 0; axis < 2; ++axis)
      {
        Math::SpatialPoint point{Real(0.3) * width, Real(0.3) * width, z};
        for (const Real fraction : {Real(0.13), Real(0.41), Real(0.77)})
        {
          point[axis] = fraction * width;
          const Real actual = element.getBasis(basis)(point);
          Real interpolated = 0;
          for (size_t node = 0; node <= K; ++node)
          {
            const Real nodeFraction = Real(node) / Real(K);
            Real cardinal = 1;
            for (size_t other = 0; other <= K; ++other)
            {
              if (other != node)
              {
                cardinal *= (fraction - Real(other) / Real(K)) /
                  (nodeFraction - Real(other) / Real(K));
              }
            }
            point[axis] = nodeFraction * width;
            interpolated += cardinal * element.getBasis(basis)(point);
          }
          EXPECT_NEAR(interpolated, actual, 5e-10) << "K=" << K << " basis=" << basis;
        }
      }
    }
  };
  check.template operator()<1>();
  check.template operator()<2>();
  check.template operator()<3>();
  check.template operator()<4>();
}

TEST(Rodin_Variational_Degree, IntegrationFallbackUsesKnownStructure)
{
  struct UnknownElement
  {
    Optional<size_t> getDegree() const
    {
      return std::nullopt;
    }
  };
  for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Pyramid})
  {
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    const auto cell = *mesh.getCell();
    const H1Element<2, Real> element(geometry);
    EXPECT_EQ(
      Integrator::getDefaultOrder(UnknownElement{}, cell), Integrator::DefaultOrder);
    EXPECT_EQ(Integrator::getDefaultOrder(element, cell),
      geometry == Polytope::Type::Pyramid ? 4 : 2);
  }
}
