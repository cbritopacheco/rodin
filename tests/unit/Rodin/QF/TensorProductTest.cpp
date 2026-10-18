/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Certifies tensor-product ordering, ownership and invariant factor access.
 */
#include <gtest/gtest.h>
#include <vector>
#include "Rodin/QF/TensorProduct.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Counts factor queries without changing their immutable numerical data.
   *
   * These samples need not integrate a particular polynomial exactly: the
   * tested contract is Cartesian concatenation and ordered multiplication.
   */
  class ObservedProductFactor final : public QF::QuadratureFormulaBase
  {
    public:
      /**
       * @brief Owns samples and weights and initializes query counts.
       * @param points Factor coordinates.
       * @param weights Factor weights.
       */
      ObservedProductFactor(
        std::vector<Math::SpatialPoint> points, std::vector<Real> weights)
        : m_points(std::move(points)),
          m_weights(std::move(weights)),
          m_pointReads(0),
          m_weightReads(0),
          m_sizeReads(0)
      {}

      /// @brief Gets the factor size and records the query.
      size_t getSize() const override
      {
        ++m_sizeReads;
        return m_points.size();
      }

      /**
       * @brief Reads a weight and records the query.
       * @param index Factor sample index.
       * @returns Immutable factor weight.
       */
      Real getWeight(size_t index) const override
      {
        ++m_weightReads;
        return m_weights.at(index);
      }

      /**
       * @brief Reads coordinates and records the query.
       * @param index Factor sample index.
       * @returns Immutable factor coordinates.
       */
      const Math::SpatialPoint& getPoint(size_t index) const override
      {
        ++m_pointReads;
        return m_points.at(index);
      }

      /// @brief Copies immutable samples with fresh counters and identity.
      ObservedProductFactor* copy() const noexcept override
      {
        return new ObservedProductFactor(m_points, m_weights);
      }

      /// @brief Number of point queries since construction.
      size_t pointReads() const
      {
        return m_pointReads;
      }

      /// @brief Number of weight queries since construction.
      size_t weightReads() const
      {
        return m_weightReads;
      }

      /// @brief Number of size queries since construction.
      size_t sizeReads() const
      {
        return m_sizeReads;
      }

    private:
      std::vector<Math::SpatialPoint> m_points;
      std::vector<Real> m_weights;
      mutable size_t m_pointReads;
      mutable size_t m_weightReads;
      mutable size_t m_sizeReads;
  };

  /**
   * @brief Factor descriptors are read once at their invariant loop level.
   *
   * Zero-dimensional factors and unequal coordinate dimensions are included.
   * Product values and ordering are compared exactly, independently of the
   * production query counts. There is no wall-clock performance assertion.
   */
  TEST(TensorProduct, InvariantAccessAndExactConcatenation)
  {
    for (const auto& [ldim, rdim] :
      {std::pair{0u, 0u}, {0u, 1u}, {1u, 0u}, {1u, 1u}, {1u, 2u}, {2u, 1u}})
    {
      SCOPED_TRACE(ldim);
      SCOPED_TRACE(rdim);
      const size_t nl = ldim == 0 ? 1 : 2, nr = rdim == 0 ? 1 : 3;
      std::vector<Math::SpatialPoint> leftPoints, rightPoints;
      std::vector<Real> leftWeights, rightWeights;
      for (size_t i = 0; i < nl; ++i)
      {
        Math::SpatialPoint point(ldim);
        for (size_t j = 0; j < ldim; ++j)
          point(j) = Real(i + j + 1) / 8;
        leftPoints.push_back(point);
        leftWeights.push_back(Real(i + 1) / 7);
      }
      for (size_t i = 0; i < nr; ++i)
      {
        Math::SpatialPoint point(rdim);
        for (size_t j = 0; j < rdim; ++j)
          point(j) = Real(i + j + 1) / 9;
        rightPoints.push_back(point);
        rightWeights.push_back(Real(i + 2) / 11);
      }
      const ObservedProductFactor left(leftPoints, leftWeights);
      const ObservedProductFactor right(rightPoints, rightWeights);
      using G = Geometry::Polytope::Type;
      const auto geometry = ldim + rdim == 0 ? G::Point
        : ldim + rdim == 1                   ? G::Segment
        : ldim + rdim == 2                   ? G::Quadrilateral
                                             : G::Hexahedron;
      const QF::TensorProduct product(geometry, left, right);
      EXPECT_EQ(left.sizeReads(), 1u);
      EXPECT_EQ(right.sizeReads(), 1u);
      EXPECT_EQ(left.pointReads(), nl);
      EXPECT_EQ(right.pointReads(), nl * nr);
      EXPECT_EQ(left.weightReads(), nl);
      EXPECT_EQ(right.weightReads(), nl * nr);
      ASSERT_EQ(product.getSize(), nl * nr);
      EXPECT_EQ(product.getGeometry(), geometry);
      for (size_t i = 0; i < nl; ++i)
      {
        for (size_t j = 0; j < nr; ++j)
        {
          const size_t q = i * nr + j;
          EXPECT_EQ(product.getWeight(q), leftWeights[i] * rightWeights[j]);
          ASSERT_EQ(size_t(product.getPoint(q).size()), ldim + rdim);
          for (size_t k = 0; k < ldim; ++k)
            EXPECT_EQ(product.getPoint(q)(k), leftPoints[i](k));
          for (size_t k = 0; k < rdim; ++k)
            EXPECT_EQ(product.getPoint(q)(ldim + k), rightPoints[j](k));
        }
      }
    }
  }

  /// @brief Empty products do not query unavailable factor samples.
  TEST(TensorProduct, EmptyFactors)
  {
    for (const bool emptyLeft : {false, true})
    {
      const std::vector<Math::SpatialPoint> points{Math::SpatialPoint{Real(1) / 3}};
      const ObservedProductFactor left(
        emptyLeft ? std::vector<Math::SpatialPoint>{} : points,
        emptyLeft ? std::vector<Real>{} : std::vector<Real>{1});
      const ObservedProductFactor right(
        emptyLeft ? points : std::vector<Math::SpatialPoint>{},
        emptyLeft ? std::vector<Real>{1} : std::vector<Real>{});
      const QF::TensorProduct product(
        Geometry::Polytope::Type::Quadrilateral, left, right);
      EXPECT_EQ(product.getSize(), 0u);
      EXPECT_EQ(left.pointReads() + right.pointReads(), 0u);
      EXPECT_EQ(left.weightReads() + right.weightReads(), 0u);
    }
  }

  /// @brief Construction and copying retain zeroed inactive fixed capacity.
  TEST(TensorProduct, InactiveCapacityIsZero)
  {
    for (const auto& [leftDimension, rightDimension] : {std::pair{0u, 0u}, {0u, 1u},
           {1u, 0u}, {0u, 2u}, {1u, 1u}, {2u, 0u}, {1u, 2u}, {2u, 1u}})
    {
      SCOPED_TRACE(leftDimension);
      SCOPED_TRACE(rightDimension);
      const size_t dimension = leftDimension + rightDimension;
      Math::SpatialPoint leftPoint(leftDimension), rightPoint(rightDimension);
      for (size_t j = 0; j < leftDimension; ++j)
        leftPoint(j) = Real(j + 1) / 7;
      for (size_t j = 0; j < rightDimension; ++j)
        rightPoint(j) = Real(j + 2) / 11;
      const ObservedProductFactor left({leftPoint}, {1}), right({rightPoint}, {1});
      using G = Geometry::Polytope::Type;
      const auto geometry = dimension == 0 ? G::Point
        : dimension == 1                   ? G::Segment
        : dimension == 2                   ? G::Quadrilateral
                                           : G::Hexahedron;
      const QF::TensorProduct product(geometry, left, right);
      const QF::TensorProduct copy(product);
      for (const auto* rule : {&product, &copy})
      {
        ASSERT_EQ(rule->getSize(), 1u);
        const auto& point = rule->getPoint(0);
        ASSERT_EQ(point.size(), dimension);
        for (size_t j = dimension; j < Math::SpatialPoint::MaxSize; ++j)
          EXPECT_EQ(point.getData()(j), Real(0));
      }
    }
  }

  /// @brief Aliased factors preserve three-factor order and owned copies.
  TEST(TensorProduct, AliasedThreeFactorProductAndOwnedCopy)
  {
    const std::vector<Math::SpatialPoint> points{
      Math::SpatialPoint{Real(1) / 4}, Math::SpatialPoint{Real(3) / 4}};
    const std::vector<Real> weights{Real(1) / 3, Real(2) / 3};
    const auto product = [&] {
      const ObservedProductFactor factor(points, weights);
      return QF::TensorProduct(
        Geometry::Polytope::Type::Hexahedron, factor, factor, factor);
    }();
    const QF::TensorProduct copy(product);
    EXPECT_NE(copy.getCacheIdentity(), product.getCacheIdentity());
    ASSERT_EQ(product.getSize(), 8u);
    ASSERT_EQ(copy.getSize(), product.getSize());
    for (size_t i = 0; i < 2; ++i)
    {
      for (size_t j = 0; j < 2; ++j)
      {
        for (size_t k = 0; k < 2; ++k)
        {
          const size_t q = (i * 2 + j) * 2 + k;
          EXPECT_EQ(product.getWeight(q), (weights[i] * weights[j]) * weights[k]);
          EXPECT_EQ(copy.getWeight(q), product.getWeight(q));
          EXPECT_EQ(product.getPoint(q)(0), points[i](0));
          EXPECT_EQ(product.getPoint(q)(1), points[j](0));
          EXPECT_EQ(product.getPoint(q)(2), points[k](0));
          for (size_t direction = 0; direction < 3; ++direction)
            EXPECT_EQ(copy.getPoint(q)(direction), product.getPoint(q)(direction));
        }
      }
    }
  }
}
