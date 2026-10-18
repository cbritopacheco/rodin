/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <algorithm>
#include <cstring>
#include <vector>
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Geometry.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Variational/P1.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Observes metadata reads and logical sample visits without mutation.
   *
   * ## Architecture
   *
   * An owned dispatcher supplies immutable coordinates and weights. An optional
   * empty view supplies no samples. Observations change neither payload nor
   * logical identity. Copy construction owns independent formula contents and
   * resets observations; no borrowed fixture or address is used as an identity key.
   */
  class ConstructionFormula final : public QF::QuadratureFormulaBase
  {
    public:
      /**
       * @brief Constructs an observed immutable formula.
       * @param geometry Reference-cell family.
       * @param order Requested polynomial degree.
       * @param empty Whether the observed formula has no samples.
       */
      ConstructionFormula(Geometry::Polytope::Type geometry, size_t order, bool empty)
        : m_formula(order, geometry),
          m_empty(empty),
          m_sizeQueries(0),
          m_weightQueries(0),
          m_sampleVisits()
      {}

      /**
       * @brief Copies the payload under a fresh identity with no observations.
       * @param other Formula whose immutable payload is copied.
       */
      ConstructionFormula(const ConstructionFormula& other)
        : QuadratureFormulaBase(other),
          m_formula(other.m_formula),
          m_empty(other.m_empty),
          m_sizeQueries(0),
          m_weightQueries(0),
          m_sampleVisits()
      {}

      /**
       * @brief Counts a query of the immutable sample count.
       * @returns Number of samples, zero for the empty view.
       */
      size_t getSize() const override
      {
        ++m_sizeQueries;
        return m_empty ? 0 : m_formula.getSize();
      }

      /**
       * @brief Observes a weight query, which geometric construction does not need.
       * @param i Logical sample index.
       * @returns Stored reference weight.
       */
      Real getWeight(size_t i) const override
      {
        if (m_empty)
          Alert::MemberFunctionException(*this, __func__)
            << "An empty construction formula has no weights" << Alert::Raise;
        ++m_weightQueries;
        return m_formula.getWeight(i);
      }

      /**
       * @brief Records a logical sample read without changing its coordinates.
       * @param i Logical sample index.
       * @returns Stored reference coordinates.
       */
      const Math::SpatialPoint& getPoint(size_t i) const override
      {
        if (m_empty)
          Alert::MemberFunctionException(*this, __func__)
            << "An empty construction formula has no points" << Alert::Raise;
        m_sampleVisits.push_back(i);
        return m_formula.getPoint(i);
      }

      /**
       * @brief Clears observations without changing the formula's logical version.
       */
      void resetObservations() const
      {
        m_sizeQueries = 0;
        m_weightQueries = 0;
        m_sampleVisits.clear();
      }

      /**
       * @brief Borrows the unobserved payload used for exact input comparisons.
       * @returns The owned dispatcher; the empty flag is handled by the fixture.
       */
      const QF::PolytopeQuadratureFormula& getPayload() const
      {
        return m_formula;
      }

      /**
       * @brief Returns the observed metadata-query count.
       * @returns Number of calls to getSize().
       */
      size_t getSizeQueries() const
      {
        return m_sizeQueries;
      }

      /**
       * @brief Returns the observed weight-query count.
       * @returns Number of calls to getWeight().
       */
      size_t getWeightQueries() const
      {
        return m_weightQueries;
      }

      /**
       * @brief Borrows sample visits in their original call order.
       * @returns Logical indices passed to getPoint().
       */
      const std::vector<size_t>& getSampleVisits() const
      {
        return m_sampleVisits;
      }

      /**
       * @brief Copies the observed formula polymorphically.
       * @returns A newly owned copy with fresh observations and logical identity.
       */
      ConstructionFormula* copy() const noexcept override
      {
        return new ConstructionFormula(*this);
      }

    private:
      QF::PolytopeQuadratureFormula m_formula;
      bool m_empty;
      mutable size_t m_sizeQueries;
      mutable size_t m_weightQueries;
      mutable std::vector<size_t> m_sampleVisits;
  };

  /**
   * @brief Checks construction metadata and sample order on every native P1 chart.
   */
  class PolytopeQuadratureConstructionTest
    : public ::testing::TestWithParam<Geometry::Polytope::Type>
  {
    protected:
      /**
       * @brief Checks bounded metadata work, sample order and copied Jacobians.
       *
       * A cold reference fill may read size once; construction adds at most one
       * further non-diagnostic query. Assertion-enabled builds additionally
       * read size once for the container's consistency check and once for
       * each point's index check. Those diagnostic reads are counted explicitly;
       * the bound on production metadata reads is unchanged.
       * Perturbed vertex controls exercise sample-dependent
       * native tensor, wedge and pyramid Jacobians. Retained matrices must
       * reproduce the complete bulk output, including inactive storage bits.
       * @param warm Whether the native element's reference table is prefilled.
       * @param empty Whether the formula deliberately supplies no samples.
       */
      void verifyConstruction(bool warm, bool empty)
      {
        using namespace Geometry;
        const Polytope::Traits traits(GetParam());
        const size_t dimension = traits.getDimension();
        for (size_t embedding = std::max(size_t(1), dimension); embedding <= 3;
             ++embedding)
        {
          Mesh<Context::Local>::Builder builder;
          builder.initialize(embedding).nodes(traits.getVertexCount());
          IndexArray vertices(traits.getVertexCount());
          for (size_t a = 0; a < traits.getVertexCount(); ++a)
          {
            Math::SpatialPoint point(embedding);
            for (size_t j = 0; j < embedding; ++j)
            {
              point(j) = j < dimension ? traits.getVertex(a)(j) : Real(j + 1) / 4;
              point(j) += Real((a + 2) * (j + 1) % 5) / 64;
            }
            builder.vertex(point);
            vertices(a) = a;
          }
          if (dimension > 0)
            builder.polytope(GetParam(), vertices);
          auto mesh = builder.finalize();
          const auto cell = *mesh.getPolytope(dimension, 0);
          const Variational::RealP1Element element(GetParam());
          for (const size_t order : {0, 2, 16})
          {
            SCOPED_TRACE(::testing::Message()
              << "embedding=" << embedding << ", order=" << order << ", warm=" << warm
              << ", empty=" << empty);
            ConstructionFormula formula(GetParam(), order, empty);
            const size_t expectedSize = empty ? 0 : formula.getPayload().getSize();
            const size_t identity = formula.getCacheIdentity();
            if (warm)
              element.getTabulation(formula);
            formula.resetObservations();
            const PolytopeQuadrature quadrature(cell, formula);
            size_t sizeQueryBudget = size_t(1) + !warm;
#ifndef NDEBUG
            sizeQueryBudget += expectedSize + 1;
#endif
            EXPECT_LE(formula.getSizeQueries(), sizeQueryBudget);
            EXPECT_EQ(formula.getWeightQueries(), 0u);
            EXPECT_EQ(formula.getCacheIdentity(), identity);
            EXPECT_EQ(quadrature.getQuadratureFormulaIdentity(), identity);
            ASSERT_EQ(quadrature.getSize(), expectedSize);
            const auto& visits = formula.getSampleVisits();
            ASSERT_EQ(visits.size(), expectedSize * (warm ? 1 : 2));
            for (size_t read = 0; read < visits.size(); ++read)
              EXPECT_EQ(visits[read], read % expectedSize);
            std::vector<Math::SpatialMatrix<Real>> jacobians;
            cell.getTransformation().jacobian(jacobians, formula);
            ASSERT_EQ(jacobians.size(), expectedSize);
            for (size_t qp = 0; qp < expectedSize; ++qp)
            {
              const auto& actual = quadrature.getPoint(qp).getReferenceCoordinates();
              const auto& expected = formula.getPayload().getPoint(qp);
              ASSERT_EQ(actual.size(), dimension);
              ASSERT_EQ(actual.size(), expected.size());
              for (size_t j = 0; j < dimension; ++j)
                EXPECT_EQ(actual(j), expected(j));
              const auto& matrix = quadrature.getPoint(qp).getJacobian();
              ASSERT_EQ(matrix.rows(), embedding);
              ASSERT_EQ(matrix.cols(), dimension);
              EXPECT_EQ(
                std::memcmp(matrix.getData().data(), jacobians[qp].getData().data(),
                  matrix.getData().size() * sizeof(Real)),
                0);
            }
          }
        }
      }
  };

  /// @brief Cold construction needs at most one size read per owning operation.
  TEST_P(PolytopeQuadratureConstructionTest, ColdMetadataReadsAreConstant)
  {
    verifyConstruction(false, false);
  }

  /// @brief Warm construction must not repeat metadata reads for every sample.
  TEST_P(PolytopeQuadratureConstructionTest, WarmMetadataReadsAreConstant)
  {
    verifyConstruction(true, false);
  }

  /// @brief Empty formulas retain their identity without reading any sample.
  TEST_P(PolytopeQuadratureConstructionTest, EmptyFormulaCopiesNoSamples)
  {
    verifyConstruction(false, true);
    verifyConstruction(true, true);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PolytopeQuadratureConstructionTest,
    ::testing::Values(Geometry::Polytope::Type::Point, Geometry::Polytope::Type::Segment,
      Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Quadrilateral,
      Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Pyramid,
      Geometry::Polytope::Type::Hexahedron, Geometry::Polytope::Type::Wedge));
}
