/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Stable quadrature ownership through cache growth and empty-dimensional metrics. */
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include "Rodin/Geometry.h"
#include "Rodin/QF/Centroid.h"
#include "Rodin/QF/QuadratureFormula.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  /// @brief External PointBase implementation owning its reference sample.
  class OwnedReferencePoint final : public PointBase
  {
    public:
      OwnedReferencePoint(const Polytope& cell, const Math::SpatialPoint& reference)
        : PointBase(cell),
          m_reference(reference)
      {}
      OwnedReferencePoint(const Math::SpatialMatrix<Real>& jacobian, const Polytope& cell,
        const Math::SpatialPoint& reference)
        : PointBase(jacobian, cell),
          m_reference(reference)
      {}
      OwnedReferencePoint(const Math::SpatialMatrix<Real>& jacobian, const Polytope& cell,
        const QF::QuadratureFormula& formula)
        : PointBase(jacobian, cell, formula, 0),
          m_reference(formula.getPoint(0))
      {}
      const Math::SpatialPoint& getReferenceCoordinates() const override
      {
        return m_reference;
      }

    private:
      Math::SpatialPoint m_reference;
  };
}

/// @brief External point implementations own seeded data and invalidate all geometry on rebinding.
TEST(Geometry_CacheGrowth, PointBaseExtensionOwnsSamplesAndRebinds)
{
  constexpr Real GeometryTolerance = 1e-12;
  for (const auto geometry :
    {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    const Polytope::Traits traits(geometry);
    const size_t dimension = traits.getDimension(), vertices = traits.getVertexCount();
    for (size_t embedding = dimension; embedding <= 3; ++embedding)
    {
      SCOPED_TRACE(::testing::Message() << int(geometry) << ':' << embedding);
      LocalMesh::Builder builder;
      builder.initialize(embedding).nodes(2 * vertices);
      for (Index cell = 0; cell < 2; ++cell)
      {
        IndexArray connectivity(vertices);
        for (size_t local = 0; local < vertices; ++local)
        {
          Math::SpatialPoint physical(embedding);
          for (size_t row = 0; row < embedding; ++row)
          {
            physical[row] = Real(10 * cell + row + 1);
            for (size_t column = 0; column < dimension; ++column)
            {
              const Real entry = row == column ? Real(2 + cell) : Real(row >= dimension);
              physical[row] += entry * traits.getVertex(local)[column];
            }
          }
          builder.vertex(physical);
          connectivity[local] = cell * vertices + local;
        }
        if (dimension > 0)
          builder.polytope(geometry, connectivity);
      }
      auto mesh = builder.finalize();
      const auto first = *mesh.getPolytope(dimension, 0);
      const auto second = *mesh.getPolytope(dimension, 1);
      const auto reference = traits.getCentroid();
      OwnedReferencePoint lazy(first, reference);
      std::optional<OwnedReferencePoint> seeded, sampled;
      {
        Math::SpatialMatrix<Real> jacobian(embedding, dimension);
        for (size_t row = 0; row < embedding; ++row)
        {
          for (size_t column = 0; column < dimension; ++column)
            jacobian(row, column) = row == column ? Real(2) : Real(row >= dimension);
        }
        seeded.emplace(jacobian, first, reference);
        const auto& formula = QF::QuadratureFormula::get<QF::Centroid>(geometry);
        sampled.emplace(jacobian, first, formula);
        jacobian.setZero();
      }
      const auto check = [&](const PointBase& point, size_t cell) {
        const Real scale = Real(2 + cell);
        for (size_t row = 0; row < embedding; ++row)
        {
          Real expected = Real(10 * cell + row + 1);
          for (size_t column = 0; column < dimension; ++column)
          {
            const Real entry = row == column ? scale : Real(row >= dimension);
            expected += entry * reference[column];
            EXPECT_NEAR(point.getJacobian()(row, column), entry, GeometryTolerance);
          }
          EXPECT_NEAR(point.getPhysicalCoordinates()[row], expected, GeometryTolerance);
        }
        // det(s^2 I + e 11^T) = s^(2d) (1 + e d / s^2).
        const Real density = std::pow(scale, Real(dimension)) *
          std::sqrt(1 + Real((embedding - dimension) * dimension) / (scale * scale));
        EXPECT_NEAR(point.getJacobianDeterminant(), density, GeometryTolerance);
        EXPECT_NEAR(point.getDistortion(), density, GeometryTolerance);
        const auto inverse = point.getJacobianInverse();
        EXPECT_EQ(inverse.rows(), dimension);
        EXPECT_EQ(inverse.cols(), embedding);
        EXPECT_LT((inverse * point.getJacobian() -
                    Math::SpatialMatrix<Real>::Identity(dimension, dimension))
                    .norm(),
          GeometryTolerance);
      };
      for (auto* point : {&lazy, &*seeded, &*sampled})
      {
        check(*point, 0);
        OwnedReferencePoint copied(*point);
        OwnedReferencePoint moved(std::move(copied));
        check(moved, 0);
        moved.setPolytope(second);
        check(moved, 1);
      }
    }
  }
}

/// @brief Reserved slot and dimension growth preserve cached objects and their physical maps.
TEST(Geometry_CacheGrowth, QuadratureSlotsAndDimensionBuckets)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {5});
  const auto& centroid =
    QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Segment);
  const auto& gaussian = QF::QuadratureFormula::get(size_t{3}, Polytope::Type::Segment);
  const std::array<const QF::QuadratureFormula*, 2> formulas{&centroid, &gaussian};
  std::array<const PolytopeQuadrature*, 2> retained;
  PolytopeQuadratureIndex index;
  index.initialize(1);
  size_t factoryCalls = 0;
  auto get = [&](const auto& cache, Index cell, size_t formula) -> const auto& {
    return cache.get({1, cell}, mesh.getCellCount(), *formulas[formula], [&]() {
      ++factoryCalls;
      return std::make_unique<PolytopeQuadrature>(
        *mesh.getCell(cell), *formulas[formula]);
    });
  };
  for (size_t formula = 0; formula < formulas.size(); ++formula)
    retained[formula] = &get(index, 0, formula);
  // reserve-like growth must not discard slots, even when a smaller extent is requested later.
  index.resize(1, 2 * mesh.getCellCount());
  index.resize(1, mesh.getCellCount());
  index.initialize(3);
  PolytopeQuadratureIndex moved(std::move(index));
  constexpr Real CoordinateTolerance = 1e-12;
  for (size_t repetition = 0; repetition < 2; ++repetition)
  {
    for (size_t formula = 0; formula < formulas.size(); ++formula)
    {
      EXPECT_EQ(&get(moved, 0, formula), retained[formula]);
      for (Index cell = 0; cell < mesh.getCellCount(); ++cell)
      {
        const auto& quadrature = get(moved, cell, formula);
        EXPECT_EQ(&quadrature.getQuadratureFormula(), formulas[formula]);
        for (size_t qp = 0; qp < formulas[formula]->getSize(); ++qp)
        {
          const auto& point = quadrature.getPoint(qp);
          EXPECT_NEAR(point(0), Real(cell) + formulas[formula]->getPoint(qp)[0],
            CoordinateTolerance);
          EXPECT_EQ(point.getJacobianDeterminant(), 1);
        }
      }
    }
  }
  EXPECT_EQ(factoryCalls, formulas.size() * mesh.getCellCount());
  moved.clear();
  for (size_t formula = 0; formula < formulas.size(); ++formula)
    EXPECT_EQ(get(moved, 0, formula).getPoint(0).getPolytope().getIndex(), 0u);
  EXPECT_EQ(factoryCalls, formulas.size() * (mesh.getCellCount() + 1));
}

/// @brief A point in zero ambient dimensions has the empty identity map and counting measure one.
TEST(Geometry_CacheGrowth, ZeroAmbientPointMetrics)
{
  const Math::SpatialPoint empty(0);
  auto mesh = LocalMesh::Builder().initialize(0).nodes(1).vertex(empty).finalize();
  ASSERT_EQ(mesh.getSpaceDimension(), 0u);
  ASSERT_EQ(mesh.getDimension(), 0u);
  ASSERT_EQ(mesh.getVertexCount(), 1u);
  const auto cell = *mesh.getCell(0);
  const Point point(cell, empty);
  EXPECT_EQ(point.getPhysicalCoordinates().size(), 0u);
  EXPECT_EQ(point.getReferenceCoordinates().size(), 0u);
  EXPECT_EQ(point.getJacobian().rows(), 0u);
  EXPECT_EQ(point.getJacobian().cols(), 0u);
  EXPECT_EQ(point.getJacobianInverse().rows(), 0u);
  EXPECT_EQ(point.getJacobianInverse().cols(), 0u);
  EXPECT_EQ(point.getJacobianDeterminant(), 1);
  EXPECT_EQ(point.getDistortion(), 1);
  EXPECT_EQ(mesh.getMeasure(0), 1);
}
