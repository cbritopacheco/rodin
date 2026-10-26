/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Affine chart ownership and tabulated or direct element contracts. */
#include <gtest/gtest.h>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>
#include "Rodin/Geometry.h"
#include "Rodin/QF/QuadratureFormula.h"
#include "Rodin/Variational/P1.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  /// @brief A scalar element exposing the documented pointwise interface only.
  struct DirectElement
  {
    using RangeType = Real;
    Variational::RealP1Element element;
    auto getGeometry() const
    {
      return element.getGeometry();
    }
    size_t getCount() const
    {
      return element.getCount();
    }
    size_t getOrder() const
    {
      return element.getOrder();
    }
    const auto& getBasis(size_t local) const
    {
      return element.getBasis(local);
    }
  };

  class ChartOwnership : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      template <class Element>
      void check(Element element)
      {
        const auto geometry = element.getGeometry();
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        using Chart = ParametricTransformation<Element>;
        const Variational::RealP1Element nodes(geometry);
        for (size_t embedding = std::max(size_t(1), dimension); embedding <= 3;
             ++embedding)
        {
          SCOPED_TRACE(embedding);
          Math::SpatialMatrix<Real> linear(embedding, dimension);
          Math::SpatialPoint offset(embedding);
          for (size_t j = 0; j < embedding; ++j)
          {
            offset[j] = Real(j + 1);
            for (size_t i = 0; i < dimension; ++i)
              linear(j, i) = j == i ? Real(i + 2) : Real(j + 1);
          }
          PointCloud controls(embedding, element.getCount());
          for (size_t local = 0; local < element.getCount(); ++local)
          {
            Math::SpatialPoint reference(dimension);
            for (size_t i = 0; i < dimension; ++i)
              reference[i] = nodes.getNode(local)[i];
            const Math::SpatialPoint physical = offset + linear * reference;
            for (size_t j = 0; j < embedding; ++j)
              controls(j, local) = physical[j];
          }
          const auto& formula = QF::QuadratureFormula::get(size_t{4}, geometry);
          constexpr Real CoordinateTolerance = 1e-11;
          auto inspect = [&](const PolytopeTransformation& chart) {
            EXPECT_EQ(chart.getOrder(), element.getOrder());
            EXPECT_EQ(chart.getReferenceDimension(), dimension);
            EXPECT_EQ(chart.getPhysicalDimension(), embedding);
            if constexpr (std::is_same_v<Element, DirectElement>)
            {
              EXPECT_EQ(chart.getFactorOrder(), element.getOrder());
            }
            std::vector<Math::SpatialMatrix<Real>> jacobians;
            chart.jacobian(jacobians, formula);
            ASSERT_EQ(jacobians.size(), formula.getSize());
            for (size_t qp = 0; qp < formula.getSize(); ++qp)
            {
              Math::SpatialPoint physical;
              const auto& reference = formula.getPoint(qp);
              chart.transform(physical, reference);
              const Math::SpatialPoint expected = offset + linear * reference;
              EXPECT_LT((physical - expected).norm(), CoordinateTolerance);
              EXPECT_LT((jacobians[qp] - linear).norm(), CoordinateTolerance);
            }
          };
          std::unique_ptr<Chart> survivor;
          {
            Chart original(controls, element);
            Chart copied(original);
            Chart moved(std::move(copied));
            Chart movedControls{PointCloud(controls), element};
            Chart movedElement{controls, Element(element)};
            Chart movedBoth{PointCloud(controls), Element(element)};
            inspect(original);
            inspect(moved);
            inspect(movedControls);
            inspect(movedElement);
            inspect(movedBoth);
            survivor.reset(original.copy());
          }
          inspect(*survivor);
          EXPECT_EQ(survivor->getPointMatrix().rows(), embedding);
          EXPECT_EQ(survivor->getPointMatrix().cols(), element.getCount());
          for (size_t local = 0; local < element.getCount(); ++local)
          {
            for (size_t j = 0; j < embedding; ++j)
              EXPECT_EQ(survivor->getPointMatrix()(j, local), controls(j, local));
          }

          LocalMesh::Builder builder;
          builder.initialize(embedding).nodes(nodes.getCount());
          IndexArray vertices(nodes.getCount());
          for (size_t local = 0; local < nodes.getCount(); ++local)
          {
            builder.vertex(controls[local]);
            vertices[local] = local;
          }
          if (dimension > 0)
            builder.polytope(geometry, vertices);
          auto mesh = builder.finalize();
          mesh.setPolytopeTransformation({dimension, 0}, survivor.release());
          const auto cell = *mesh.getPolytope(dimension, 0);
          const PolytopeQuadrature quadrature(cell, formula);
          for (size_t qp = 0; qp < formula.getSize(); ++qp)
          {
            const auto& point = quadrature.getPoint(qp);
            const Math::SpatialPoint expected = offset + linear * formula.getPoint(qp);
            EXPECT_LT(
              (point.getPhysicalCoordinates() - expected).norm(), CoordinateTolerance);
            const Point seeded(linear, cell, formula.getPoint(qp));
            EXPECT_EQ((seeded.getJacobian() - linear).norm(), 0);
            EXPECT_LT(
              (seeded.getPhysicalCoordinates() - expected).norm(), CoordinateTolerance);
          }
        }
      }
  };

  /// @brief Tabulated charts own coefficients through all constructors, copies and moves.
  TEST_P(ChartOwnership, TabulatedElement)
  {
    check(Variational::RealP1Element(GetParam()));
  }

  /// @brief Elements without table or factor-order APIs retain the pointwise contract.
  TEST_P(ChartOwnership, DirectElement)
  {
    check(DirectElement{Variational::RealP1Element(GetParam())});
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, ChartOwnership,
    ::testing::Values(Polytope::Type::Point, Polytope::Type::Segment,
      Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge));
}
