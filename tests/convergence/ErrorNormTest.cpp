/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Norm integration lifetime and unchanged-summation regressions. */
#include <gtest/gtest.h>
#include "Convergence.h"
#include "CurvedGeometry.h"
#include "LiftedErrorNorm.h"
#include "SineMap.h"

namespace Rodin::Tests::Convergence
{
  namespace
  {
    /** Counts mesh-cache requests without changing the cache's semantics. */
    class ObservedMesh : public Geometry::LocalMesh
    {
      public:
        explicit ObservedMesh(Geometry::LocalMesh&& mesh)
          : Geometry::LocalMesh(std::move(mesh))
        {}

        const Geometry::PolytopeQuadrature& getQuadrature(
          size_t d, Index i, const QF::QuadratureFormulaBase& qf) const override
        {
          ++requests;
          return Geometry::LocalMesh::getQuadrature(d, i, qf);
        }

        mutable size_t requests = 0;
    };
  }

  TEST(ErrorNormTest, CellLocalQuadraturePreservesNormsWithoutMeshCacheRequests)
  {
    using namespace Variational;
    using Type = Geometry::Polytope::Type;
    for (const auto geometry : {Type::Segment, Type::Triangle, Type::Quadrilateral,
           Type::Tetrahedron, Type::Pyramid, Type::Hexahedron, Type::Wedge})
    {
      SCOPED_TRACE(UniformGrid::getGeometryName(geometry));
      ObservedMesh mesh(UniformGrid(geometry).makeMesh(2));
      const size_t dim = mesh.getDimension();
      P1<Real, Geometry::LocalMesh> scalarSpace(mesh);
      P1<Math::SpatialVector<Real>, Geometry::LocalMesh> vectorSpace(mesh, dim);
      const GridFunction scalar(scalarSpace);
      const GridFunction vector(vectorSpace);
      const auto exact = [](const Geometry::Point& p) {
        Real value = 1;
        for (size_t j = 0; j < p.getDimension(); ++j)
          value += p(j);
        return value;
      };
      const auto gradient = [dim](const Geometry::Point&) {
        Math::SpatialVector<Real> value(dim);
        for (size_t j = 0; j < dim; ++j)
          value(j) = 1;
        return value;
      };
      const auto exactVector = [](const Geometry::Point& p) {
        return p.getPhysicalCoordinates();
      };
      const auto jacobian = [dim](const Geometry::Point&) {
        Math::SpatialMatrix<Real> value(dim, dim);
        value.setIdentity();
        return value;
      };

      for (const size_t order : {size_t(4), size_t(16), size_t(18)})
      {
        SCOPED_TRACE(order);
        mesh.requests = 0;
        const auto value = ErrorNorm::computeL2(mesh, scalar, exact, order);
        const auto scalarError = ErrorNorm::compute(mesh, scalar, exact, gradient, order);
        const auto vectorError =
          ErrorNorm::computeVector(mesh, vector, exactVector, jacobian, order);
        const auto divergence = ErrorNorm::computeDivergenceL2(mesh, vector, order);
        EXPECT_EQ(mesh.requests, 0u);

        // Original mesh-cached path: identical points, weights and sum order.
        Real scalarSquared = 0, vectorSquared = 0, derivativeSquared = 0;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const auto& qf = QF::PolytopeQuadratureFormula::get(order, geometry);
          const auto& q = cell->getQuadrature(qf);
          for (size_t qp = 0; qp < q.getSize(); ++qp)
          {
            const auto& p = q.getPoint(qp);
            const Real weight = qf.getWeight(qp) * p.getDistortion();
            scalarSquared += weight * ErrorNorm::squaredMagnitude(exact(p));
            vectorSquared += weight * ErrorNorm::squaredMagnitude(exactVector(p));
            derivativeSquared += weight * Real(dim);
          }
        }
        EXPECT_EQ(value, std::sqrt(scalarSquared));
        EXPECT_EQ(scalarError.getL2(), value);
        EXPECT_EQ(scalarError.getH1Seminorm(), std::sqrt(derivativeSquared));
        EXPECT_EQ(vectorError.getL2(), std::sqrt(vectorSquared));
        EXPECT_EQ(vectorError.getH1Seminorm(), std::sqrt(derivativeSquared));
        EXPECT_EQ(divergence, 0);
        EXPECT_EQ(mesh.requests, mesh.getPolytopeCount(dim));
      }
    }
  }

  TEST(ErrorNormTest, LiftedNormsDoNotRetainEitherChartQuadrature)
  {
    using Type = Geometry::Polytope::Type;
    for (const auto geometry : {Type::Segment, Type::Triangle, Type::Quadrilateral,
           Type::Tetrahedron, Type::Pyramid, Type::Hexahedron, Type::Wedge})
    {
      SCOPED_TRACE(UniformGrid::getGeometryName(geometry));
      const auto check = [&]<class Range>() {
        using Scalar = typename FormLanguage::Traits<Range>::ScalarType;
        struct ConstantData
        {
            size_t dim;
            Range getSolution(const Math::SpatialPoint&) const
            {
              if constexpr (FormLanguage::IsVectorRange<Range>::Value)
              {
                Range value(dim);
                for (size_t j = 0; j < dim; ++j)
                  value(j) = phase();
                return value;
              }
              else
                return phase();
            }
            Scalar phase() const
            {
              if constexpr (std::is_same_v<Scalar, Complex>)
                return Scalar(1, 2);
              else
                return Scalar(1);
            }
            auto getGradient(const Math::SpatialPoint&) const
            {
              Math::SpatialVector<Scalar> value(dim);
              value.setZero();
              return value;
            }
            auto getJacobian(const Math::SpatialPoint&) const
            {
              Math::SpatialMatrix<Scalar> value(dim, dim);
              value.setZero();
              return value;
            }
        };

        ObservedMesh reference(UniformGrid(geometry).makeMesh(2));
        ObservedMesh represented{Geometry::LocalMesh(reference)};
        CurvedGeometry mapping(represented, CurvedGeometry<ObservedMesh>::Map::Sine);
        mapping.template install<2>();
        const size_t dim = reference.getDimension();
        auto space = [&] {
          if constexpr (FormLanguage::IsVectorRange<Range>::Value)
            return Variational::P1<Range, Geometry::LocalMesh>(represented, dim);
          else
            return Variational::P1<Range, Geometry::LocalMesh>(represented);
        }();
        const Variational::GridFunction state(space);
        const ConstantData data{dim};
        const SineMap exactMap;
        for (const size_t order : {size_t(8), size_t(18)})
        {
          SCOPED_TRACE(order);
          // Independent cached reference path for the constant field defect.
          // Its gradient and geometry defects vanish on every regular chart.
          Real expectedSquared = 0;
          for (auto cell = reference.getCell(); cell; ++cell)
          {
            const auto& qf = QF::PolytopeQuadratureFormula::get(order, geometry);
            const auto& quadrature = cell->getQuadrature(qf);
            for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
            {
              const auto& p = quadrature.getPoint(qp);
              const Real weight = qf.getWeight(qp) * p.getDistortion() *
                exactMap.getJacobian(p.getPhysicalCoordinates()).determinant();
              expectedSquared += weight *
                ErrorNorm::squaredMagnitude(data.getSolution(p.getPhysicalCoordinates()));
            }
          }
          reference.requests = 0;
          const auto result = LiftedErrorNorm::compute(
            reference, represented, state, data, exactMap, order);
          EXPECT_EQ(reference.requests, 0u);
          EXPECT_EQ(represented.requests, 0u);
          EXPECT_EQ(result.field.getL2(), std::sqrt(expectedSquared));
          EXPECT_EQ(result.total.getL2(), result.field.getL2());
          EXPECT_EQ(result.geometry.getL2(), 0);
          for (const auto& error : {result.field, result.geometry, result.total})
            EXPECT_EQ(error.getH1Seminorm(), 0);
        }
      };
      check.operator()<Real>();
      check.operator()<Complex>();
      check.operator()<Math::SpatialVector<Real>>();
      check.operator()<Math::SpatialVector<Complex>>();
    }
  }
}
