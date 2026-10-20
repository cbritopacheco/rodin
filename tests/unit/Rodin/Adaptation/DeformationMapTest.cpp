/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Adaptation/DeformationMap.h"
#include "Rodin/Geometry.h"
#include "Rodin/Location/AABB.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/Variational.h"
#include "../../../convergence/Convergence.h"

using namespace Rodin;
using namespace Rodin::Adaptation;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    struct CountingLocator
    {
      Location::AABB<LocalMesh> locator;
      mutable std::size_t calls = 0;

      Optional<Point> locate(std::size_t dimension, const Math::SpatialPoint& x) const
      {
        ++calls;
        return locator.locate(dimension, x);
      }
    };

    template <std::size_t Dimension, std::size_t Order>
    void checkPointEvaluation()
    {
      auto mesh = [] {
        if constexpr (Dimension == 2)
          return LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
        else
          return LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
      }();
      for (std::size_t from = 0; from <= Dimension; ++from)
      {
        for (std::size_t to = 0; to <= Dimension; ++to)
        {
          if (from != to)
            mesh.getConnectivity().compute(from, to);
        }
      }
      auto space = [&] {
        if constexpr (Order == 1)
          return P1<Math::SpatialVector<Real>, LocalMesh>(mesh, Dimension);
        else
          return H1(std::integral_constant<std::size_t, Order>{}, mesh, Dimension);
      }();
      constexpr Real Expansion = Real(0.02);
      constexpr Real Tolerance = Real(1e-11);
      GridFunction displacement(space);
      displacement = VectorFunction(Dimension, [Expansion](const Point& point) {
        return Math::SpatialVector<Real>(Expansion * point.getPhysicalCoordinates());
      });
      CountingLocator locator{Location::AABB<LocalMesh>(mesh)};
      DeformationMap deformation(displacement, locator);
      auto face = mesh.getFace();
      ASSERT_TRUE(face);
      const auto& firstRule = QF::PolytopeQuadratureFormula::get(2, face->getGeometry());
      const auto& secondRule = QF::PolytopeQuadratureFormula::get(6, face->getGeometry());
      const auto& firstQuadrature = face->getQuadrature(firstRule);
      const auto& secondQuadrature = face->getQuadrature(secondRule);
      const Point& first = firstQuadrature.getPoint(0);
      const Point& second = secondQuadrature.getPoint(0);
      ASSERT_GT((first.getPhysicalCoordinates() - second.getPhysicalCoordinates()).norm(),
        Tolerance);
      const auto check = [&](const IntegrationPoint& ip, Real expansion) {
        const auto expected = Math::SpatialPoint(
          (Real(1) + expansion) * ip.getPoint().getPhysicalCoordinates());
        EXPECT_LT(
          (deformation.setIntegrationPoint(ip).getPoint().getPhysicalCoordinates() -
            expected)
            .norm(),
          Tolerance);
      };

      // Generic points share index zero, but neither their basis nor position.
      const IntegrationPoint genericFirst(first), genericSecond(second);
      EXPECT_LT((deformation.getDisplacementValue(first, genericFirst) -
                  Expansion * first.getPhysicalCoordinates())
                  .norm(),
        Tolerance);
      EXPECT_LT((deformation.getDisplacementValue(second, genericSecond) -
                  Expansion * second.getPhysicalCoordinates())
                  .norm(),
        Tolerance);
      check(genericFirst, Expansion);
      check(genericSecond, Expansion);
      check(genericFirst, Expansion);
      EXPECT_EQ(locator.calls, 3u);

      deformation.invalidate();
      const IntegrationPoint indexedFirst(firstQuadrature, 0);
      const IntegrationPoint indexedSecond(secondQuadrature, 0);
      deformation.setPolytope(*face).setIntegrationPoint(indexedFirst);
      const std::size_t beforeRepeat = locator.calls;
      const auto expectedFirst =
        Math::SpatialPoint((Real(1) + Expansion) * first.getPhysicalCoordinates());
      EXPECT_LT((deformation.getPoint().getPhysicalCoordinates() - expectedFirst).norm(),
        Tolerance);
      deformation.getPoint();
      EXPECT_EQ(locator.calls, beforeRepeat);
      const auto copy = deformation;
      EXPECT_THROW(copy.getPoint(), Alert::Exception);
      check(indexedSecond, Expansion);
      EXPECT_EQ(locator.calls, beforeRepeat + 1);

      displacement *= Real(2);
      deformation.invalidate();
      check(indexedSecond, Real(2) * Expansion);
      EXPECT_EQ(locator.calls, beforeRepeat + 2);
    }
  }

  TEST(Rodin_Adaptation_DeformationMap, P1TrianglePointAndQuadratureIdentity)
  {
    checkPointEvaluation<2, 1>();
  }

  TEST(Rodin_Adaptation_DeformationMap, P2TrianglePointAndQuadratureIdentity)
  {
    checkPointEvaluation<2, 2>();
  }

  TEST(Rodin_Adaptation_DeformationMap, P3TrianglePointAndQuadratureIdentity)
  {
    checkPointEvaluation<2, 3>();
  }

  TEST(Rodin_Adaptation_DeformationMap, P1TetrahedronPointAndQuadratureIdentity)
  {
    checkPointEvaluation<3, 1>();
  }

  TEST(Rodin_Adaptation_DeformationMap, P2TetrahedronPointAndQuadratureIdentity)
  {
    checkPointEvaluation<3, 2>();
  }

  TEST(Rodin_Adaptation_DeformationMap, P3TetrahedronPointAndQuadratureIdentity)
  {
    checkPointEvaluation<3, 3>();
  }

  /// @brief Failed location diagnoses unsupported source dimensions and clears binding state.
  TEST(Rodin_Adaptation_DeformationMap, FallbackRequiresCellOrFacet)
  {
    struct MissingLocator
    {
      Optional<Point> locate([[maybe_unused]] size_t dimension,
        [[maybe_unused]] const Math::SpatialPoint& physical) const
      {
        return {};
      }
    };
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = *mesh.getCell();
    const auto vertex = *mesh.getPolytope(0, 0);
    VectorFunction zero(
      size_t(2), [](const Point&) { return Math::SpatialVector<Real>{0, 0}; });
    MissingLocator locator;
    DeformationMap deformation(zero, locator);
    const Point source(cell, Math::SpatialPoint{Real(0.25), Real(0.25)});
    const IntegrationPoint cellPoint(source);
    deformation.setPolytope(cell).setIntegrationPoint(cellPoint);
    EXPECT_EQ(
      (deformation.getPoint().getPhysicalCoordinates() - source.getPhysicalCoordinates())
        .norm(),
      0);
    const Point vertexSource(vertex, Math::SpatialPoint(0));
    const IntegrationPoint vertexPoint(vertexSource);
    EXPECT_THROW(deformation.setIntegrationPoint(vertexPoint), Alert::Exception);
    EXPECT_THROW(deformation.getPoint(), Alert::Exception);
  }

  /**
   * @brief Located-point reuse follows reference association and local state.
   *
   * Zero displacement isolates the cache from displacement updates. Formula
   * assignment and geometry dilation preserve entity indices but must trigger
   * relocation. Unproven point/formula/index bindings always relocate, even
   * when their originating index is unchanged. A counting cell locator tests
   * this contract independently of bounding-volume hierarchy rebuilding.
   */
  TEST(Rodin_Adaptation_DeformationMap, ReferenceAssociationAndStateTransitions)
  {
    struct CellLocator
    {
      explicit CellLocator(const LocalMesh& mesh)
        : mesh(mesh),
          calls(0),
          fail(false)
      {}
      Optional<Point> locate(size_t dimension, const Math::SpatialPoint& physical) const
      {
        ++calls;
        if (fail)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Locator failure." << Alert::Raise;
        }
        const auto cell = mesh.getPolytope(dimension, 0);
        Math::SpatialPoint reference;
        cell->getTransformation().inverse(reference, physical);
        return Point(*cell, reference, physical);
      }
      const LocalMesh& mesh;
      mutable size_t calls;
      bool fail;
    };
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      auto mesh = Convergence::UniformGrid(geometry).makeMesh(2);
      const auto cell = *mesh.getCell();
      const size_t dimension = mesh.getDimension();
      VectorFunction zero(dimension, [dimension](const Point&) {
        Math::SpatialVector<Real> value(dimension);
        value.setZero();
        return value;
      });
      CellLocator locator(mesh);
      DeformationMap deformation(zero, locator);
      const auto& formula =
        QF::QuadratureFormula::get<QF::GaussLegendre>(geometry, size_t{2});
      const PolytopeQuadrature first(cell, formula);
      const IntegrationPoint original(first, 0);
      EXPECT_THROW(deformation.getPoint(), Alert::Exception);
      deformation.setPolytope(cell).setIntegrationPoint(original);
      deformation.getPoint();
      deformation.getPoint();
      EXPECT_EQ(locator.calls, 1);
      const auto& replacement =
        QF::QuadratureFormula::get<QF::GaussLegendre>(geometry, size_t{4});
      const PolytopeQuadrature second(cell, replacement);
      const IntegrationPoint reassigned(second, 0);
      const auto secondPhysical =
        deformation.setIntegrationPoint(reassigned).getPoint().getPhysicalCoordinates();
      EXPECT_EQ(locator.calls, 2);
      for (size_t j = 0; j < dimension; ++j)
        EXPECT_EQ(secondPhysical(j), reassigned.getPoint().getPhysicalCoordinates()(j));
      mesh.scale(2);
      const PolytopeQuadrature third(cell, replacement);
      const IntegrationPoint rescaled(third, 0);
      const auto thirdPhysical =
        deformation.setIntegrationPoint(rescaled).getPoint().getPhysicalCoordinates();
      EXPECT_EQ(locator.calls, 3);
      for (size_t j = 0; j < dimension; ++j)
        EXPECT_EQ(thirdPhysical(j), rescaled.getPoint().getPhysicalCoordinates()(j));
      for (const size_t sample : {size_t(0), size_t(1), size_t(0)})
      {
        const IntegrationPoint unproven(third.getPoint(sample), &formula, 0);
        const auto actual =
          deformation.setIntegrationPoint(unproven).getPoint().getPhysicalCoordinates();
        for (size_t j = 0; j < dimension; ++j)
          EXPECT_EQ(actual(j), unproven.getPoint().getPhysicalCoordinates()(j));
      }
      EXPECT_EQ(locator.calls, 6);
      locator.fail = true;
      EXPECT_THROW(deformation.setIntegrationPoint(rescaled), Alert::Exception);
      EXPECT_THROW(deformation.getPoint(), Alert::Exception);
      locator.fail = false;
      EXPECT_NO_THROW(deformation.setIntegrationPoint(rescaled));
      EXPECT_NO_THROW(deformation.getPoint());
    }
  }
}
