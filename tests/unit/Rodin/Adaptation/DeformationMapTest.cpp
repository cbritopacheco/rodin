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
#include "Rodin/Variational.h"

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
          (deformation.getMovedPoint(ip).getPhysicalCoordinates() - expected).norm(),
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
      const IntegrationPoint indexedFirst(first, &firstRule, 0);
      const IntegrationPoint indexedSecond(second, &secondRule, 0);
      check(indexedFirst, Expansion);
      const std::size_t beforeRepeat = locator.calls;
      check(indexedFirst, Expansion);
      EXPECT_EQ(locator.calls, beforeRepeat);
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
}
