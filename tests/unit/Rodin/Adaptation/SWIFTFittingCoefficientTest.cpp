/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Reentrant fitting coefficients retain independent deformed points. */
#include <gtest/gtest.h>
#include <functional>
#include "Rodin/Adaptation/SWIFT/FittingForce.h"
#include "Rodin/Adaptation/SWIFT/FittingTensor.h"
#include "Rodin/Geometry.h"
#include "Rodin/Location/AABB.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    /**
     * @brief Checks nested evaluation on the same coefficient against analytic values.
     * @tparam Tensor Whether the fitting tensor rather than the force is evaluated.
     * @tparam Order Displacement finite element order.
     */
    template <bool Tensor, size_t Order>
    void checkReentrantEvaluation()
    {
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 2);
      mesh.getConnectivity().compute(1, 0);
      auto space = [&] {
        if constexpr (Order == 1)
          return P1<Math::SpatialVector<Real>>(mesh, 2);
        else
        {
          return H1<Order, Math::SpatialVector<Real>>(
            std::integral_constant<size_t, Order>{}, mesh, 2);
        }
      }();
      GridFunction displacement(space);
      constexpr Real Expansion = Real(0.02);
      displacement = VectorFunction(2, [Expansion](const Point& point) {
        return Math::SpatialVector<Real>(Expansion * point.getPhysicalCoordinates());
      });
      const Location::AABB<LocalMesh> locator(mesh);
      const auto cell = mesh.getCell();
      const Point first(*cell, Math::SpatialPoint{Real(0.15), Real(0.15)});
      const Point second(*cell, Math::SpatialPoint{Real(0.55), Real(0.15)});
      const IntegrationPoint outer(first), inner(second);
      constexpr Real Normalization = Real(0.7);
      constexpr Real Tolerance = Real(1e-12);
      const Adaptation::SWIFT::Loss loss(2);
      Adaptation::SWIFT::Parameters parameters;
      parameters.model.fit = Real(1.3);
      std::function<void()> nested;
      bool evaluatingNested = false;
      const auto reenter = [&] {
        if (nested && !evaluatingNested)
        {
          evaluatingNested = true;
          nested();
          evaluatingNested = false;
        }
      };
      const auto phi = RealFunction([&](const Point& point) {
        if constexpr (!Tensor)
          reenter();
        return Real(1) + point.getPhysicalCoordinates().squaredNorm();
      });
      const auto gradient = VectorFunction(2, [&](const Point& point) {
        if constexpr (Tensor)
          reenter();
        return Math::SpatialVector<Real>(Real(2) * point.getPhysicalCoordinates());
      });
      const auto coefficient = [&] {
        if constexpr (Tensor)
        {
          return Adaptation::SWIFT::FittingTensor(
            gradient, displacement, locator, parameters, Normalization, size_t(2));
        }
        else
        {
          return Adaptation::SWIFT::FittingForce(
            phi, gradient, displacement, locator, loss, Normalization, size_t(2));
        }
      }();
      const auto expected = [&](const Point& point) {
        const Math::SpatialPoint physical =
          (Real(1) + Expansion) * point.getPhysicalCoordinates();
        const Math::SpatialVector<Real> g = Real(2) * physical;
        if constexpr (Tensor)
        {
          Math::SpatialMatrix<Real> value(2, 2);
          for (size_t row = 0; row < 2; ++row)
          {
            for (size_t col = 0; col < 2; ++col)
              value(row, col) = parameters.model.fit * Normalization * g(row) * g(col);
          }
          return value;
        }
        else
        {
          const Real residual = Real(1) + physical.squaredNorm();
          return Math::SpatialVector<Real>(
            -Normalization * loss.getWeight(residual) * residual * g);
        }
      };
      nested = [&] {
        const auto actual = coefficient.getValue(inner);
        EXPECT_LT((actual - expected(second)).norm(), Tolerance);
      };
      EXPECT_LT((coefficient.getValue(outer) - expected(first)).norm(), Tolerance);
      // A second pass also observes changed borrowed displacement coefficients.
      displacement *= Real(0);
      nested = {};
      const auto actual = coefficient.getValue(outer);
      const auto directGradient =
        Math::SpatialVector<Real>(Real(2) * first.getPhysicalCoordinates());
      if constexpr (Tensor)
      {
        const Math::SpatialMatrix<Real> direct = parameters.model.fit * Normalization *
          (directGradient.getData().head(2) *
            directGradient.getData().head(2).transpose());
        EXPECT_LT((actual - direct).norm(), Tolerance);
      }
      else
      {
        const Real residual = Real(1) + first.getPhysicalCoordinates().squaredNorm();
        const Math::SpatialVector<Real> direct =
          -Normalization * loss.getWeight(residual) * residual * directGradient;
        EXPECT_LT((actual - direct).norm(), Tolerance);
      }
    }
  }

  /// @brief Nested scalar target evaluation cannot replace the outer force's point.
  TEST(Rodin_Adaptation_SWIFTFittingCoefficient, ReentrantForceP1P2)
  {
    checkReentrantEvaluation<false, 1>();
    checkReentrantEvaluation<false, 2>();
  }

  /// @brief Nested gradient evaluation cannot replace the outer tensor's point.
  TEST(Rodin_Adaptation_SWIFTFittingCoefficient, ReentrantTensorP1P2)
  {
    checkReentrantEvaluation<true, 1>();
    checkReentrantEvaluation<true, 2>();
  }
}
