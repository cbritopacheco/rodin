/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/// @file @brief Independent directional and linearity checks for the Holzapfel-Ogden tangent.
#include <gtest/gtest.h>
#include "Rodin/Solid/Constitutive/HolzapfelOgden.h"

using namespace Rodin;
using namespace Rodin::Solid;

namespace
{
  struct MaterialResponse
  {
    Real energy;
    Math::SpatialMatrix<Real> stress, tangent;
  };

  MaterialResponse evaluate(const HolzapfelOgden& law,
    const Math::SpatialMatrix<Real>& gradient, const Math::SpatialMatrix<Real>& increment)
  {
    KinematicState state(gradient.rows());
    state.setDisplacementGradient(gradient);
    ConstitutivePoint point(state);
    Math::SpatialVector<Real> fiber(gradient.rows());
    for (size_t i = 0; i < fiber.size(); ++i)
      fiber[i] = 1 + i;
    point.set<Tags::FiberDirection>(fiber);
    auto prepared = law;
    prepared.setConstitutivePoint(point);
    MaterialResponse response;
    response.energy = prepared.getStrainEnergyDensity();
    prepared.getFirstPiolaKirchhoffStress(response.stress);
    prepared.getMaterialTangent(response.tangent, increment);
    return response;
  }

  Real inner(const Math::SpatialMatrix<Real>& lhs, const Math::SpatialMatrix<Real>& rhs)
  {
    Real result = 0;
    for (size_t i = 0; i < lhs.rows(); ++i)
    {
      for (size_t j = 0; j < lhs.cols(); ++j)
        result += lhs(i, j) * rhs(i, j);
    }
    return result;
  }
}

TEST(HolzapfelOgden, AnalyticTangentAndDirectionScaling)
{
  const HolzapfelOgden law(2, 0.7, 0.4, 0.3, 0.6, 0.2, 11);
  for (const size_t dimension : {2u, 3u})
  {
    Math::SpatialMatrix<Real> gradient(dimension, dimension),
      direction(dimension, dimension);
    for (size_t i = 0; i < dimension; ++i)
    {
      for (size_t j = 0; j < dimension; ++j)
      {
        gradient(i, j) = i == j ? 0.1 * (i + 1) : 0.02 * (Real(i) - Real(j));
        direction(i, j) = 0.13 * (i + 1) - 0.07 * j;
      }
    }
    const auto response = evaluate(law, gradient, direction);
    const Real step = 1e-6;
    const auto plus = evaluate(law, gradient + step * direction, direction);
    const auto minus = evaluate(law, gradient - step * direction, direction);
    const Math::SpatialMatrix<Real> difference =
      (Real(0.5) / step) * (plus.stress - minus.stress);
    EXPECT_LT(
      (difference - response.tangent).norm() / std::max(Real(1), response.tangent.norm()),
      2e-8);
    EXPECT_NEAR(
      (plus.energy - minus.energy) / (2 * step), inner(response.stress, direction), 2e-8);
    for (const Real scale :
      {Real(0), Real(1e-10), Real(1e-4), Real(1), Real(1e4), Real(1e8)})
    {
      const auto actual = evaluate(law, gradient, scale * direction);
      const Math::SpatialMatrix<Real> expected = scale * response.tangent;
      EXPECT_LT(
        (actual.tangent - expected).norm() / std::max(Real(1), expected.norm()), 5e-13)
        << "scale=" << scale;
    }
    const Math::SpatialMatrix<Real> second = direction.transpose();
    const auto other = evaluate(law, gradient, second);
    EXPECT_NEAR(inner(second, response.tangent), inner(direction, other.tangent), 1e-12);
  }
}
