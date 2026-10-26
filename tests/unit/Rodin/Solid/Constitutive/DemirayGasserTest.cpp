/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Independent energy, stress, tangent, and objectivity tests.
#include <gtest/gtest.h>
#include <Eigen/Geometry>
#include "Rodin/Solid/Constitutive/DemirayGasser.h"

using namespace Rodin;
using namespace Rodin::Solid;

namespace
{
  struct Response
  {
    Real energy;
    Math::SpatialMatrix<Real> stress, tangent;
  };

  Response evaluate(const DemirayGasser& material, const Eigen::Matrix3d& f,
    const Eigen::Matrix3d& direction)
  {
    KinematicState state(3);
    state.setDisplacementGradient(
      Math::SpatialMatrix<Real>(f - Eigen::Matrix3d::Identity()));
    ConstitutivePoint point(state);
    Math::SpatialVector<Real> fiber(3);
    fiber(0) = 1;
    fiber(1) = 0;
    fiber(2) = 0;
    point.set<Tags::FiberDirection>(fiber);
    auto prepared = material;
    prepared.setConstitutivePoint(point);
    Response result;
    result.energy = prepared.getStrainEnergyDensity();
    prepared.getFirstPiolaKirchhoffStress(result.stress);
    prepared.getMaterialTangent(result.tangent, Math::SpatialMatrix<Real>(direction));
    return result;
  }
}

TEST(Rodin_Solid_DemirayGasser, IndependentEnergyAndAllStressDerivatives)
{
  const DemirayGasser::Parameters parameters{3, 2, 4, 5, 17};
  const DemirayGasser material(parameters);
  Eigen::Matrix3d f;
  f << 1.12, .08, -.02, .03, .94, .04, -.01, .02, 1.05;
  Eigen::Matrix3d direction;
  direction << .12, -.03, .08, .02, -.09, .04, -.05, .07, .03;
  const auto actual = evaluate(material, f, direction);
  const Real scale = std::pow(f.determinant(), -2.0 / 3);
  const Real first = scale * f.squaredNorm();
  const Real fourth = scale * f.col(0).squaredNorm();
  const Real expected = .75 * std::expm1(2 * (first - 3)) +
    .4 * std::expm1(5 * (fourth - 1) * (fourth - 1)) +
    17 * std::pow(std::log(f.determinant()), 2);
  EXPECT_NEAR(actual.energy, expected, 1e-12);
  constexpr Real step = 1e-6;
  for (Index i = 0; i < 3; ++i)
  {
    for (Index j = 0; j < 3; ++j)
    {
      Eigen::Matrix3d plus = f, minus = f;
      plus(i, j) += step;
      minus(i, j) -= step;
      EXPECT_NEAR(actual.stress(i, j),
        (evaluate(material, plus, direction).energy -
          evaluate(material, minus, direction).energy) /
          (2 * step),
        2e-7);
    }
  }
  const auto plus = evaluate(material, f + step * direction, direction);
  const auto minus = evaluate(material, f - step * direction, direction);
  EXPECT_LT((actual.tangent - (plus.stress - minus.stress) / (2 * step)).norm(), 2e-7);
}

TEST(Rodin_Solid_DemirayGasser, ReferenceStateVolumeAndZeroExponentLimits)
{
  DemirayGasser::Parameters parameters{3, 0, 4, 0, 17};
  const DemirayGasser material(parameters);
  const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
  const auto rest = evaluate(material, identity, identity);
  EXPECT_DOUBLE_EQ(rest.energy, 0);
  EXPECT_LT(rest.stress.norm(), 1e-13);
  constexpr Real stretch = 1.1;
  const auto volume = evaluate(material, stretch * identity, identity);
  const Real logJ = 3 * std::log(stretch);
  EXPECT_NEAR(volume.energy, 17 * logJ * logJ, 1e-12);
  EXPECT_LT(
    (volume.stress - (34 * logJ / stretch) * Math::SpatialMatrix<Real>(identity)).norm(),
    1e-12);
  Eigen::Matrix3d f = identity;
  f(0, 1) = .2;
  const Real first = f.squaredNorm(), fourth = f.col(0).squaredNorm();
  EXPECT_NEAR(evaluate(material, f, identity).energy,
    1.5 * (first - 3) + 2 * (fourth - 1) * (fourth - 1), 1e-12);
}

TEST(Rodin_Solid_DemirayGasser, RotationCovarianceAndTangentLinearity)
{
  const DemirayGasser material(DemirayGasser::Parameters{});
  Eigen::Matrix3d f = Eigen::Matrix3d::Identity();
  f(0, 0) = 1.08;
  f(1, 2) = .04;
  const Eigen::Matrix3d q =
    Eigen::AngleAxisd(.6, Eigen::Vector3d(1, 2, 3).normalized()).toRotationMatrix();
  const Eigen::Matrix3d h = .1 * Eigen::Matrix3d::Ones();
  const auto a = evaluate(material, f, h), b = evaluate(material, q * f, q * h);
  EXPECT_NEAR(a.energy, b.energy, 1e-11);
  EXPECT_LT((b.stress - Math::SpatialMatrix<Real>(q) * a.stress).norm(), 1e-10);
  EXPECT_LT((b.tangent - Math::SpatialMatrix<Real>(q) * a.tangent).norm(), 1e-9);
  EXPECT_LT((evaluate(material, f, 2 * h).tangent - 2 * a.tangent).norm(), 1e-12);
}

TEST(Rodin_Solid_DemirayGasser, InvalidCoefficientsAndUnboundState)
{
  auto parameters = DemirayGasser::Parameters{};
  parameters.matrixScale = -1;
  EXPECT_THROW(DemirayGasser material(parameters), std::exception);
  const DemirayGasser material(DemirayGasser::Parameters{});
  EXPECT_THROW(material.getStrainEnergyDensity(), std::exception);
  KinematicState two(2);
  ConstitutivePoint point(two);
  auto other = material;
  EXPECT_THROW(other.setConstitutivePoint(point), std::exception);
}

TEST(Rodin_Solid_DemirayGasser, OwnedSnapshotAndChangedFiberBinding)
{
  DemirayGasser material(DemirayGasser::Parameters{});
  KinematicState state(3);
  Math::SpatialMatrix<Real> h(3, 3);
  h.setZero();
  h(0, 0) = .1;
  state.setDisplacementGradient(h);
  ConstitutivePoint point(state);
  Math::SpatialVector<Real> fiber(3);
  fiber.setZero();
  fiber(0) = 1;
  point.set<Tags::FiberDirection>(fiber);
  material.setConstitutivePoint(point);
  Math::SpatialMatrix<Real> before, after;
  material.getFirstPiolaKirchhoffStress(before);
  auto copy = material;
  fiber.setZero();
  fiber(1) = 1;
  point.set<Tags::FiberDirection>(fiber);
  material.setConstitutivePoint(point);
  material.getFirstPiolaKirchhoffStress(after);
  EXPECT_GT((before - after).norm(), 1e-3);
  state.setDisplacementGradient(Math::SpatialMatrix<Real>(Eigen::Matrix3d::Zero()));
  copy.getFirstPiolaKirchhoffStress(after);
  EXPECT_DOUBLE_EQ((before - after).norm(), 0);
}
