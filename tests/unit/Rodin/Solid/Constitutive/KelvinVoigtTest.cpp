/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Finite-strain viscosity, power, and consistent time tangent tests.
#include <gtest/gtest.h>
#include <Eigen/Geometry>
#include "Rodin/Solid/Constitutive/DemirayGasser.h"
#include "Rodin/Solid/Constitutive/NeoHookean.h"
#include "Rodin/Solid/Constitutive/KelvinVoigt.h"

using namespace Rodin;
using namespace Rodin::Solid;

namespace
{
  struct Response
  {
    Real dissipation;
    Math::SpatialMatrix<Real> stress, tangent;
  };

  Response evaluate(Real viscosity, const Eigen::Matrix3d& f, const Eigen::Matrix3d& g,
    const Eigen::Matrix3d& direction, Real coefficient)
  {
    KinematicState state(3);
    state.setDisplacementGradient(
      Math::SpatialMatrix<Real>(f - Eigen::Matrix3d::Identity()));
    ConstitutivePoint point(state);
    point.set<Tags::VelocityGradient>(Math::SpatialMatrix<Real>(g));
    point.set<Tags::VelocityGradientDerivative>(coefficient);
    KelvinVoigt material(DemirayGasser(DemirayGasser::Parameters{}), viscosity);
    material.setConstitutivePoint(point);
    Response result;
    material.getFirstPiolaKirchhoffStress(result.stress);
    material.getMaterialTangent(result.tangent, Math::SpatialMatrix<Real>(direction));
    result.dissipation = material.getDissipationRate();
    return result;
  }
}

TEST(Rodin_Solid_KelvinVoigt, IndependentStressPowerAndFullTimeTangent)
{
  Eigen::Matrix3d f;
  f << 1.1, .03, 0, .02, .95, .04, 0, .01, 1.02;
  Eigen::Matrix3d g;
  g << .2, .1, -.03, .04, -.1, .02, .03, .07, .08;
  const Eigen::Matrix3d h = .1 * Eigen::Matrix3d::Ones();
  constexpr Real viscosity = 2.3, coefficient = 17, step = 1e-6;
  const auto response = evaluate(viscosity, f, g, h, coefficient);
  const auto elastic = evaluate(0, f, g, h, coefficient);
  const Eigen::Matrix3d rate = .5 * (f.transpose() * g + g.transpose() * f);
  const Eigen::Matrix3d expected = viscosity * f * rate;
  EXPECT_LT(
    (response.stress - elastic.stress - Math::SpatialMatrix<Real>(expected)).norm(),
    1e-12);
  EXPECT_NEAR(response.dissipation, viscosity * rate.squaredNorm(), 1e-12);
  EXPECT_GE(response.dissipation, 0);
  const auto viscous = response.stress - elastic.stress;
  EXPECT_NEAR(viscous.dot(Math::SpatialMatrix<Real>(g)), response.dissipation, 1e-12);
  const auto plus =
    evaluate(viscosity, f + step * h, g + step * coefficient * h, h, coefficient);
  const auto minus =
    evaluate(viscosity, f - step * h, g - step * coefficient * h, h, coefficient);
  EXPECT_LT((response.tangent - (plus.stress - minus.stress) / (2 * step)).norm(), 2e-6);
}

TEST(Rodin_Solid_KelvinVoigt, RigidRotationHasNoViscousStress)
{
  const Eigen::Matrix3d f =
    Eigen::AngleAxisd(.7, Eigen::Vector3d::UnitY()).toRotationMatrix();
  Eigen::Matrix3d spin;
  spin << 0, -2, 1, 2, 0, -3, -1, 3, 0;
  const auto actual = evaluate(3, f, spin * f, Eigen::Matrix3d::Zero(), 0);
  const auto elastic = evaluate(0, f, spin * f, Eigen::Matrix3d::Zero(), 0);
  EXPECT_LT((actual.stress - elastic.stress).norm(), 1e-12);
  EXPECT_LT(actual.dissipation, 1e-26);
}

TEST(Rodin_Solid_KelvinVoigt, InvalidViscosityAndMissingRateInputs)
{
  const DemirayGasser elastic(DemirayGasser::Parameters{});
  EXPECT_THROW((KelvinVoigt<DemirayGasser>(elastic, -1)), std::exception);
  KelvinVoigt material(elastic, 1);
  KinematicState state(3);
  ConstitutivePoint point(state);
  EXPECT_THROW(material.setConstitutivePoint(point), std::exception);
  KelvinVoigt zero(elastic, 0);
  EXPECT_NO_THROW(zero.setConstitutivePoint(point));
  EXPECT_DOUBLE_EQ(zero.getDissipationRate(), 0);
}

TEST(Rodin_Solid_KelvinVoigt, SuperposedObserverSpinPreservesDissipation)
{
  Eigen::Matrix3d f = Eigen::Matrix3d::Identity();
  f(0, 0) = 1.08;
  f(1, 2) = .03;
  Eigen::Matrix3d g = .1 * Eigen::Matrix3d::Identity();
  g(0, 1) = .04;
  const Eigen::Matrix3d q =
    Eigen::AngleAxisd(.6, Eigen::Vector3d::UnitZ()).toRotationMatrix();
  Eigen::Matrix3d spin;
  spin << 0, -2, 0, 2, 0, -1, 0, 1, 0;
  const Eigen::Matrix3d h = Eigen::Matrix3d::Zero();
  const auto a = evaluate(3, f, g, h, 0),
             b = evaluate(3, q * f, q * g + spin * q * f, h, 0);
  EXPECT_NEAR(a.dissipation, b.dissipation, 1e-12);
  EXPECT_LT((b.stress - Math::SpatialMatrix<Real>(q) * a.stress).norm(), 1e-10);
}

TEST(Rodin_Solid_KelvinVoigt, RebindingRatesAndFailureInvalidatePreparedState)
{
  KelvinVoigt material(DemirayGasser(DemirayGasser::Parameters{}), 2);
  KinematicState state(3);
  ConstitutivePoint point(state);
  point.set<Tags::VelocityGradient>(Math::SpatialMatrix<Real>(Eigen::Matrix3d::Zero()));
  point.set<Tags::VelocityGradientDerivative>(1);
  material.setConstitutivePoint(point);
  EXPECT_DOUBLE_EQ(material.getDissipationRate(), 0);
  point.set<Tags::VelocityGradient>(Math::SpatialMatrix<Real>::Identity(3, 3));
  material.setConstitutivePoint(point);
  EXPECT_NEAR(material.getDissipationRate(), 6, 1e-12);
  point.set<Tags::VelocityGradientDerivative>(-1);
  EXPECT_THROW(material.setConstitutivePoint(point), std::exception);
  EXPECT_THROW(material.getStrainEnergyDensity(), std::exception);
}

TEST(Rodin_Solid_KelvinVoigt, TwoDimensionalRateAndTimeTangent)
{
  KelvinVoigt material(NeoHookean(1, 1), 2);
  KinematicState state(2);
  ConstitutivePoint point(state);
  const auto identity = Math::SpatialMatrix<Real>::Identity(2, 2);
  point.set<Tags::VelocityGradient>(identity);
  point.set<Tags::VelocityGradientDerivative>(3);
  material.setConstitutivePoint(point);
  Math::SpatialMatrix<Real> stress, tangent;
  material.getFirstPiolaKirchhoffStress(stress);
  material.getMaterialTangent(tangent, identity);
  EXPECT_LT((stress - 2 * identity).norm(), 1e-13);
  EXPECT_LT((tangent - 14 * identity).norm(), 1e-13);
  EXPECT_DOUBLE_EQ(material.getDissipationRate(), 4);
}
