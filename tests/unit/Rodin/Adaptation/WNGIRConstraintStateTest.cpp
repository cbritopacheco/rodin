#include <gtest/gtest.h>

#include "Rodin/Adaptation/WNGIRConstraintState.h"

using namespace Rodin;
using namespace Rodin::Adaptation;

namespace Rodin::Tests::Unit
{
  namespace
  {
    CellDeformation compressedCell()
    {
      Math::SpatialMatrix<Real> f = Math::SpatialMatrix<Real>::Identity(2, 2);
      f(0, 0) = Real(0.1);
      CellDeformation deformation(2);
      deformation.setDeformationGradient(f);
      return deformation;
    }

    WNGIRParameters jacobianOnlyParameters()
    {
      WNGIRParameters parameters;
      parameters.gammaJ = 1;
      parameters.gammaQ = 0;
      parameters.gammaQual = 0;
      parameters.gammaSize = 0;
      parameters.jSafe = Real(0.01);
      parameters.s0J = Real(0.25);
      return parameters;
    }
  }

  TEST(Rodin_Adaptation_WNGIRConstraintState, SignBlind_UsesBothDirections)
  {
    const auto deformation = compressedCell();
    auto parameters = jacobianOnlyParameters();
    parameters.constraintFormulation =
      WNGIRConstraintFormulation::SignBlindMetric;

    Detail::WNGIRConstraintState state(deformation, parameters);
    EXPECT_GT(state.getJacobianWeight(), 0);
  }

  TEST(Rodin_Adaptation_WNGIRConstraintState, Directional_RejectsImprovingAction)
  {
    const auto deformation = compressedCell();
    auto parameters = jacobianOnlyParameters();
    parameters.constraintFormulation =
      WNGIRConstraintFormulation::DirectionalMetric;
    Math::SpatialMatrix<Real> improving = Math::SpatialMatrix<Real>::Identity(2, 2);
    improving.setZero();
    improving(0, 0) = 1;
    Math::SpatialMatrix<Real> consuming = -improving;

    Detail::WNGIRConstraintState improvingState(
      deformation, parameters, &improving);
    Detail::WNGIRConstraintState consumingState(
      deformation, parameters, &consuming);
    EXPECT_EQ(improvingState.getJacobianWeight(), 0);
    EXPECT_GT(consumingState.getJacobianWeight(), 0);
  }

  TEST(Rodin_Adaptation_WNGIRConstraintState, Fractional_ActivatesAtFractionalMargin)
  {
    const auto deformation = compressedCell();
    auto parameters = jacobianOnlyParameters();
    parameters.constraintFormulation =
      WNGIRConstraintFormulation::FractionalMarginMetric;
    parameters.marginFraction = Real(0.5);
    Math::SpatialMatrix<Real> below = Math::SpatialMatrix<Real>::Identity(2, 2);
    below.setZero();
    below(0, 0) = Real(-0.01);
    Math::SpatialMatrix<Real> above = Math::SpatialMatrix<Real>::Identity(2, 2);
    above.setZero();
    above(0, 0) = Real(-0.1);

    Detail::WNGIRConstraintState belowState(deformation, parameters, &below);
    Detail::WNGIRConstraintState aboveState(deformation, parameters, &above);
    EXPECT_EQ(belowState.getJacobianWeight(), 0);
    EXPECT_GT(aboveState.getJacobianWeight(), 0);
    EXPECT_NEAR(aboveState.getJacobianTarget(), Real(0.045), 1e-14);
  }

  TEST(Rodin_Adaptation_WNGIRConstraintState, Safeguarded_UsesFractionalLocalModel)
  {
    const auto deformation = compressedCell();
    auto parameters = jacobianOnlyParameters();
    parameters.constraintFormulation =
      WNGIRConstraintFormulation::SafeguardedMarginMetric;
    parameters.marginFraction = Real(0.25);
    Math::SpatialMatrix<Real> predictor =
      Math::SpatialMatrix<Real>::Identity(2, 2);
    predictor.setZero();
    predictor(0, 0) = Real(-0.1);

    Detail::WNGIRConstraintState state(deformation, parameters, &predictor);
    EXPECT_GT(state.getJacobianWeight(), 0);
    EXPECT_NEAR(state.getJacobianTarget(), Real(0.0225), 1e-14);
  }
}
