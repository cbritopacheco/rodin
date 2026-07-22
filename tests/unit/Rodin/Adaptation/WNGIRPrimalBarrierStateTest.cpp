/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>

#include "Rodin/Adaptation/WNGIRPrimalBarrierState.h"

using namespace Rodin;
using namespace Rodin::Adaptation;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Adaptation_WNGIRPrimalBarrierState, ComputesNewtonCoefficients)
  {
    CellDeformation deformation(2);
    deformation.setDeformationGradient(
      Math::SpatialMatrix<Real>::Identity(2, 2));
    Math::SpatialMatrix<Real> inner =
      Math::SpatialMatrix<Real>::Identity(2, 2);
    inner *= Real(-0.1);
    WNGIRParameters parameters;
    parameters.jSafe = Real(0.1);
    parameters.qMax = Real(3);
    parameters.gammaJ = Real(2);
    parameters.gammaQ = 0;
    parameters.primalBarrierMu = Real(0.01);

    Detail::WNGIRPrimalBarrierState state(
      deformation, inner, parameters, Real(0.01));
    ASSERT_TRUE(state.isFeasible());
    EXPECT_NEAR(state.getJacobianAction(), Real(0.2), 1e-14);
    EXPECT_NEAR(state.getJacobianSlack(), Real(0.7), 1e-14);
    EXPECT_NEAR(state.getJacobianHessian(), Real(0.02 / 0.49), 1e-14);
    EXPECT_NEAR(state.getJacobianForce(), Real(-0.02 / 0.98), 1e-14);
  }

  TEST(Rodin_Adaptation_WNGIRPrimalBarrierState, RejectsNonpositiveTrialSlack)
  {
    CellDeformation deformation(2);
    deformation.setDeformationGradient(
      Math::SpatialMatrix<Real>::Identity(2, 2));
    Math::SpatialMatrix<Real> inner =
      -Math::SpatialMatrix<Real>::Identity(2, 2);
    WNGIRParameters parameters;
    parameters.jSafe = Real(0.1);

    Detail::WNGIRPrimalBarrierState state(
      deformation, inner, parameters, Real(0.01));
    EXPECT_FALSE(state.isFeasible());
  }
}
