/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <vector>

#include <gtest/gtest.h>

#include "TMOPFDTestHelpers.h"

using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  TEST(Rodin_Adaptation_TargetMatrixOptimization_PerClass,
       TMOPResidualIsPositiveEnergyGradient)
  {
    auto mesh = Rodin::Tests::Unit::TMOPFD::makeTwoTriangleSquareWithVerticalInterface();
    P1 space(mesh, 2);
    GridFunction u(space);
    Rodin::Tests::Unit::TMOPFD::fillDisplacement(u.getData(), Rodin::Real(0.002));
    TestFunction v(space);

    auto quality = Rodin::Tests::Unit::TMOPFD::makeQualityTerm();
    auto deviation = Rodin::Tests::Unit::TMOPFD::makeDeviationTerm();
    auto fit = Rodin::Tests::Unit::TMOPFD::makeAnalyticFitTerm();
    auto phase = Rodin::Tests::Unit::TMOPFD::makePhaseTerm();

    auto assembleResidual = [&]()
    {
      LinearForm r(v);
      r = quality.residual(u, v)
        + deviation.residual(u, v)
        + fit.residual(u, v)
        + phase.residual(u, v);
      r.assemble();
      return r.getVector();
    };
    auto energy = [&]()
    {
      return quality.energy(u)
           + deviation.energy(u)
           + fit.energy(u)
           + phase.energy(u);
    };

    const auto u0 = u.getData();
    auto direction = u0;
    Rodin::Tests::Unit::TMOPFD::fillDirection(direction);
    const auto g = assembleResidual();
    const Rodin::Real directionalResidual = g.dot(direction);

    const Rodin::Real eps = 1e-6;
    u.getData() = u0 + eps * direction;
    const Rodin::Real ePlus = energy();
    u.getData() = u0 - eps * direction;
    const Rodin::Real eMinus = energy();
    u.getData() = u0;

    const Rodin::Real fd = (ePlus - eMinus) / (2 * eps);
    const Rodin::Real denom =
      std::max<Rodin::Real>({1, std::abs(fd), std::abs(directionalResidual)});
    EXPECT_NEAR(fd, directionalResidual, 1e-7 * denom);
    EXPECT_GT(fd * directionalResidual, 0);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization_PerClass,
       TMOPMinimizerAcceptsOnlyMonotoneEnergySteps)
  {
    auto mesh = Rodin::Tests::Unit::TMOPFD::makeTwoTriangleSquareWithVerticalInterface();
    RealH1Element<2> fe(Polytope::Type::Triangle);
    upgradeTransformations(mesh, fe, Rodin::Tests::Unit::TMOPFD::Interface);

    VectorH1<2, LocalMesh> space(
        std::integral_constant<size_t, 2>{}, mesh, 2);
    GridFunction u(space);
    u.getData().setZero();
    TrialFunction p(space);
    TestFunction v(space);

    CurvedQualityTargetJacobian target(mesh, 0.25);
    QualityTerm quality(ShapeSizeBlendMetric(0.5), target, 0.5);
    quality.setQuadratureOrder(4);
    DeviationTerm deviation(0.2);
    auto fit = Rodin::Tests::Unit::TMOPFD::makeAnalyticFitTerm();
    auto phase = Rodin::Tests::Unit::TMOPFD::makePhaseTerm();

    auto makeResidual = [&]()
    {
      return quality.residual(u, v)
           + deviation.residual(u, v)
           + fit.residual(u, v)
           + phase.residual(u, v);
    };
    auto energy = [&]()
    {
      return quality.energy(u)
           + deviation.energy(u)
           + fit.energy(u)
           + phase.energy(u);
    };
    auto admissible = [&]()
    {
      return isTargetAdmissible(
          mesh,
          u,
          fe,
          target,
          Rodin::Real(0),
          Polytope::Type::Triangle,
          4);
    };

    std::vector<Rodin::Real> acceptedEnergies;
    IsoparametricTMOPMinimizerParameters params;
    params.maxIterations = 8;
    params.gradientTolerance = 1e-14;
    params.stepTolerance = 1e-14;
    params.energyTolerance = 0;
    params.preconditionerLength = 0.25;
    params.acceptedEnergyMonitor = [&](Rodin::Real e)
    {
      acceptedEnergies.push_back(e);
    };

    const auto report = minimizeIsoparametricTMOP(
        mesh,
        fe,
        u,
        p,
        v,
        makeResidual,
        energy,
        admissible,
        Rodin::Tests::Unit::TMOPFD::Interface,
        params);

    ASSERT_GT(report.acceptedSteps, 0);
    ASSERT_EQ(acceptedEnergies.size(), static_cast<size_t>(report.acceptedSteps));
    Rodin::Real previous = report.initialEnergy;
    for (const Rodin::Real e : acceptedEnergies)
    {
      EXPECT_LE(e, previous + 1e-12);
      previous = e;
    }
    EXPECT_LE(report.finalEnergy, report.initialEnergy);
    EXPECT_TRUE(std::isfinite(report.finalGradientNorm));
  }
}
