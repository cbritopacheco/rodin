/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <sstream>

#include "Rodin/Adaptation.h"
#include "Rodin/Assembly.h"
#include "Rodin/Geometry.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Solver/CG.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Adaptation;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    TEST(Rodin_Adaptation_WNGIRSolver, CalibratedDefaults)
    {
      const WNGIRParameters parameters;
      EXPECT_EQ(parameters.kappaBulk, Real(1e-4));
      EXPECT_EQ(parameters.rigidStabilisationLevel, Real(0.1));
      EXPECT_EQ(parameters.muHat, Real(0.9));
      EXPECT_EQ(parameters.primalBarrierIterations, 15);
      EXPECT_EQ(parameters.primalBarrierRelativeTolerance, Real(1e-3));
      EXPECT_EQ(parameters.cgMaxIterations, 1000);
      EXPECT_EQ(parameters.maxIterations, 30);
      EXPECT_EQ(parameters.geometricValidationOrder, 0);
      EXPECT_EQ(wngirInterfaceQuadratureOrder(1), 4);
      EXPECT_EQ(wngirInterfaceQuadratureOrder(2), 6);
      EXPECT_EQ(wngirInterfaceQuadratureOrder(3), 8);
      EXPECT_EQ(wngirGeometricValidationOrder(1), 6);
      EXPECT_EQ(wngirGeometricValidationOrder(2), 8);
      EXPECT_EQ(wngirGeometricValidationOrder(3), 10);
    }

    TEST(Rodin_Adaptation_WNGIRSolver, GeometricValidationUsesPhysicalDistance)
    {
      WNGIRReport report;
      constexpr Real h = Real(0.02);
      report.effectiveTauRmsH = Real(0.01) * h;
      report.geometricRMS = Real(0.009) * h * h;
      EXPECT_NEAR(report.getGeometricRMSTolerance(h), Real(0.01) * h * h, Real(1e-16));
      EXPECT_TRUE(report.hasGeometricRMSConverged(h));

      report.geometricRMS = Real(0.011) * h * h;
      EXPECT_FALSE(report.hasGeometricRMSConverged(h));
    }

    constexpr Attribute Interface = 10;
    constexpr Attribute FixedBoundary = 11;

    struct SolveState
    {
        Math::Vector<Real> displacement;
        WNGIRReport report;
        Real fixedBoundaryMaximum = 0;
    };

    SolveState solveTranslatedLine(
      Real levelSetScale, Real robustScale = 0, bool freezeBoundary = false, bool trace = false)
    {
      constexpr std::size_t n = 5;
      constexpr Real h = Real(1) / Real(n - 1);
      LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
      mesh.scale(h);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      mesh.getConnectivity().compute(1, 2);

      std::vector<Index> interfaceFacets;
      for (auto face = mesh.getFace(); face; ++face)
      {
        bool onInterface = true;
        for (const Index vertex : face->getVertices())
          onInterface &=
            std::abs(mesh.getVertexCoordinates(vertex)(0) - Real(0.5)) < Real(1e-12);
        if (onInterface)
        {
          interfaceFacets.push_back(face->getIndex());
          mesh.setAttribute({1, face->getIndex()}, Interface);
        }
        bool onBoundary = false;
        for (std::size_t component = 0; component < 2; ++component)
        {
          bool onLowerSide = true;
          bool onUpperSide = true;
          for (const Index vertex : face->getVertices())
          {
            const Real coordinate = mesh.getVertexCoordinates(vertex)(component);
            onLowerSide &= std::abs(coordinate) < Real(1e-12);
            onUpperSide &= std::abs(coordinate - Real(1)) < Real(1e-12);
          }
          onBoundary |= onLowerSide || onUpperSide;
        }
        if (onBoundary)
          mesh.setAttribute({1, face->getIndex()}, FixedBoundary);
      }
      EXPECT_FALSE(interfaceFacets.empty());

      P1<Math::SpatialVector<Real>, LocalMesh> fes(mesh, 2);
      TrialFunction trial(fes);
      TestFunction test(fes);
      WNGIR solver(trial, test);
      WNGIRParameters parameters;
      parameters.h = h;
      parameters.trace = trace;
      parameters.robustScale = robustScale;
      parameters.hasInterfaceAttribute = true;
      parameters.interfaceAttribute = Interface;
      if (freezeBoundary)
      {
        parameters.fixedBoundaryAttributes = {FixedBoundary};
        parameters.rigidStabilisationLevel = 0;
      }
      parameters.maxIterations = 12;
      parameters.tauRms = 0;
      parameters.tauInf = 0;
      // Tight scale-aware tolerances, with the normal-jump raise switched off
      // by zero factors, so the solve is not stopped by the geometric criterion.
      parameters.tauRmsHFloor = Real(1e-3);
      parameters.tauInfHFloor = Real(1e-3);
      parameters.tauJumpRms = 0;
      parameters.tauJumpInf = 0;
      parameters.acceptedStepOverHTol = 0;
      parameters.energyStagTol = 0;
      parameters.quadratureOrder = 2;
      // A tight linear solve, so that the geometric-invariance assertion below
      // measures the invariance and not the conjugate-gradient round-off.
      parameters.cgRelativeTolerance = Real(1e-12);
      solver.setParameters(parameters);

      RealFunction phi([levelSetScale](const Point& point) {
        return levelSetScale * (point.x() - Real(0.55));
      });
      AnalyticVectorFunction grad(
        [levelSetScale](const Point&) {
          Math::SpatialVector<Real> value(2);
          value(0) = levelSetScale;
          value(1) = 0;
          return value;
        },
        2);

      const WNGIRReport report = solver.solve(mesh, interfaceFacets, phi, grad);
      Real fixedBoundaryMaximum = 0;
      if (freezeBoundary)
      {
        const auto& displacement = trial.getSolution().getData();
        for (std::size_t j = 0; j < n; ++j)
        {
          for (std::size_t i = 0; i < n; ++i)
          {
            if (i != 0 && i != n - 1 && j != 0 && j != n - 1)
              continue;
            const auto& dofs = fes.getDOFs(0, j * n + i);
            for (const auto dof : dofs)
              fixedBoundaryMaximum =
                std::max(fixedBoundaryMaximum, std::abs(displacement(dof)));
          }
        }
      }
      return {trial.getSolution().getData(), report, fixedBoundaryMaximum};
    }
  }

  /// @brief Essential constraints keep every boundary displacement exactly zero.
  TEST(Rodin_Adaptation_WNGIRSolver, FreezesMarkedBoundary)
  {
    const SolveState state = solveTranslatedLine(Real(1), Real(0), true);
    EXPECT_LT(state.fixedBoundaryMaximum, Real(1e-12));
  }

  /// @brief The default primal-barrier solve reduces fit while preserving geometry.
  TEST(Rodin_Adaptation_WNGIRSolver, PrimalBarrierFitsAdmissibleTranslation)
  {
    const SolveState state = solveTranslatedLine(Real(1));
    EXPECT_LT(state.report.activeRMS, Real(0.05));
    EXPECT_GT(state.report.minJ, Real(1e-2));
    EXPECT_LT(state.report.maxQRel, Real(10));
    EXPECT_GE(state.report.maxJ, state.report.minJ);
    EXPECT_GT(state.report.activeFraction, Real(0));
    EXPECT_TRUE(std::isfinite(state.report.geometricRMS));
    EXPECT_TRUE(std::isfinite(state.report.geometricSup));
    EXPECT_LT(state.report.normalRMS, Real(1e-12));
    EXPECT_EQ(state.report.rigidModeDimension, 3);
    // The interface is a straight line, so the translation along it is invisible
    // to the rank-one observation and the rigid-mode coercivity constant
    // vanishes. The solve still reaches the fit above because the rigid-mode
    // stabilisation bounds that mode from below.
    EXPECT_LT(state.report.rigidModeCoercivity, Real(1e-12));
    EXPECT_GE(state.report.rigidModeCoercivity, Real(0));
    EXPECT_GT(state.report.iterations, 0);
    EXPECT_GT(state.report.linearSolveCount, 0);
    EXPECT_LE(state.report.maxLinearIterations, 1000);
    EXPECT_TRUE(std::isfinite(state.report.energy));
  }

  /// @brief Tracing preserves numerics and records each accepted geometry.
  TEST(Rodin_Adaptation_WNGIRSolver, TraceRecordsInnerAndAcceptedGeometryWithoutChangingSolve)
  {
    const SolveState base = solveTranslatedLine(Real(1));
    testing::internal::CaptureStdout();
    const SolveState traced = solveTranslatedLine(Real(1), Real(0), false, true);
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(base.report.iterations, traced.report.iterations);
    EXPECT_EQ(base.report.primalBarrierIterations, traced.report.primalBarrierIterations);
    EXPECT_EQ(base.report.energy, traced.report.energy);
    EXPECT_EQ(base.report.geometricSup, traced.report.geometricSup);
    EXPECT_EQ((base.displacement - traced.displacement).norm(), Real(0));
    EXPECT_NE(output.find("barrier inner="), std::string::npos);
    EXPECT_NE(output.find("cg_it="), std::string::npos);
    EXPECT_NE(output.find("rel="), std::string::npos);
    std::istringstream lines(output);
    std::string line;
    std::size_t accepted = 0, initial = 0, final = 0;
    while (std::getline(lines, line))
    {
      if (line.find("wngir geometry:") == std::string::npos)
        continue;
      accepted += line.find("phase=accepted") != std::string::npos;
      initial += line.find("phase=initial") != std::string::npos;
      final += line.find("phase=final") != std::string::npos;
      EXPECT_NE(line.find("inner_total="), std::string::npos);
      EXPECT_NE(line.find("seconds="), std::string::npos);
      EXPECT_NE(line.find("max_qrel="), std::string::npos);
      const std::size_t start = line.find("geom_sup=");
      ASSERT_NE(start, std::string::npos);
      const Real distance = std::stod(line.substr(start + 9));
      EXPECT_TRUE(std::isfinite(distance));
      if (line.find("phase=final") != std::string::npos)
        EXPECT_EQ(distance, traced.report.geometricSup);
    }
    EXPECT_EQ(initial, 1);
    EXPECT_EQ(final, 1);
    EXPECT_EQ(accepted, traced.report.iterations);
  }

  /// @brief Rescaling a level set leaves the geometric displacement unchanged.
  TEST(Rodin_Adaptation_WNGIRSolver, LevelSetScalingIsGeometricallyInvariant)
  {
    const SolveState base = solveTranslatedLine(Real(1));
    const SolveState scaled = solveTranslatedLine(Real(7));
    ASSERT_EQ(base.displacement.size(), scaled.displacement.size());
    EXPECT_NEAR((base.displacement - scaled.displacement).norm(), Real(0), Real(1e-8));
    EXPECT_NEAR(scaled.report.sigma, Real(7) * base.report.sigma, Real(1e-12));
    EXPECT_NEAR(scaled.report.levelSetGradientScale,
      Real(7) * base.report.levelSetGradientScale, Real(1e-12));
    EXPECT_NEAR(
      scaled.report.rigidModeCoercivity, base.report.rigidModeCoercivity, Real(1e-10));
    EXPECT_NEAR(scaled.report.geometricRMS, base.report.geometricRMS, Real(1e-10));
    EXPECT_NEAR(scaled.report.geometricSup, base.report.geometricSup, Real(1e-10));
    EXPECT_NEAR(scaled.report.normalRMS, base.report.normalRMS, Real(1e-12));
    EXPECT_EQ(base.report.iterations, scaled.report.iterations);
    EXPECT_STREQ(base.report.exitReason, scaled.report.exitReason);
  }

  /// @brief A vanishing sampled target gradient is reported before assembly.
  TEST(Rodin_Adaptation_WNGIRSolver, RejectsDegenerateTargetGradient)
  {
    const SolveState state = solveTranslatedLine(Real(0));
    EXPECT_STREQ(state.report.exitReason, "degenerate-target-gradient");
    EXPECT_EQ(state.report.iterations, 0);
    EXPECT_EQ(state.report.levelSetGradientScale, Real(0));
  }

  /// @brief Complete robust saturation is a degeneracy, not a zero residual.
  TEST(Rodin_Adaptation_WNGIRSolver, RejectsEmptyRobustActiveSet)
  {
    const SolveState state = solveTranslatedLine(Real(1), Real(1e-6));
    EXPECT_STREQ(state.report.exitReason, "observation-degenerate-active-set");
    EXPECT_EQ(state.report.iterations, 0);
    EXPECT_EQ(state.report.activeFraction, Real(0));
    EXPECT_TRUE(std::isinf(state.report.activeRMS));
    EXPECT_TRUE(std::isinf(state.report.activeSup));
  }
}
