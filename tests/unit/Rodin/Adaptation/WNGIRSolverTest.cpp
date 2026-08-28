/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>

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
    constexpr Attribute Interface = 10;

    struct SolveState
    {
        Math::Vector<Real> displacement;
        WNGIRReport report;
    };

    SolveState solveTranslatedLine(Real levelSetScale, Real robustScale = 0)
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
      }
      EXPECT_FALSE(interfaceFacets.empty());

      P1<Math::SpatialVector<Real>, LocalMesh> fes(mesh, 2);
      TrialFunction trial(fes);
      TestFunction test(fes);
      WNGIR solver(trial, test);
      WNGIRParameters parameters;
      parameters.h = h;
      parameters.robustScale = robustScale;
      parameters.hasInterfaceAttribute = true;
      parameters.interfaceAttribute = Interface;
      parameters.maxIterations = 12;
      parameters.activeRMSTol = 0;
      parameters.activeSupTol = 0;
      parameters.activeRMSOverHTol = Real(1e-3);
      parameters.activeSupOverHTol = Real(1e-3);
      parameters.geometryAwareTolerances = false;
      parameters.acceptedStepOverHTol = 0;
      parameters.energyStagTol = 0;
      parameters.quadratureOrder = 2;
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
      return {trial.getSolution().getData(), report};
    }
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
    EXPECT_EQ(state.report.rigidModeDimension, 3);
    EXPECT_GT(state.report.rigidModeCoercivity, Real(0));
    EXPECT_GT(state.report.rigidModeCoercivityRatio, Real(0));
    EXPECT_GT(state.report.iterations, 0);
    EXPECT_TRUE(std::isfinite(state.report.energy));
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
