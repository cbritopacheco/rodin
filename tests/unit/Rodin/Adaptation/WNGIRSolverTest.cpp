/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <functional>
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
    TEST(Rodin_Adaptation_WNGIRSolver, CanonicalDefaults)
    {
      const WNGIRParameters parameters;
      EXPECT_EQ(parameters.kappaF, Real(1));
      EXPECT_EQ(parameters.muHat, Real(90));
      EXPECT_EQ(parameters.primalBarrierIterations, 15);
      EXPECT_EQ(parameters.primalBarrierRelativeTolerance, Real(1e-3));
      EXPECT_EQ(parameters.cgMaxIterations, 1000);
      EXPECT_EQ(parameters.kappaS, Real(1));
      EXPECT_EQ(parameters.kappaD, Real(1));
      EXPECT_EQ(parameters.directSolverThreads, 0u);
      EXPECT_EQ(parameters.maxIterations, 30);
      EXPECT_TRUE(parameters.directionalNewton);
      EXPECT_EQ(parameters.directionalNewtonMaxStepOverH, Real(1));
      EXPECT_EQ(parameters.qualityGuard, Real(0.1));
      EXPECT_EQ(parameters.geometricSupTolerance, Real(0));
      EXPECT_EQ(parameters.geometricValidationOrder, 0);
      EXPECT_EQ(wngirInterfaceQuadratureOrder(1), 12);
      EXPECT_EQ(wngirInterfaceQuadratureOrder(2), 12);
      EXPECT_EQ(wngirInterfaceQuadratureOrder(3), 12);
      EXPECT_EQ(wngirInterfaceQuadratureOrder(7), 16);
      EXPECT_EQ(wngirGeometricValidationOrder(1), 14);
      EXPECT_EQ(wngirGeometricValidationOrder(2), 14);
      EXPECT_EQ(wngirGeometricValidationOrder(3), 14);
    }

    TEST(Rodin_Adaptation_WNGIRSolver, CurvedSurfaceQuadratureResolvesNonpolynomialFit)
    {
      const auto integrate = [](size_t order) {
        const auto& rule =
          QF::PolytopeQuadratureFormula::get(order, Polytope::Type::Segment);
        const WNGIRLoss loss(Real(0.31927236562898026));
        Math::Vector<Real> integral = Math::Vector<Real>::Zero(3);
        for (size_t q = 0; q < rule.getSize(); ++q)
        {
          const Real t = rule.getPoint(q)(0);
          const Real x = Real(0.73092457) * (1 - t) + Real(0.60939475) * t - Real(0.5);
          const Real y = Real(0.59087505) * (1 - t) + Real(0.60939475) * t - Real(0.5);
          const Real r = std::hypot(x, y), theta = std::atan2(y, x);
          const Real residual = r - Real(0.24) - Real(0.08) * std::cos(Real(4) * theta);
          const Real angular = Real(0.32) * std::sin(Real(4) * theta) / (r * r);
          integral(0) += rule.getWeight(q) * loss.getValue(residual);
          integral(1) +=
            rule.getWeight(q) * loss.getInfluence(residual) * (x / r - angular * y);
          integral(2) +=
            rule.getWeight(q) * loss.getInfluence(residual) * (y / r + angular * x);
        }
        return integral;
      };
      const auto reference = integrate(24),
                 actual = integrate(wngirInterfaceQuadratureOrder(1));
      EXPECT_LT((actual - reference).norm() / reference.norm(), Real(1e-4));
      EXPECT_LT(std::abs(actual(0) - reference(0)) / reference(0), Real(1e-4));
      EXPECT_GT(std::abs(integrate(4)(0) - reference(0)) / reference(0), Real(0.1));
    }

    TEST(Rodin_Adaptation_WNGIRSolver, DirectionalNewtonParameters)
    {
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
      P1<Math::SpatialVector<Real>, LocalMesh> fes(mesh, 2);
      TrialFunction trial(fes);
      TestFunction test(fes);
      WNGIR solver(trial, test);
      WNGIRParameters p;
      EXPECT_NO_THROW(solver.setParameters(p));
      for (const Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(),
             std::numeric_limits<Real>::quiet_NaN()})
      {
        p.directionalNewtonMaxStepOverH = invalid;
        EXPECT_THROW(solver.setParameters(p), Alert::Exception);
      }
      p.directionalNewtonMaxStepOverH = Real(0.5);
      EXPECT_NO_THROW(solver.setParameters(p));
      p.directionalNewton = false;
      EXPECT_NO_THROW(solver.setParameters(p));
    }

    constexpr Attribute Interface = 10;

    struct SolveState
    {
        Math::Vector<Real> displacement;
        WNGIRReport report;
        Real physicalNorm = 0;
    };

    template <std::size_t Order = 1>
    SolveState solveTranslatedLine(Real levelSetScale, Real robustScale = 0,
      bool trace = false, Real target = 0, bool partialGradient = false,
      std::size_t cgCap = 1000, bool strictCG = false,
      WNGIRParameters::DirectSolver directSolver =
        WNGIRParameters::DirectSolver::SparseLU,
      bool flat = false, Real innerTolerance = Real(1e-3), Real muHat = Real(90),
      const std::function<void(WNGIRParameters&)>& configure = {})
    {
      constexpr std::size_t n = 5;
      constexpr Real h = Real(1) / Real(n - 1);
      LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
      mesh.scale(h);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      mesh.getConnectivity().compute(1, 2);
      if constexpr (Order > 1)
        for (std::size_t from = 1; from <= 2; ++from)
          for (std::size_t to = 0; to <= 2; ++to)
            if (from != to)
              mesh.getConnectivity().compute(from, to);

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

      auto fes = [&] {
        if constexpr (Order == 1)
          return P1<Math::SpatialVector<Real>, LocalMesh>(mesh, 2);
        else
          return H1(std::integral_constant<std::size_t, Order>{}, mesh, 2);
      }();
      TrialFunction trial(fes);
      TestFunction test(fes);
      WNGIR solver(trial, test);
      WNGIRParameters parameters;
      parameters.h = h;
      parameters.trace = trace;
      parameters.traceQualityWitness = trace;
      parameters.geometricSupTolerance = target == Real(0) ? Real(1e-3) : target;
      parameters.directSolver = directSolver;
      parameters.directSolverThreads = 2;
      parameters.primalBarrierRelativeTolerance = innerTolerance;
      parameters.muHat = muHat;
      parameters.robustScale = robustScale;
      parameters.hasInterfaceAttribute = true;
      parameters.interfaceAttribute = Interface;
      parameters.maxIterations = 12;
      parameters.acceptedStepOverHTol = 0;
      parameters.energyStagTol = 0;
      parameters.quadratureOrder = 2;
      // A tight linear solve, so that the geometric-invariance assertion below
      // measures the invariance and not the conjugate-gradient round-off.
      parameters.cgRelativeTolerance = Real(1e-10);
      parameters.cgMaxIterations = cgCap;
      if (configure)
        configure(parameters);
      solver.setParameters(parameters);

      RealFunction phi([levelSetScale, flat](const Point& point) {
        return levelSetScale *
          (point.x() - Real(0.55) +
            (flat ? Real(0) : Real(0.02) * std::sin(Real(6) * point.y())));
      });
      AnalyticVectorFunction grad(
        [levelSetScale, partialGradient, flat](const Point& point) {
          Math::SpatialVector<Real> value(2);
          value(0) = partialGradient && point.y() < Real(0.4) ? Real(0) : levelSetScale;
          value(1) = flat || (partialGradient && point.y() < Real(0.4))
            ? Real(0)
            : levelSetScale * Real(0.12) * std::cos(Real(6) * point.y());
          return value;
        },
        2);

      const WNGIRReport report = solver.solve(mesh, interfaceFacets, phi, grad);
      Real physicalNorm = 0;
      const auto& solution = trial.getSolution();
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const Polytope::Traits traits(cell->getGeometry());
        const auto measure = [&](const Point& point) {
          const auto value = solution.getValue(point);
          for (std::size_t component = 0; component < 2; ++component)
            physicalNorm = std::max(physicalNorm, std::abs(value(component)));
        };
        for (std::size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
          measure(Point(*cell, traits.getVertex(vertex)));
        const auto& qf = QF::PolytopeQuadratureFormula::get(
          wngirGeometricValidationOrder(Order), cell->getGeometry());
        const auto& quadrature = cell->getQuadrature(qf);
        for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          measure(quadrature.getPoint(q));
      }
      return {trial.getSolution().getData(), report, physicalNorm};
    }
  }

  TEST(Rodin_Adaptation_WNGIRSolver, SmallDirectionCanProduceLargeAcceptedStep)
  {
    constexpr Real tolerance = Real(0.01);
    const auto state = solveTranslatedLine(Real(1), 0, false, 0, false, 1000,
      true, WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), Real(0),
      [](WNGIRParameters& p) {
        p.kappaF = Real(100);
        p.stepTol = tolerance;
        p.maxIterations = 1;
      });
    EXPECT_EQ(state.report.iterations, 1u);
    EXPECT_LE(state.report.lastAlpha, Real(1));
    EXPECT_GT(state.report.predictorScale, Real(1));
    EXPECT_LT(state.report.acceptedStep / state.report.predictorScale, tolerance);
    EXPECT_GT(state.report.acceptedStep, tolerance);
    EXPECT_STRNE(state.report.exitReason, "best-effort-step-stagnation");
    EXPECT_GT(state.report.minJ, Real(0.01));
    EXPECT_LT(state.report.maxQRel, Real(10));
  }

  TEST(Rodin_Adaptation_WNGIRSolver, AbsoluteStagnationUsesAcceptedStep)
  {
    constexpr Real tolerance = Real(1);
    const auto state = solveTranslatedLine(Real(1), 0, false, 0, false, 1000,
      true, WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), Real(0),
      [](WNGIRParameters& p) {
        p.kappaF = Real(100);
        p.stepTol = tolerance;
        p.stagnationIterations = 1;
        p.geometricSupTolerance = Real(1e-12);
      });
    EXPECT_EQ(state.report.iterations, 1u);
    EXPECT_GT(state.report.acceptedStep, Real(0));
    EXPECT_LE(state.report.acceptedStep, tolerance);
    EXPECT_STREQ(state.report.exitReason, "best-effort-step-stagnation");
    EXPECT_EQ((state.displacement.cwiseAbs().maxCoeff()), state.report.acceptedStep);
  }

  TEST(Rodin_Adaptation_WNGIRSolver, GeometricTargetPrecedesAbsoluteStagnation)
  {
    const auto state = solveTranslatedLine(Real(1), 0, false, Real(0.06), false,
      1000, true, WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), Real(0),
      [](WNGIRParameters& p) {
        p.kappaF = Real(100);
        p.stepTol = Real(1);
      });
    EXPECT_LE(state.report.iterations, 1u);
    EXPECT_LE(state.report.acceptedStep, Real(1));
    EXPECT_LE(state.report.geometricSup, Real(0.06));
    EXPECT_STREQ(state.report.exitReason, "full-interface-geometric-sup-converged");
  }

  TEST(Rodin_Adaptation_WNGIRSolver, SquaredObservationDropsOnlyLevelSetHessian)
  {
    for (const size_t dimension : {2u, 3u})
    {
      auto mesh = dimension == 2
        ? LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3})
        : LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
      mesh.scale(Real(0.5));
      P1<Math::SpatialVector<Real>, LocalMesh> fes(mesh, dimension);
      GridFunction current(fes);
      current.getData().setZero();
      const Location::AABB<LocalMesh> locator(mesh);
      auto cell = mesh.getCell();
      const auto& qf = QF::PolytopeQuadratureFormula::get(2, cell->getGeometry());
      const auto& point = cell->getQuadrature(qf).getPoint(0);
      const IntegrationPoint ip(point, &qf, 0);
      WNGIRParameters p;
      for (const Real levelSetScale : {Real(1), Real(7)})
      {
        RealFunction phi([levelSetScale](const Point& x) {
          return levelSetScale * (Real(2) + x.getCoordinates().squaredNorm());
        });
        AnalyticVectorFunction grad(
          [levelSetScale](const Point& x) {
            return Math::SpatialVector<Real>(
              Real(2) * levelSetScale * x.getCoordinates());
          },
          dimension);
        const Real normalization = Real(1) / (levelSetScale * levelSetScale);
        Detail::WNGIRObservationCoefficient coefficient(
          grad, current, locator, p, normalization, dimension);
        const auto actual = coefficient.getValue(ip);
        const auto x = point.getCoordinates();
        Math::SpatialMatrix<Real> expected(dimension, dimension);
        for (size_t i = 0; i < dimension; ++i)
          for (size_t j = 0; j < dimension; ++j)
            expected(i, j) = Real(4) * x(i) * x(j);
        EXPECT_LT((actual - expected).norm(), Real(1e-12));
        Math::SpatialVector<Real> v(dimension);
        for (size_t i = 0; i < dimension; ++i)
          v(i) = Real(1);
        Math::SpatialVector<Real> z = v;
        z(0) = Real(0.3);
        const auto energy = [&](const Math::SpatialVector<Real>& y) {
          const Real residual = levelSetScale * (Real(2) + y.squaredNorm());
          return Real(0.5) * normalization * residual * residual;
        };
        constexpr Real eps = Real(1e-4);
        const Real mixed = (energy(Math::SpatialVector<Real>(x + eps * v + eps * z)) -
                             energy(Math::SpatialVector<Real>(x + eps * v - eps * z)) -
                             energy(Math::SpatialVector<Real>(x - eps * v + eps * z)) +
                             energy(Math::SpatialVector<Real>(x - eps * v - eps * z))) /
          (Real(4) * eps * eps);
        const Real omitted = Real(2) * (Real(2) + x.squaredNorm()) * v.dot(z);
        EXPECT_NEAR(v.dot(actual * z), mixed - omitted, Real(2e-6));
        p.kappaF = Real(0.5);
        EXPECT_LT((coefficient.getValue(ip) - Real(0.5) * expected).norm(), Real(1e-12));
        p.kappaF = Real(1);
      }
    }
  }

  /// @brief The default primal-barrier solve reduces fit while preserving geometry.
  TEST(Rodin_Adaptation_WNGIRSolver, FitsAdmissibleWavyTarget)
  {
    const SolveState state = solveTranslatedLine(Real(1));
    EXPECT_LT(state.report.activeRMS, Real(0.05));
    EXPECT_GT(state.report.minJ, Real(1e-2));
    EXPECT_LT(state.report.maxQRel, Real(10));
    EXPECT_GE(state.report.maxJ, state.report.minJ);
    EXPECT_GT(state.report.activeFraction, Real(0));
    EXPECT_TRUE(std::isfinite(state.report.geometricRMS));
    EXPECT_TRUE(std::isfinite(state.report.geometricSup));
    EXPECT_TRUE(std::isfinite(state.report.normalRMS));
    EXPECT_GT(state.report.iterations, 0);
    EXPECT_GT(state.report.linearSolveCount, 0);
    EXPECT_LE(state.report.maxLinearIterations, 1000);
    EXPECT_TRUE(std::isfinite(state.report.energy));
  }

  /// @brief Tracing preserves numerics and records each accepted geometry.
  TEST(Rodin_Adaptation_WNGIRSolver,
    TraceRecordsInnerAndAcceptedGeometryWithoutChangingSolve)
  {
    const SolveState base = solveTranslatedLine(Real(1));
    testing::internal::CaptureStdout();
    const SolveState traced = solveTranslatedLine(Real(1), Real(0), true);
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(base.report.iterations, traced.report.iterations);
    EXPECT_EQ(base.report.primalBarrierIterations, traced.report.primalBarrierIterations);
    EXPECT_EQ(base.report.energy, traced.report.energy);
    EXPECT_EQ(base.report.geometricSup, traced.report.geometricSup);
    EXPECT_EQ((base.displacement - traced.displacement).norm(), Real(0));
    EXPECT_TRUE(output.find("barrier inner=") != std::string::npos ||
      output.find("barrier skip:") != std::string::npos);
    EXPECT_NE(output.find("wngir directional:"), std::string::npos);
    EXPECT_NE(output.find("wngir quality witness:"), std::string::npos);
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

  /// @brief Target stopping uses the whole interface, not only the robust active set.

  /// @brief Target stopping uses the whole interface, not only the robust active set.
  TEST(Rodin_Adaptation_WNGIRSolver, FullInterfaceTarget)
  {
    const auto initial = solveTranslatedLine(Real(1), 0, false, Real(0.1));
    EXPECT_EQ(initial.report.iterations, 0);
    EXPECT_LT(initial.report.geometricSup, Real(0.1));
    EXPECT_STREQ(initial.report.exitReason, "full-interface-geometric-sup-converged");
    const auto fitted = solveTranslatedLine(Real(1), 0, false, Real(1e-3));
    EXPECT_GT(fitted.report.iterations, 0);
    EXPECT_LE(fitted.report.geometricSup, Real(1e-3));
    EXPECT_GT(fitted.report.minJ, Real(0.01));
    EXPECT_LT(fitted.report.maxQRel, Real(10));
    EXPECT_STREQ(fitted.report.exitReason, "full-interface-geometric-sup-converged");
    EXPECT_THROW(solveTranslatedLine(Real(1), 0, false, Real(-1)), Alert::Exception);
    const auto invalid = solveTranslatedLine(Real(1), 0, false, Real(0.1), true);
    // Invalid gradients on part of the initial interface forbid an initial hit.
    EXPECT_FALSE(invalid.report.geometricTargetReached);
    EXPECT_FALSE(std::isfinite(invalid.report.geometricSup));
  }

  TEST(Rodin_Adaptation_WNGIRSolver, GlobalizedInnerMeritDecreases)
  {
    {
      testing::internal::CaptureStdout();
      const SolveState state =
        solveTranslatedLine(Real(1), Real(0), true, Real(0), false, 1000, true);
      const std::string output = testing::internal::GetCapturedStdout();
      EXPECT_LT(state.report.geometricSup, Real(0.05));
      EXPECT_GT(state.report.minJ, Real(0.01));
      EXPECT_LT(state.report.maxQRel, Real(10));
      EXPECT_GT(state.report.tPrimalBarrierLineSearch, Real(0));
      std::istringstream lines(output);
      std::string line;
      std::size_t merits = 0;
      while (std::getline(lines, line))
      {
        if (line.find("barrier merit:") == std::string::npos)
          continue;
        ++merits;
        const Real before = std::stod(line.substr(line.find("before=") + 7));
        const Real after = std::stod(line.substr(line.find("after=") + 6));
        EXPECT_TRUE(std::isfinite(before));
        EXPECT_TRUE(std::isfinite(after));
        EXPECT_LE(after, before + Real(1e-12) * std::max(std::abs(before), Real(1e-30)));
        EXPECT_NE(line.find("accepted=1"), std::string::npos);
      }
      EXPECT_GT(merits, 0u);
    }
  }

  TEST(Rodin_Adaptation_WNGIRSolver, CanonicalMetricIsFrozenWithoutCompletion)
  {
    const auto state =
      solveTranslatedLine(Real(1), Real(0), false, Real(0), false, 1000, true);
    EXPECT_LT(state.report.geometricSup, Real(0.05));
    EXPECT_GT(state.report.minJ, Real(0.01));
    EXPECT_LT(state.report.maxQRel, Real(10));
    EXPECT_TRUE(state.report.primalBarrierConverged);
  }

  /// @brief Sparse augmented LU retains the same low-rank rigid operator as CG.
  TEST(Rodin_Adaptation_WNGIRSolver, DirectAugmentedSolvePreservesSimilarityGauge)
  {
    const auto solve = [](WNGIRParameters::DirectSolver backend) {
      return solveTranslatedLine(
        Real(1), Real(0), false, Real(0), false, 1000, true, backend);
    };
    const auto cg = solve(WNGIRParameters::DirectSolver::CG);
    const auto lu = solve(WNGIRParameters::DirectSolver::SparseLU);
    EXPECT_EQ(cg.report.iterations, lu.report.iterations);
    EXPECT_STREQ(cg.report.exitReason, lu.report.exitReason);
    EXPECT_NEAR((cg.displacement - lu.displacement).norm(), Real(0), Real(1e-8));
    EXPECT_EQ(lu.report.linearIterations, 0u);
    EXPECT_GT(lu.report.linearSolveCount, 0u);
    EXPECT_LE(lu.report.linearError, Real(1e-10));
  }

#ifdef RODIN_USE_MUMPS
  TEST(Rodin_Adaptation_WNGIRSolver, InactiveHingesSkipZeroCorrections)
  {
    const auto state = solveTranslatedLine(Real(1), 0, false, 0, false, 1000,
      true, WNGIRParameters::DirectSolver::MUMPS, false, Real(1e-3), Real(0));
    EXPECT_GT(state.report.inactiveHingeSkips, 0u);
    EXPECT_EQ(state.report.primalBarrierIterations, 0u);
    EXPECT_EQ(state.report.tPrimalBarrierAssembly, Real(0));
    EXPECT_EQ(state.report.tPrimalBarrierSolve, Real(0));
    EXPECT_TRUE(state.report.primalBarrierConverged);
    EXPECT_LT(state.report.geometricSup, Real(0.05));
  }

  TEST(Rodin_Adaptation_WNGIRSolver, RejectsDisabledInnerResidualTest)
  {
    EXPECT_THROW(solveTranslatedLine(Real(1), 0, false, 0, false, 1000,
      true, WNGIRParameters::DirectSolver::MUMPS, false, Real(0), Real(0)), Alert::Exception);
  }

  TEST(Rodin_Adaptation_WNGIRSolver, MUMPSAugmentedSolvePreservesSimilarityGauge)
  {
    const auto solve = [](WNGIRParameters::DirectSolver directSolver) {
      return solveTranslatedLine(
        Real(1), Real(0), false, Real(0), false, 1000, true, directSolver);
    };
    const auto lu = solve(WNGIRParameters::DirectSolver::SparseLU);
    const auto mumps = solve(WNGIRParameters::DirectSolver::MUMPS);
    EXPECT_EQ(lu.report.iterations, mumps.report.iterations);
    EXPECT_STREQ(lu.report.exitReason, mumps.report.exitReason);
    EXPECT_NEAR((lu.displacement - mumps.displacement).norm(), Real(0), Real(1e-8));
    EXPECT_EQ(mumps.report.linearIterations, 0u);
    EXPECT_GT(mumps.report.linearSolveCount, 0u);
    EXPECT_LE(mumps.report.linearError, Real(1e-10));
    EXPECT_GT(mumps.report.directFactorizations, 0u);
    EXPECT_LT(mumps.report.directAnalyses, mumps.report.directFactorizations);
  }
#endif

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
  TEST(Rodin_Adaptation_WNGIRSolver, UnresolvedTranslationHasNoInertiaGate)
  {
    const auto state = solveTranslatedLine(Real(1), 0, false, 0, false, 1000, true,
      WNGIRParameters::DirectSolver::SparseLU, true);
    EXPECT_STRNE(state.report.exitReason, "metric-inertia-unresolved");
    EXPECT_GT(state.report.linearSolveCount, 0u);
    EXPECT_GT(state.report.minJ, Real(0.01));
    EXPECT_LT(state.report.maxQRel, Real(10));
    EXPECT_GE(state.report.unresolvedSimilarityModes, 1u);
    EXPECT_TRUE(state.report.geometricTargetReached);
    // A plane does not distinguish translation from uniform dilation about it.
    // The coefficient-minimum solution may combine them, but must remain a similarity.
    EXPECT_NEAR(state.report.minJ, state.report.maxJ, Real(1e-10));
    EXPECT_NEAR(state.report.maxQRel, Real(1), Real(1e-10));
  }

#ifdef RODIN_USE_MUMPS
  TEST(Rodin_Adaptation_WNGIRSolver, SimilarityGaugePreservesPlaneTranslation3D)
  {
    constexpr size_t n = 3;
    constexpr Real h = Real(0.5);
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(h);
    for (size_t from = 1; from <= 3; ++from)
      for (size_t to = 0; to <= 3; ++to)
        if (from != to)
          mesh.getConnectivity().compute(from, to);
    std::vector<Index> facets;
    for (auto face = mesh.getFace(); face; ++face)
    {
      bool onInterface = true;
      for (const Index vertex : face->getVertices())
        onInterface &=
          std::abs(mesh.getVertexCoordinates(vertex)(0) - Real(0.5)) < Real(1e-12);
      if (onInterface)
      {
        facets.push_back(face->getIndex());
        mesh.setAttribute({2, face->getIndex()}, Interface);
      }
    }
    ASSERT_FALSE(facets.empty());
    P1<Math::SpatialVector<Real>, LocalMesh> fes(mesh, 3);
    TrialFunction trial(fes);
    TestFunction test(fes);
    WNGIR solver(trial, test);
    WNGIRParameters p;
    p.h = h;
    p.directSolver = WNGIRParameters::DirectSolver::MUMPS;
    p.directSolverThreads = 1;
    p.maxIterations = 12;
    p.hasInterfaceAttribute = true;
    p.interfaceAttribute = Interface;
    p.geometricSupTolerance = Real(1e-6);
    p.cgRelativeTolerance = Real(1e-10);
    solver.setParameters(p);
    RealFunction phi([](const Point& point) { return point.x() - Real(0.55); });
    Math::Vector<Real> normal(3);
    normal << Real(1), Real(0), Real(0);
    const VectorFunction gradient(normal);
    const auto report = solver.solve(mesh, facets, phi, gradient);
    EXPECT_TRUE(report.geometricTargetReached);
    EXPECT_EQ(report.unresolvedSimilarityModes, 4u);
    EXPECT_LE(report.linearError, p.cgRelativeTolerance);
    EXPECT_NEAR(report.minJ, Real(1), Real(1e-10));
    EXPECT_NEAR(report.maxJ, Real(1), Real(1e-10));
    EXPECT_NEAR(report.maxQRel, Real(1), Real(1e-10));
    GridFunction expected(fes);
    normal *= Real(0.05);
    expected = VectorFunction(normal);
    EXPECT_LT((trial.getSolution().getData() - expected.getData()).norm(),
      Real(1e-6) * std::sqrt(Real(fes.getSize())));
  }
#endif

  TEST(Rodin_Adaptation_WNGIRSolver, CommonMetricScalePreservesTheHingeModelP1P2)
  {
    const auto check = []<size_t Order>() {
      const auto solve = [](Real scale) {
        return solveTranslatedLine<Order>(1, 0, false, 0, false, 1000, true,
          WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-6), Real(1000),
          [=](WNGIRParameters& p) {
            p.kappaF = scale;
            p.kappaS = scale * Real(0.1);
            p.kappaD = scale;
            p.qualityGuard = Real(0.9);
            p.maxIterations = 3;
          });
      };
      const auto reference = solve(1);
      ASSERT_GT(reference.report.iterations, 0u);
      ASSERT_GT(reference.report.primalBarrierIterations, 0u);
      for (const Real scale : {Real(1e-4), Real(1e3)})
      {
        const auto scaled = solve(scale);
        EXPECT_EQ(reference.report.iterations, scaled.report.iterations);
        EXPECT_EQ(reference.report.primalBarrierIterations,
          scaled.report.primalBarrierIterations);
        EXPECT_STREQ(reference.report.exitReason, scaled.report.exitReason);
        EXPECT_NEAR((reference.displacement - scaled.displacement).norm(), 0, 1e-8);
        EXPECT_NEAR(reference.report.primalBarrierCoefficient,
          scaled.report.primalBarrierCoefficient, 1e-9);
      }
    };
    check.template operator()<1>();
    check.template operator()<2>();
  }

  TEST(Rodin_Adaptation_WNGIRSolver, InnerResidualUsesTheFixedForceScale)
  {
    testing::internal::CaptureStdout();
    const auto state = solveTranslatedLine(1, 0, true, 0, false, 1000, true,
      WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), Real(1000),
      [](WNGIRParameters& p) {
        p.qualityGuard = Real(0.95);
        p.maxIterations = 1;
      });
    const auto output = testing::internal::GetCapturedStdout();
    EXPECT_GT(state.report.primalBarrierIterations, 0u);
    std::istringstream lines(output);
    std::string line;
    size_t checked = 0;
    while (std::getline(lines, line))
      if (line.find("barrier residual:") != std::string::npos)
      {
        const Real force = std::stod(line.substr(line.find("force_norm=") + 11));
        const Real threshold = std::stod(line.substr(line.find("tolerance=") + 10));
        EXPECT_NEAR(threshold, Real(1e-12) + Real(1e-3) * force, Real(1e-9) * force);
        ++checked;
      }
    EXPECT_GT(checked, 1u);
  }

  TEST(Rodin_Adaptation_WNGIRSolver, RejectsInvalidCanonicalWeights)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1<Math::SpatialVector<Real>, LocalMesh> fes(mesh, 2);
    TrialFunction trial(fes);
    TestFunction test(fes);
    WNGIR solver(trial, test);
    for (const Real invalid : {Real(-1), std::numeric_limits<Real>::quiet_NaN()})
    {
      WNGIRParameters p;
      p.kappaS = invalid;
      EXPECT_THROW(solver.setParameters(p), Alert::Exception);
    }
    for (const Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::infinity()})
    {
      WNGIRParameters p;
      p.kappaD = invalid;
      EXPECT_THROW(solver.setParameters(p), Alert::Exception);
      p = WNGIRParameters{};
      p.kappaF = invalid;
      EXPECT_THROW(solver.setParameters(p), Alert::Exception);
    }
    for (const Real invalid : {Real(0), Real(1), std::numeric_limits<Real>::quiet_NaN()})
    {
      WNGIRParameters p;
      p.qualityGuard = invalid;
      EXPECT_THROW(solver.setParameters(p), Alert::Exception);
    }
  }

  TEST(Rodin_Adaptation_WNGIRSolver, StrictCGCapIsNotBypassed)
  {
    const auto state = solveTranslatedLine(
      Real(1), 0, false, 0, false, 1, true, WNGIRParameters::DirectSolver::CG);
    EXPECT_STREQ(state.report.exitReason, "solve-predictor-failed");
    EXPECT_EQ(state.report.linearSolveCount, 1u);
    EXPECT_LE(state.report.maxLinearIterations, 1u);
    EXPECT_EQ(state.report.iterations, 0u);
  }

  TEST(Rodin_Adaptation_WNGIRSolver, RejectsMismatchedInterface)
  {
    EXPECT_THROW(solveTranslatedLine(Real(1), 0, false, 0, false, 1000,
      true, WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), Real(90),
      [](WNGIRParameters& p) { p.interfaceAttribute = 999; }), Alert::Exception);
  }

  TEST(Rodin_Adaptation_WNGIRSolver, SmallStepsRequireTheConfiguredPersistence)
  {
    const auto state = solveTranslatedLine(Real(1), 0, false, Real(1e-12), false, 1000,
      true, WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), Real(0),
      [](WNGIRParameters& p) {
        p.kappaF = Real(100);
        p.stepTol = Real(1);
        p.stagnationIterations = 3;
      });
    EXPECT_EQ(state.report.iterations, 3u);
    EXPECT_STREQ(state.report.exitReason, "best-effort-step-stagnation");
    EXPECT_FALSE(state.report.geometricTargetReached);
  }

  TEST(Rodin_Adaptation_WNGIRSolver, InterfaceVerticesContributeToSupremum)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.scale(Real(0.5));
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(1, 2);
    std::vector<Index> facets;
    for (auto face = mesh.getFace(); face; ++face)
    {
      bool marked = true;
      for (const Index vertex : face->getVertices())
        marked &= std::abs(mesh.getVertexCoordinates(vertex)(0) - Real(0.5)) < Real(1e-12);
      if (marked)
      {
        facets.push_back(face->getIndex());
        mesh.setAttribute({1, face->getIndex()}, Interface);
      }
    }
    P1<Math::SpatialVector<Real>, LocalMesh> fes(mesh, 2);
    TrialFunction trial(fes);
    TestFunction test(fes);
    WNGIR solver(trial, test);
    WNGIRParameters p;
    p.h = Real(0.5);
    p.hasInterfaceAttribute = true;
    p.interfaceAttribute = Interface;
    p.geometricSupTolerance = Real(0.1);
    solver.setParameters(p);
    RealFunction phi([](const Point& point) {
      const Real y = Real(2) * point.y() - Real(1);
      return point.x() - Real(0.5) + Real(0.05) * y * y;
    });
    AnalyticVectorFunction grad([](const Point& point) {
      return Math::SpatialVector<Real>{Real(1), Real(0.2) * (Real(2) * point.y() - Real(1))};
    }, 2);
    const auto report = solver.solve(mesh, facets, phi, grad);
    EXPECT_EQ(report.iterations, 0u);
    EXPECT_NEAR(report.geometricSup, Real(0.05) / std::sqrt(Real(1.04)), Real(1e-12));
  }

  TEST(Rodin_Adaptation_WNGIRSolver, InvalidGeometryCannotReachAutomaticTarget)
  {
    const auto state = solveTranslatedLine(Real(1), 0, false, 0, true, 1000,
      true, WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), Real(90),
      [](WNGIRParameters& p) { p.geometricSupTolerance = 0; });
    EXPECT_STREQ(state.report.exitReason, "geometric-validation-failed");
    EXPECT_FALSE(state.report.geometricTargetReached);
    EXPECT_TRUE(std::isinf(state.report.geometricSup));
  }

  TEST(Rodin_Adaptation_WNGIRSolver, ResidualCertifiesInactiveAndActiveInnerSolves)
  {
    for (const Real mu : {Real(0), Real(90)})
    {
      const auto state = solveTranslatedLine(Real(1), 0, false, 0, false, 1000,
        true, WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), mu,
        [](WNGIRParameters& p) { p.maxIterations = 1; });
      ASSERT_TRUE(state.report.primalBarrierConverged);
      EXPECT_LE(state.report.primalBarrierResidual, state.report.primalBarrierResidualTolerance);
      EXPECT_LE(state.report.maxPrimalBarrierIterations, 15u);
    }
  }

  TEST(Rodin_Adaptation_WNGIRSolver, P2AcceptedStepUsesPhysicalField)
  {
    const auto state = solveTranslatedLine<2>(Real(1), 0, false, Real(1e-12), false,
      1000, true, WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), Real(0),
      [](WNGIRParameters& p) {
        p.maxIterations = 1;
        p.cgRelativeTolerance = Real(1e-8);
        p.kappaS = 0;
        p.kappaD = Real(1e-4);
      });
    ASSERT_EQ(state.report.iterations, 1u);
    EXPECT_GT(state.report.acceptedStep, Real(0));
    EXPECT_NEAR(state.report.acceptedStep, state.physicalNorm, Real(1e-12));
  }

  TEST(Rodin_Adaptation_WNGIRSolver, P2PhysicalNormDetectsAnInteriorMaximum)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.scale(Real(0.5));
    for (std::size_t from = 1; from <= 2; ++from)
      for (std::size_t to = 0; to <= 2; ++to)
        if (from != to)
          mesh.getConnectivity().compute(from, to);
    H1 fes(std::integral_constant<std::size_t, 2>{}, mesh, 2);
    GridFunction field(fes);
    const auto cell = mesh.getCell();
    const auto& qf = QF::PolytopeQuadratureFormula::get(8, cell->getGeometry());
    const Real peak = cell->getQuadrature(qf).getPoint(0).x();
    field = VectorFunction(std::size_t(2), [peak](const Point& point) {
      const Real x = point.x() - peak;
      return Math::SpatialVector<Real>{Real(1) - x * x, Real(0)};
    });
    std::vector<Index> cells;
    for (auto current = mesh.getCell(); current; ++current)
      cells.push_back(current->getIndex());
    const Real norm = Detail::wngirPhysicalDisplacementNorm(mesh, fes, cells, field, 8);
    EXPECT_NEAR(norm, Real(1), Real(1e-12));
    EXPECT_LT(field.getData().cwiseAbs().maxCoeff(), norm - Real(1e-8));
  }

  TEST(Rodin_Adaptation_WNGIRSolver, RejectsInvalidHingeWeightsAndControls)
  {
    for (const auto member : {&WNGIRParameters::kappaJ, &WNGIRParameters::kappaQ,
           &WNGIRParameters::primalBarrierRelativeTolerance, &WNGIRParameters::energyStagTol})
      for (const Real invalid : {Real(-1), std::numeric_limits<Real>::infinity(),
             std::numeric_limits<Real>::quiet_NaN()})
        EXPECT_THROW(solveTranslatedLine(Real(1), 0, false, 0, false, 1000, true,
          WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), Real(90),
          [=](WNGIRParameters& p) { p.*member = invalid; }), Alert::Exception);
    EXPECT_THROW(solveTranslatedLine(Real(1), 0, false, 0, false, 1000, true,
      WNGIRParameters::DirectSolver::SparseLU, false, Real(1e-3), Real(90),
      [](WNGIRParameters& p) { p.h = 0; }), Alert::Exception);
  }

}
