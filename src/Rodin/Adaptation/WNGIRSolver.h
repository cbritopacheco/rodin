/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRSOLVER_H
#define RODIN_ADAPTATION_WNGIRSOLVER_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <Eigen/Eigenvalues>

#include "Rodin/Solver/CG.h"
#include "Rodin/Variational/Integrator.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/Trace.h"
#include "Rodin/Variational/IdentityMatrix.h"

#include "Rodin/Location.h"

#include "CellDeformation.h"
#include "WNGIRLoss.h"
#include "WNGIRParameters.h"
#include "WNGIRPrimalBarrierForce.h"
#include "WNGIRPrimalBarrierMetric.h"
#include "WNGIRReport.h"
#include "WNGIRNormalOffsetCoefficient.h"
#include "WNGIRNormalOffsetForceCoefficient.h"
#include "WNGIRObservationCoefficient.h"
#include "WNGIRSurfaceForceCoefficient.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Backend-independent WNGIR mesh-fitting solver.
   *
   * Fits the interface skeleton of a mesh to the zero level set of @f$\phi@f$
   * by repeatedly solving a regularized first-derivative metric problem and
   * line-searching the resulting displacement on the true geometry.
   *
   * Constructed from backend-matched trial and test functions, the symmetric
   * metric path uses the ordinary Rodin form-language deduction path for
   * @ref Rodin::Variational::Problem, @ref Rodin::Variational::BilinearForm,
   * and @ref Rodin::Solver::CG.
   */
  template <class TrialFunctionType, class TestFunctionType>
  class WNGIR
  {
      using Displacement = std::remove_reference_t<
        decltype(std::declval<TrialFunctionType&>().getSolution())>;
      using FESType = std::remove_reference_t<
        decltype(std::declval<Displacement&>().getFiniteElementSpace())>;
      using ProblemType = std::decay_t<decltype(Variational::Problem(
        std::declval<TrialFunctionType&>(), std::declval<TestFunctionType&>()))>;
      using LinearSystemType = typename ProblemType::LinearSystemType;
      using StepSolverType = Solver::CG<LinearSystemType>;
      using BilinearFormType = std::decay_t<decltype(Variational::BilinearForm(
        std::declval<TrialFunctionType&>(), std::declval<TestFunctionType&>()))>;

      using SpatialVec = Math::SpatialVector<Real>;
      using SpatialMat = Math::SpatialMatrix<Real>;

    public:
      /// @brief Mesh-validity summary of a displacement.
      struct AdmissibilityState
      {
        /// @brief Smallest Jacobian over all cell quadrature points.
          Real minJ = std::numeric_limits<Real>::infinity();

        /// @brief Largest Jacobian over all cell quadrature points.
          Real maxJ = -std::numeric_limits<Real>::infinity();

        /// @brief Largest relative distortion over the non-inverted points.
          Real maxQ = 0;

        /// @brief Number of quadrature points at or below the Jacobian floor.
          std::size_t inadmissibleCount = 0;
      };

      /// @brief Interface-fit summary of a displacement.
      struct SurfaceState
      {
        /// @brief Robust energy of the interface residual.
          Real energy = 0;

        /// @brief Measure carrying a non-negligible robust weight.
          Real activeLen = 0;

        /// @brief Measure of the whole interface.
          Real totalLen = 0;

        /// @brief Root-mean-square residual over the active part.
          Real activeRMS = 0;

        /// @brief Supremum residual over the active part.
          Real activeSup = 0;
      };

      /// @brief Observation coercivity restricted to rigid displacements.
      struct RigidModeState
      {
        /// @brief Smallest generalized observation eigenvalue.
          Real minimum = std::numeric_limits<Real>::quiet_NaN();

        /// @brief Ratio of the smallest to largest generalized eigenvalue.
          Real ratio = std::numeric_limits<Real>::quiet_NaN();

        /// @brief Dimension of the rigid-motion space tested.
          std::size_t dimension = 0;
      };

      /**
       * @brief Kink of the discrete interface across facet boundaries.
       *
       * The interface is only piecewise smooth: adjacent facets meet at an
       * angle that does not vanish under refinement, and no displacement can
       * fit a feature finer than that kink. The statistics therefore floor the
       * residual tolerances.
       */
      struct NormalJump
      {
        /// @brief Measure-weighted RMS angle between adjacent facet normals.
          Real rms = 0;

        /// @brief Largest angle between adjacent facet normals.
          Real max = 0;
      };

      /// @brief Constructs the WNGIR solver from trial and test functions.
      WNGIR(TrialFunctionType& du, TestFunctionType& v)
        : m_u(&du.getSolution()),
          m_duStep(du.getFiniteElementSpace()),
          m_vStep(v.getFiniteElementSpace()),
          m_stepProblem(m_duStep, m_vStep),
          m_stepSolver(m_stepProblem),
          m_bulkForm(m_duStep, m_vStep)
      {}

      /// @brief Sets WNGIR runtime parameters.
      WNGIR& setParameters(const WNGIRParameters& parameters)
      {
        m_parameters = parameters;
        const std::size_t defaultMaxIterations = std::min<std::size_t>(2000,
          std::max<std::size_t>(100, 2 * m_duStep.getFiniteElementSpace().getSize()));
        const std::size_t maxIterations = m_parameters.cgMaxIterations > 0
          ? m_parameters.cgMaxIterations
          : defaultMaxIterations;
        if constexpr (requires(StepSolverType& solver, Real tolerance) {
                        solver.setTolerance(tolerance);
                      })
          m_stepSolver.setTolerance(m_parameters.cgRelativeTolerance);
        if constexpr (requires(StepSolverType& solver, std::size_t iterations) {
                        solver.setMaxIterations(iterations);
                      })
          m_stepSolver.setMaxIterations(maxIterations);
        if constexpr (requires(
                        StepSolverType& solver, Real tolerance, std::size_t iterations) {
                        solver.setTolerances(tolerance, tolerance, tolerance, iterations);
                      })
        {
          m_stepSolver.setTolerances(m_parameters.cgRelativeTolerance,
            StepSolverType::DEFAULT_ABSTOL, StepSolverType::DEFAULT_DTOL, maxIterations);
        }
        m_bulkFormAssembled = false;
        return *this;
      }

      /// @brief Returns the current WNGIR parameters.
      const WNGIRParameters& getParameters() const
      {
        return m_parameters;
      }

      /// @brief Returns diagnostics from the most recent solve.
      const WNGIRReport& getReport() const
      {
        return m_report;
      }

      /**
       * @brief Solves the WNGIR fitting problem on marked interface facets.
       *
       * The supplied sensitivity @p grad must equal the derivative of @p phi
       * at the moved quadrature points for the assembled force to be the exact
       * first variation of the line-search energy. An independently supplied
       * sensitivity is supported, but then defines a pseudo-gradient.
       */
      template <class Mesh, class PhiDerived, class GradDerived>
      WNGIRReport solve(const Mesh& mesh,
        const std::vector<Rodin::Index>& interfaceFacets,
        const Variational::RealFunctionBase<PhiDerived>& phi,
        const Variational::VectorFunctionBase<Real, GradDerived>& grad)
      {
        using Rodin::Index;
        const WNGIRParameters& p = m_parameters;
        Displacement& u = *m_u;
        const auto& fes = u.getFiniteElementSpace();
        const std::size_t meshDim = mesh.getDimension();
        const Real d = static_cast<Real>(meshDim);

        WNGIRReport rep;
        const Real h = p.h;
        const bool usePrimalBarrier = p.kappaJ > Real(0) || p.kappaQ > Real(0);
        const bool checkGeometry = p.admissibilityChecks || usePrimalBarrier;
        const Real acceptedJacobianFloor =
          usePrimalBarrier ? std::max(p.jLineSearchRatio, p.jSafe) : p.jLineSearchRatio;
        const Real bulkDeviatoricCoefficient = h * p.kappaBulk;
        const Real bulkDivergenceCoefficient =
          bulkDeviatoricCoefficient * p.rDiv;
        const Real stepTol = p.stepTol > Real(0) ? p.stepTol : Real(1e-4) * h;
        using Clock = std::chrono::steady_clock;
        auto secondsSince = [](Clock::time_point t0) -> Real {
          return std::chrono::duration<Real>(Clock::now() - t0).count();
        };
        auto setupTic = Clock::now();

        // ============================================================
        // Per-frame geometry tabulation (only u changes per iteration).
        // Field values come from GridFunction::getValue (cached), so the
        // tables store quadrature geometry only.
        // ============================================================
        const auto normalJump = getNormalJump(mesh, fes, interfaceFacets, meshDim);
        rep.normalJumpRMS = normalJump.rms;
        rep.normalJumpMax = normalJump.max;
        std::vector<Index> validationCells;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const Index cellIndex = cell->getIndex();
          if constexpr (requires { mesh.getShard().isOwned(meshDim, cellIndex); })
          {
            if (!mesh.getShard().isOwned(meshDim, cellIndex))
              continue;
          }
          validationCells.push_back(cellIndex);
        }
        const Real floorRMS = meshDim == 3 ? p.rmsFloor3D : p.rmsFloor2D;
        const Real floorSup = meshDim == 3 ? p.supFloor3D : p.supFloor2D;
        Real activeRMSOverHTol = p.activeRMSOverHTol;
        Real activeSupOverHTol = p.activeSupOverHTol;
        if (p.geometryAwareTolerances)
        {
          activeRMSOverHTol = std::max(std::max(activeRMSOverHTol, floorRMS),
            p.rmsNormalJumpFactor * rep.normalJumpRMS);
          activeSupOverHTol = std::max(std::max(activeSupOverHTol, floorSup),
            p.supNormalJumpFactor * rep.normalJumpMax);
        }
        rep.effectiveRMSOverHTol = activeRMSOverHTol;
        rep.effectiveSupOverHTol = activeSupOverHTol;

        if (!p.hasInterfaceAttribute)
        {
          rep.exitReason = "missing-interface-attribute";
          m_report = rep;
          return rep;
        }

        // One bounding-volume locator per solve: the background mesh is fixed
        // for the frame, and the index builds lazily on first query.
        const Location::AABB<Mesh> locator(mesh);
        if (interfaceFacets.empty())
        {
          rep.exitReason = "empty-interface";
          m_report = rep;
          return rep;
        }
        const auto robustScale = getRobustScale(mesh, fes, phi, grad, interfaceFacets, h);
        if (!(robustScale.gradientScale > Real(0)) ||
          !std::isfinite(robustScale.gradientScale))
        {
          rep.exitReason = "degenerate-target-gradient";
          m_report = rep;
          return rep;
        }
        const Real sigma = robustScale.sigma;
        const Real sigma2 = sigma * sigma;
        const Real levelSetMeshScale = h * robustScale.gradientScale;
        const Real dataNormalization =
          Real(1) / (robustScale.gradientScale * robustScale.gradientScale);
        const Real activeRMSTol =
          p.activeRMSTol > Real(0) ? p.activeRMSTol : Real(4) * h * levelSetMeshScale;
        const Real activeSupTol =
          p.activeSupTol > Real(0) ? p.activeSupTol : Real(10) * h * levelSetMeshScale;
        rep.effectiveRMSTol = activeRMSTol;
        rep.effectiveSupTol = activeSupTol;
        const WNGIRLoss loss(sigma);
        rep.sigma = sigma;
        rep.levelSetGradientScale = robustScale.gradientScale;
        const Real domainMeasure =
          usePrimalBarrier ? getDomainMeasure(mesh, fes, validationCells) : Real(0);

        // ============================================================
        // Per-iteration field evaluations through the GridFunction.
        // ============================================================
        using FastAdm = AdmissibilityState;
        auto fastAdmissibility = [&](const Displacement& gf) {
          return getAdmissibilityState(mesh, fes, validationCells, gf, meshDim);
        };
        auto surfaceState = [&](const Displacement& gf) {
          return getSurfaceState(mesh, fes, gf, phi, interfaceFacets, loss,
            dataNormalization, meshDim, locator);
        };
        auto recordSurfaceState = [&rep](const SurfaceState& state) {
          rep.energy = state.energy;
          rep.activeRMS = state.activeRMS;
          rep.activeSup = state.activeSup;
          rep.activeMeasure = state.activeLen;
          rep.interfaceMeasure = state.totalLen;
          rep.activeFraction =
            state.totalLen > Real(0) ? state.activeLen / state.totalLen : Real(0);
        };

        rep.tSetup = secondsSince(setupTic);

        // ============================================================
        // Constant bulk metric, assembled once per frame.
        // ============================================================
        if (!m_bulkFormAssembled)
        {
          auto bulkTic = Clock::now();
          const auto epsTrial = Real(0.5) *
            (Variational::Jacobian(m_duStep) + Variational::Jacobian(m_duStep).T());
          const auto epsTest = Real(0.5) *
            (Variational::Jacobian(m_vStep) + Variational::Jacobian(m_vStep).T());
          const auto divTrial = Variational::Trace(epsTrial);
          const auto divTest = Variational::Trace(epsTest);
          const auto devTrial =
            epsTrial - (Real(1) / d) * divTrial * Variational::IdentityMatrix(meshDim);
          const auto devTest =
            epsTest - (Real(1) / d) * divTest * Variational::IdentityMatrix(meshDim);
          m_bulkForm =
            Variational::Integral(bulkDeviatoricCoefficient * devTrial, devTest) +
            Variational::Integral(bulkDivergenceCoefficient * divTrial, divTest);
          m_bulkForm.assemble();
          rep.tBulk = secondsSince(bulkTic);
          m_bulkFormAssembled = true;
        }

        Displacement vK(fes);
        Displacement predictor(fes);
        Displacement scratch(fes);
        Displacement previousU(fes);
        Displacement uTrial(fes);

        // The interface coefficients are not polynomial, so the quadrature
        // order is given to the integrator explicitly, per face.
        const Variational::Integrator::OrderType surfaceOrder =
          [&](const Geometry::Polytope& face) -> std::size_t {
          if (p.quadratureOrder > 0)
            return p.quadratureOrder;
          const auto& fe = fes.getFiniteElement(face.getDimension(), face.getIndex());
          return std::max<std::size_t>(2, 2 * fe.getOrder());
        };

        if (p.kappaInit > Real(0))
        {
          auto tic = Clock::now();
          Detail::WNGIRNormalOffsetCoefficient initCoeff(phi, grad, p, sigma2, meshDim);
          auto initMetric =
            Variational::FaceIntegral(Variational::Dot(initCoeff * m_duStep, m_vStep));
          initMetric.setOrder(surfaceOrder);
          initMetric.over(p.interfaceAttribute);
          Detail::WNGIRNormalOffsetForceCoefficient initForceCoeff(
            phi, grad, p, sigma2, meshDim);
          auto initForce = Variational::FaceIntegral(initForceCoeff, m_vStep);
          initForce.setOrder(surfaceOrder);
          initForce.over(p.interfaceAttribute);
          typename ProblemType::ProblemBodyType body(m_bulkForm);
          body = body + initMetric - initForce;
          Math::Vector<Real> zero(meshDim);
          zero.setZero();
          if (!p.dirichletAttributes.empty())
            body = body +
              Variational::DirichletBC(m_duStep, Variational::VectorFunction(zero))
                .on(p.dirichletAttributes);
          m_stepProblem = body;
          m_stepProblem.assemble();
          rep.tAssembly += secondsSince(tic);

          tic = Clock::now();
          std::size_t linearIterations = 0;
          Real linearError = 0;
          const bool solveOk = solveStep(vK, linearIterations, linearError);
          rep.tSolve += secondsSince(tic);
          rep.linearIterations += linearIterations;
          rep.linearError = linearError;
          if (!solveOk)
          {
            rep.exitReason = "initial-guess-solve-failed";
            m_report = rep;
            return rep;
          }

          tic = Clock::now();
          const bool useLineSearch = checkGeometry || p.energyLineSearch;
          Real eBase = 0;
          Real alpha = Real(1);
          bool accepted = false;
          FastAdm adm{};
          Real eTrial = std::numeric_limits<Real>::infinity();
          SurfaceState trialSurface{};
          bool trialSurfaceEvaluated = false;
          if (!useLineSearch)
          {
            u += vK;
            accepted = true;
          }
          else
          {
            previousU = u;
            if (p.energyLineSearch)
              eBase = surfaceState(previousU).energy;
            while (alpha >= p.alphaMin)
            {
              uTrial = vK;
              uTrial *= alpha;
              uTrial += previousU;
              bool jOK = true;
              bool qOK = true;
              if (checkGeometry)
              {
                adm = fastAdmissibility(uTrial);
                jOK = adm.inadmissibleCount == 0 && adm.minJ > acceptedJacobianFloor;
                qOK = adm.maxQ < p.qMax;
              }
              bool eOK = true;
              if (jOK && qOK && p.energyLineSearch)
              {
                trialSurface = surfaceState(uTrial);
                trialSurfaceEvaluated = true;
                eTrial = trialSurface.energy;
                eOK = std::isfinite(eTrial) && eTrial <= eBase;
              }
              if (jOK && qOK && eOK)
              {
                u = uTrial;
                accepted = true;
                break;
              }
              if (!jOK)
                ++rep.jacobianRejections;
              if (!qOK)
                ++rep.distortionRejections;
              if (jOK && qOK && !eOK)
                ++rep.energyRejections;
              alpha *= Real(0.5);
              ++rep.backtracks;
            }
          }
          rep.tLineSearch += secondsSince(tic);
          if (p.trace)
          {
            const auto surf =
              accepted && trialSurfaceEvaluated ? trialSurface : surfaceState(u);
            std::cout << "      wngir init=normal-offset"
                      << "  accepted=" << (accepted ? 1 : 0) << "  alpha=" << alpha
                      << "  actRMS=" << std::scientific << surf.activeRMS
                      << "  actRMS/(hG)="
                      << (levelSetMeshScale > Real(0) ? surf.activeRMS / levelSetMeshScale
                                                      : Real(0))
                      << "  min_j="
                      << (checkGeometry && accepted
                             ? adm.minJ
                             : std::numeric_limits<Real>::quiet_NaN())
                      << "  max_Q="
                      << (checkGeometry && accepted
                             ? adm.maxQ
                             : std::numeric_limits<Real>::quiet_NaN())
                      << '\n';
          }
        }

        FastAdm initialAdm{};
        if (checkGeometry)
        {
          initialAdm = fastAdmissibility(u);
          rep.minJ = initialAdm.minJ;
          rep.maxJ = initialAdm.maxJ;
          rep.maxQRel = initialAdm.maxQ;
        }
        if (usePrimalBarrier)
        {
          if (initialAdm.inadmissibleCount > 0 || initialAdm.minJ <= p.jSafe ||
            initialAdm.maxQ >= p.qMax)
          {
            rep.exitReason = "initial-state-not-strictly-feasible";
            m_report = rep;
            return rep;
          }
        }

        SurfaceState currentSurface = surfaceState(u);
        recordSurfaceState(currentSurface);
        const RigidModeState initialRigid = getRigidModeState(mesh, fes, u, phi, grad,
          interfaceFacets, sigma2, dataNormalization, meshDim, locator);
        rep.rigidModeCoercivity = initialRigid.minimum;
        rep.rigidModeCoercivityRatio = initialRigid.ratio;
        rep.rigidModeDimension = initialRigid.dimension;
        if (!(currentSurface.activeLen > Real(0)))
        {
          rep.exitReason = "observation-degenerate-active-set";
          m_report = rep;
          return rep;
        }
        if (initialRigid.dimension > 0 &&
          (!(initialRigid.minimum > Real(0)) || !std::isfinite(initialRigid.minimum)))
        {
          rep.exitReason = "observation-degenerate-rigid-modes";
          m_report = rep;
          return rep;
        }

        Real ePrev = currentSurface.energy;

        // ============================================================
        // Nonlinear iteration.
        // ============================================================
        for (; rep.iterations < p.maxIterations; ++rep.iterations)
        {
          auto tic = Clock::now();
          Detail::WNGIRObservationCoefficient obsCoeff(
            phi, grad, u, locator, p, sigma2, dataNormalization, meshDim);
          auto obsMetric =
            Variational::FaceIntegral(Variational::Dot(obsCoeff * m_duStep, m_vStep));
          obsMetric.setOrder(surfaceOrder);
          obsMetric.over(p.interfaceAttribute);
          Detail::WNGIRSurfaceForceCoefficient forceCoeff(
            phi, grad, u, locator, sigma2, dataNormalization, meshDim);
          auto surfaceForce = Variational::FaceIntegral(forceCoeff, m_vStep);
          surfaceForce.setOrder(surfaceOrder);
          surfaceForce.over(p.interfaceAttribute);
          std::size_t linearIterations = 0;
          Real linearError = std::numeric_limits<Real>::infinity();
          Math::Vector<Real> zero(meshDim);
          zero.setZero();
          bool solveOk = true;
          Real predictorAction = Real(0);
          typename ProblemType::ProblemBodyType predictorBody(m_bulkForm);
          predictorBody = predictorBody + obsMetric - surfaceForce;
          if (!p.dirichletAttributes.empty())
            predictorBody = predictorBody +
              Variational::DirichletBC(m_duStep, Variational::VectorFunction(zero))
                .on(p.dirichletAttributes);
          m_stepProblem = predictorBody;
          m_stepProblem.assemble();
          rep.tAssembly += secondsSince(tic);

          tic = Clock::now();
          std::size_t predictorIterations = 0;
          Real predictorError = std::numeric_limits<Real>::infinity();
          solveOk = solveStep(vK, predictorIterations, predictorError);
          linearIterations += predictorIterations;
          linearError = predictorError;
          rep.tSolve += secondsSince(tic);
          if (!solveOk)
          {
            rep.exitReason = "solve-predictor-failed";
            break;
          }

          predictor = vK;
          predictorAction = std::max(Real(0),
            getSurfaceForceAction(mesh, fes, u, predictor, phi, grad, interfaceFacets,
              sigma2, dataNormalization, meshDim, locator));
          rep.stationarityNorm = std::sqrt(predictorAction);
          if (p.stationarityTolerance > Real(0) &&
            rep.stationarityNorm <= p.stationarityTolerance)
          {
            rep.energy = ePrev;
            rep.exitReason = "stationary-converged";
            break;
          }

          if (usePrimalBarrier)
          {
            const Real modelDecrease = Real(0.5) * predictorAction;
            const Real barrierCoefficient = domainMeasure > Real(0)
              ? p.muHat * modelDecrease / domainMeasure
              : Real(0);
            rep.primalBarrierCoefficient = barrierCoefficient;

            uTrial *= Real(0);
            const Real predictorAlpha =
              getBarrierStepScale(mesh, fes, validationCells, u, uTrial, vK, meshDim);
            if (!(predictorAlpha > Real(0)))
            {
              solveOk = false;
            }
            else
            {
              // Start the barrier corrections from the feasible predictor,
              // not from the origin of the affine increment problem.
              vK *= predictorAlpha;
              const std::size_t innerIterations =
                std::max<std::size_t>(1, p.primalBarrierIterations);
              const Real predictorNorm =
                std::max(std::abs(predictor.max()), std::abs(predictor.min()));
              bool innerConverged = false;
              rep.lastPrimalBarrierIterations = 0;
              tic = Clock::now();
              for (std::size_t inner = 0; inner < innerIterations; ++inner)
              {
                Detail::WNGIRPrimalBarrierMetric barrierMetric(
                  m_duStep, m_vStep, u, vK, p, barrierCoefficient);
                Detail::WNGIRPrimalBarrierForce barrierForce(
                  m_vStep, u, vK, p, barrierCoefficient);
                typename ProblemType::ProblemBodyType body(m_bulkForm);
                body = body + obsMetric + barrierMetric - surfaceForce - barrierForce;
                if (!p.dirichletAttributes.empty())
                  body = body +
                    Variational::DirichletBC(m_duStep, Variational::VectorFunction(zero))
                      .on(p.dirichletAttributes);
                m_stepProblem = body;
                m_stepProblem.assemble();
                rep.tAssembly += secondsSince(tic);

                tic = Clock::now();
                std::size_t barrierIterations = 0;
                Real barrierError = std::numeric_limits<Real>::infinity();
                solveOk = solveStep(uTrial, barrierIterations, barrierError);
                linearIterations += barrierIterations;
                linearError = barrierError;
                rep.tSolve += secondsSince(tic);
                if (!solveOk)
                  break;

                scratch = uTrial;
                scratch -= vK;
                const Real correctionNorm =
                  std::max(std::abs(scratch.max()), std::abs(scratch.min()));
                const Real relativeCorrection = predictorNorm > Real(0)
                  ? correctionNorm / predictorNorm
                  : (correctionNorm == Real(0) ? Real(0)
                                               : std::numeric_limits<Real>::infinity());
                rep.primalBarrierRelativeCorrection = relativeCorrection;
                ++rep.primalBarrierIterations;
                ++rep.lastPrimalBarrierIterations;
                const Real innerAlpha = getBarrierStepScale(
                  mesh, fes, validationCells, u, vK, scratch, meshDim);
                if (!(innerAlpha > Real(0)))
                {
                  solveOk = false;
                  break;
                }
                scratch *= innerAlpha;
                vK += scratch;
                innerConverged = p.primalBarrierRelativeTolerance > Real(0) &&
                  relativeCorrection <= p.primalBarrierRelativeTolerance &&
                  innerAlpha >= Real(0.999);
                tic = Clock::now();
                if (innerConverged)
                  break;
              }
              rep.primalBarrierConverged = innerConverged;
              if (solveOk && p.requirePrimalBarrierConvergence &&
                p.primalBarrierRelativeTolerance > Real(0) && !innerConverged)
              {
                rep.exitReason = "primal-barrier-inner-not-converged";
                solveOk = false;
              }
            }
          }
          if (!solveOk)
          {
            if (rep.exitReason == std::string_view("iter-budget"))
              rep.exitReason = "solve-linear-failed";
            break;
          }

          rep.linearIterations += linearIterations;
          rep.linearError = linearError;
          if (!(std::isfinite(vK.max()) && std::isfinite(vK.min())))
          {
            rep.exitReason = "solve-nonfinite";
            break;
          }

          Real directionAction = getSurfaceForceAction(mesh, fes, u, vK, phi, grad,
            interfaceFacets, sigma2, dataNormalization, meshDim, locator);
          const Real predictorNorm =
            std::max(std::abs(predictor.max()), std::abs(predictor.min()));
          Real directionNorm = std::max(std::abs(vK.max()), std::abs(vK.min()));
          const bool insufficientDescent = p.descentFraction > Real(0) &&
            predictorAction > Real(0) &&
            directionAction < p.descentFraction * predictorAction;
          const bool excessiveDirection = p.directionNormFactor > Real(0) &&
            predictorNorm > Real(0) &&
            directionNorm > p.directionNormFactor * predictorNorm;
          if ((insufficientDescent || excessiveDirection) && usePrimalBarrier)
          {
            vK = predictor;
            if (usePrimalBarrier)
            {
              uTrial *= Real(0);
              const Real predictorAlpha =
                getBarrierStepScale(mesh, fes, validationCells, u, uTrial, vK, meshDim);
              if (!(predictorAlpha > Real(0)))
              {
                rep.exitReason = "descent-fallback-not-feasible";
                break;
              }
              vK *= predictorAlpha;
            }
            directionAction = getSurfaceForceAction(mesh, fes, u, vK, phi, grad,
              interfaceFacets, sigma2, dataNormalization, meshDim, locator);
            directionNorm = std::max(std::abs(vK.max()), std::abs(vK.min()));
          }
          rep.directionAction = directionAction;
          rep.descentRatio =
            predictorAction > Real(0) ? directionAction / predictorAction : Real(0);
          rep.directionNormRatio =
            predictorNorm > Real(0) ? directionNorm / predictorNorm : Real(0);
          if (p.armijoCoefficient > Real(0) &&
            (!(directionAction > Real(0)) ||
              (p.descentFraction > Real(0) && predictorAction > Real(0) &&
                directionAction < p.descentFraction * predictorAction)))
          {
            rep.exitReason = "no-gradient-related-feasible-direction";
            break;
          }

          const Real maxStep = std::max(std::abs(vK.max()), std::abs(vK.min()));
          if (maxStep <= stepTol)
          {
            rep.exitReason = "best-effort-step-stagnation";
            break;
          }

          // ---- Nonlinear line search on TRUE geometry ----
          tic = Clock::now();
          const bool useLineSearch = checkGeometry || p.energyLineSearch;
          Real alpha = Real(1);
          bool accepted = false;
          std::size_t backtracks = 0;
          FastAdm adm{};
          Real eTrial = std::numeric_limits<Real>::infinity();
          SurfaceState trialSurface{};
          bool trialSurfaceEvaluated = false;
          if (!useLineSearch)
          {
            u += vK;
            accepted = true;
          }
          else
          {
            previousU = u;
            while (alpha >= p.alphaMin)
            {
              // uTrial = previousU + alpha * vK
              uTrial = vK;
              uTrial *= alpha;
              uTrial += previousU;
              bool jOK = true;
              bool qOK = true;
              if (checkGeometry)
              {
                adm = fastAdmissibility(uTrial);
                jOK = adm.inadmissibleCount == 0 && adm.minJ > acceptedJacobianFloor;
                qOK = adm.maxQ < p.qMax;
              }
              bool eOK = true;
              if (jOK && qOK && p.energyLineSearch)
              {
                trialSurface = surfaceState(uTrial);
                trialSurfaceEvaluated = true;
                eTrial = trialSurface.energy;
                const Real sufficientDecrease = p.armijoCoefficient > Real(0)
                  ? p.armijoCoefficient * alpha * directionAction
                  : Real(0);
                eOK = std::isfinite(eTrial) && eTrial <= ePrev - sufficientDecrease;
              }
              if (jOK && qOK && eOK)
              {
                u = uTrial;
                accepted = true;
                break;
              }
              if (!jOK)
                ++rep.jacobianRejections;
              if (!qOK)
                ++rep.distortionRejections;
              if (jOK && qOK && !eOK)
                ++rep.energyRejections;
              alpha *= Real(0.5);
              ++backtracks;
            }
          }
          rep.tLineSearch += secondsSince(tic);
          if (!accepted)
          {
            u = previousU;
            rep.exitReason = "line-search-failure";
            break;
          }

          FastAdm acceptedAdm = adm;
          SurfaceState acceptedSurf =
            trialSurfaceEvaluated ? trialSurface : surfaceState(u);
          Real acceptedEnergy =
            p.energyLineSearch && std::isfinite(eTrial) ? eTrial : acceptedSurf.energy;

          rep.lastAlpha = alpha;
          if (useLineSearch)
          {
            // acceptedStep = max_i |u_i - previousU_i|
            scratch = u;
            scratch -= previousU;
            rep.acceptedStep = std::max(std::abs(scratch.max()), std::abs(scratch.min()));
          }
          else
            rep.acceptedStep = maxStep;
          rep.minJ =
            checkGeometry ? acceptedAdm.minJ : std::numeric_limits<Real>::quiet_NaN();
          rep.maxJ =
            checkGeometry ? acceptedAdm.maxJ : std::numeric_limits<Real>::quiet_NaN();
          rep.maxQRel =
            checkGeometry ? acceptedAdm.maxQ : std::numeric_limits<Real>::quiet_NaN();

          const auto surf = acceptedSurf;
          currentSurface = surf;
          const Real eNow = acceptedEnergy;
          rep.backtracks += backtracks;
          rep.actualPredictedDecrease = alpha * directionAction > Real(0)
            ? (ePrev - eNow) / (alpha * directionAction)
            : Real(0);
          recordSurfaceState(surf);
          if (!(surf.activeLen > Real(0)))
          {
            rep.exitReason = "observation-degenerate-active-set";
            ++rep.iterations;
            break;
          }
          if (p.trace)
            std::cout << "      wngir it=" << std::setw(3) << rep.iterations
                      << "  E=" << std::scientific << std::setprecision(3) << eNow
                      << "  actRMS=" << surf.activeRMS << "  actRMS/(hG)="
                      << (levelSetMeshScale > Real(0) ? surf.activeRMS / levelSetMeshScale
                                                      : Real(0))
                      << "  actSup=" << surf.activeSup
                      << "  actFrac=" << rep.activeFraction
                      << "  step/h=" << (h > Real(0) ? rep.acceptedStep / h : Real(0))
                      << "  linIt=" << rep.linearIterations << "  alpha=" << alpha
                      << "  muEff=" << rep.primalBarrierCoefficient
                      << "  pbIt=" << rep.lastPrimalBarrierIterations
                      << "  pbRel=" << rep.primalBarrierRelativeCorrection
                      << "  stat=" << rep.stationarityNorm << "  desc=" << rep.descentRatio
                      << "  dirNorm=" << rep.directionNormRatio
                      << "  ared/pred=" << rep.actualPredictedDecrease
                      << "  bt=" << backtracks << "  rejJ=" << rep.jacobianRejections
                      << "  rejQ=" << rep.distortionRejections
                      << "  rejE=" << rep.energyRejections << "  min_j=" << rep.minJ
                      << "  max_j=" << rep.maxJ
                      << "  max_Q=" << rep.maxQRel << '\n';

          if (surf.activeRMS <= activeRMSTol)
          {
            rep.exitReason = "numerical-rms-converged";
            ++rep.iterations;
            break;
          }
          if (surf.activeSup <= activeSupTol)
          {
            rep.exitReason = "numerical-sup-converged";
            ++rep.iterations;
            break;
          }
          if (activeRMSOverHTol > Real(0) && levelSetMeshScale > Real(0) &&
            surf.activeRMS / levelSetMeshScale <= activeRMSOverHTol)
          {
            rep.exitReason = "geometric-rms-converged";
            ++rep.iterations;
            break;
          }
          if (activeSupOverHTol > Real(0) && levelSetMeshScale > Real(0) &&
            surf.activeSup / levelSetMeshScale <= activeSupOverHTol)
          {
            rep.exitReason = "geometric-sup-converged";
            ++rep.iterations;
            break;
          }
          if (p.acceptedStepOverHTol > Real(0) && h > Real(0) &&
            rep.acceptedStep / h <= p.acceptedStepOverHTol)
          {
            rep.exitReason = "best-effort-scale-step-stagnation";
            ++rep.iterations;
            break;
          }
          const Real eRel = std::abs(ePrev - eNow) / std::max(ePrev, Real(1e-30));
          if (eRel < p.energyStagTol)
          {
            rep.exitReason = "best-effort-energy-stagnation";
            ++rep.iterations;
            break;
          }
          ePrev = eNow;
        }

        const RigidModeState finalRigid = getRigidModeState(mesh, fes, u, phi, grad,
          interfaceFacets, sigma2, dataNormalization, meshDim, locator);
        rep.rigidModeCoercivity = finalRigid.minimum;
        rep.rigidModeCoercivityRatio = finalRigid.ratio;
        rep.rigidModeDimension = finalRigid.dimension;

        m_report = rep;
        return rep;
      }

    private:
      /**
       * @brief Quadrature formula WNGIR uses on a polytope.
       *
       * Exact for the products of shape functions the forms contain, unless an
       * order is pinned in the parameters. Serves cells and interface facets
       * alike: the polytope carries its own dimension, so the element degree
       * follows from it.
       */
      template <class FES>
      const QF::QuadratureFormulaBase& getQuadrature(
        const Geometry::Polytope& polytope, const FES& fes) const
      {
        const auto& fe =
          fes.getFiniteElement(polytope.getDimension(), polytope.getIndex());
        const std::size_t order = m_parameters.quadratureOrder > 0
          ? m_parameters.quadratureOrder
          : std::max<std::size_t>(2, 2 * fe.getOrder());
        return QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
      }

      /// @brief Measure of the fixed background domain used to normalize the QP.
      template <class Mesh, class FES>
      Real getDomainMeasure(
        const Mesh& mesh, const FES& fes, const std::vector<Index>& validationCells) const
      {
        if constexpr (requires { mesh.getShard(); })
        {
          return mesh.getMeasure(mesh.getDimension());
        }
        else
        {
          Real measure = 0;
          for (const Index cellIndex : validationCells)
          {
            const auto cell = mesh.getCell(cellIndex);
            const auto& qf = getQuadrature(*cell, fes);
            const auto& quadrature = cell->getQuadrature(qf);
            for (std::size_t q = 0; q < quadrature.getSize(); ++q)
              measure += qf.getWeight(q) * quadrature.getPoint(q).getDistortion();
          }
          return measure;
        }
      }

      /// @brief Action of the negative robust-energy first variation on a direction.
      template <class Mesh, class FES, class PhiType, class GradType, class LocatorType>
      Real getSurfaceForceAction(const Mesh& mesh, const FES& fes,
        const Displacement& current, const Displacement& direction, const PhiType& phi,
        const GradType& grad, const std::vector<Index>& interfaceFacets, Real sigma2,
        Real normalization, std::size_t dimension, const LocatorType& locator) const
      {
        if constexpr (requires { current.acquire(); })
          current.acquire();
        if constexpr (requires { direction.acquire(); })
          direction.acquire();
        if constexpr (requires { phi.acquire(); })
          phi.acquire();
        if constexpr (requires { grad.acquire(); })
          grad.acquire();
        std::vector<Real> facetActions(interfaceFacets.size(), Real(0));
#ifdef RODIN_USE_OPENMP
#pragma omp parallel
#endif
        {
          Detail::WNGIRSurfaceForceCoefficient force(
            phi, grad, current, locator, sigma2, normalization, dimension);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(interfaceFacets.size()); ++i)
          {
            Real& facetAction = facetActions[static_cast<std::size_t>(i)];
            const Index facetIndex = interfaceFacets[static_cast<std::size_t>(i)];
            const auto face = mesh.getFace(facetIndex);
            const auto& qf = getQuadrature(*face, fes);
            const auto& quadrature = face->getQuadrature(qf);
            for (std::size_t q = 0; q < quadrature.getSize(); ++q)
            {
              const auto& point = quadrature.getPoint(q);
              const Variational::IntegrationPoint ip(point, &qf, q);
              facetAction += qf.getWeight(q) * point.getDistortion() *
                force.getValue(ip).dot(direction.getValue(point));
            }
          }
        }
        Real action = 0;
        for (const Real facetAction : facetActions)
          action += facetAction;
        return action;
      }

      /**
       * @brief Observation coercivity on the rigid-motion kernel of the bulk form.
       *
       * The matrices are the sampled realizations of @f$G_R@f$ and @f$H_R@f$
       * from the analytical model. A full-vector Dirichlet condition on a
       * boundary attribute removes the rigid kernel, in which case the
       * diagnostic is vacuous.
       */
      template <class Mesh, class FES, class PhiType, class GradType, class LocatorType>
      RigidModeState getRigidModeState(const Mesh& mesh, const FES& fes,
        const Displacement& current, const PhiType& phi, const GradType& grad,
        const std::vector<Index>& interfaceFacets, Real sigma2, Real normalization,
        std::size_t dimension, const LocatorType& locator) const
      {
        if (!m_parameters.dirichletAttributes.empty())
          return {std::numeric_limits<Real>::infinity(), Real(1), 0};

        const std::size_t count = dimension * (dimension + 1) / 2;
        if (count == 0)
          return {std::numeric_limits<Real>::infinity(), Real(1), 0};

        std::vector<Math::Matrix<Real>> gradients(
          count, Math::Matrix<Real>::Zero(dimension, dimension));
        std::size_t mode = dimension;
        for (std::size_t a = 0; a < dimension; ++a)
        {
          for (std::size_t b = a + 1; b < dimension; ++b)
          {
            gradients[mode](static_cast<Eigen::Index>(a),
              static_cast<Eigen::Index>(b)) = Real(-1);
            gradients[mode](static_cast<Eigen::Index>(b),
              static_cast<Eigen::Index>(a)) = Real(1);
            ++mode;
          }
        }

        auto valuesAt = [&](const Math::SpatialPoint& point) {
          std::vector<SpatialVec> values(count, SpatialVec::Zero(dimension));
          for (std::size_t axis = 0; axis < dimension; ++axis)
            values[axis](static_cast<Eigen::Index>(axis)) = Real(1);
          std::size_t rotation = dimension;
          for (std::size_t a = 0; a < dimension; ++a)
          {
            for (std::size_t b = a + 1; b < dimension; ++b)
            {
              values[rotation](static_cast<Eigen::Index>(a)) =
                -point(static_cast<Eigen::Index>(b));
              values[rotation](static_cast<Eigen::Index>(b)) =
                point(static_cast<Eigen::Index>(a));
              ++rotation;
            }
          }
          return values;
        };

        Math::Matrix<Real> hGram = Math::Matrix<Real>::Zero(count, count);
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const auto& qf = getQuadrature(*cell, fes);
          const auto& quadrature = cell->getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& point = quadrature.getPoint(q);
            const Real weight = qf.getWeight(q) * point.getDistortion();
            const auto values = valuesAt(point.getCoordinates());
            for (std::size_t i = 0; i < count; ++i)
            {
              for (std::size_t j = 0; j < count; ++j)
              {
                hGram(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) +=
                  weight * (values[i].dot(values[j]) +
                    (gradients[i].array() * gradients[j].array()).sum());
              }
            }
          }
        }

        Math::Matrix<Real> observationGram = Math::Matrix<Real>::Zero(count, count);
        Detail::WNGIRObservationCoefficient observation(
          phi, grad, current, locator, m_parameters, sigma2, normalization, dimension);
        for (const Index facetIndex : interfaceFacets)
        {
          const auto face = mesh.getFace(facetIndex);
          const auto& qf = getQuadrature(*face, fes);
          const auto& quadrature = face->getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& point = quadrature.getPoint(q);
            const Variational::IntegrationPoint ip(point, &qf, q);
            const Real weight = qf.getWeight(q) * point.getDistortion();
            const auto values = valuesAt(point.getCoordinates());
            const SpatialMat coefficient = observation.getValue(ip);
            for (std::size_t i = 0; i < count; ++i)
            {
              for (std::size_t j = 0; j < count; ++j)
              {
                observationGram(
                  static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) +=
                  weight * values[i].dot(coefficient * values[j]);
              }
            }
          }
        }

        Eigen::GeneralizedSelfAdjointEigenSolver<Math::Matrix<Real>> eig(
          observationGram, hGram, Eigen::EigenvaluesOnly);
        if (eig.info() != Eigen::Success)
          return {std::numeric_limits<Real>::quiet_NaN(),
            std::numeric_limits<Real>::quiet_NaN(), count};
        const Real minimum = std::max(Real(0), eig.eigenvalues().minCoeff());
        const Real maximum = std::max(Real(0), eig.eigenvalues().maxCoeff());
        const Real ratio = maximum > Real(0) ? minimum / maximum : Real(0);
        return {minimum, ratio, count};
      }

      /// @brief Fraction-to-boundary factor for a primal barrier Newton update.
      template <class Mesh, class FES>
      Real getBarrierStepScale(const Mesh& mesh, const FES& fes,
        const std::vector<Index>& validationCells, const Displacement& current,
        const Displacement& inner, const Displacement& increment,
        std::size_t dimension) const
      {
        const Real tau =
          std::clamp(m_parameters.thetaBoundary, Real(0), Real(0.999));
        if constexpr (requires { current.acquire(); })
          current.acquire();
        if constexpr (requires { inner.acquire(); })
          inner.acquire();
        if constexpr (requires { increment.acquire(); })
          increment.acquire();
        Real alpha = Real(1);
        int feasible = 1;
#ifdef RODIN_USE_OPENMP
#pragma omp parallel reduction(min : alpha) reduction(& : feasible)
#endif
        {
          CellDeformation deformation(dimension);
          auto currentJacobian = Variational::Jacobian(current);
          auto innerJacobian = Variational::Jacobian(inner);
          auto incrementJacobian = Variational::Jacobian(increment);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(validationCells.size()); ++i)
          {
            const Index cellIndex = validationCells[static_cast<std::size_t>(i)];
            const auto cell = mesh.getCell(cellIndex);
            const auto& qf = getQuadrature(*cell, fes);
            const auto& quadrature = cell->getQuadrature(qf);
            for (std::size_t q = 0; q < quadrature.getSize(); ++q)
            {
              const Variational::IntegrationPoint ip(quadrature.getPoint(q), &qf, q);
              deformation.setDisplacementGradient(currentJacobian.getValue(ip));
              if (!deformation.isAdmissible())
              {
                feasible = 0;
                continue;
              }
              const auto innerGradient = innerJacobian.getValue(ip);
              const auto incrementGradient = incrementJacobian.getValue(ip);
              const Real innerJ = -deformation.getJacobianAction(innerGradient);
              const Real incrementJ = -deformation.getJacobianAction(incrementGradient);
              const Real innerQ = deformation.getRelativeDistortionAction(innerGradient);
              const Real incrementQ =
                deformation.getRelativeDistortionAction(incrementGradient);
              const Real slackJ = deformation.getJacobian() - m_parameters.jSafe - innerJ;
              const Real slackQ =
                m_parameters.qMax - deformation.getRelativeDistortion() - innerQ;
              if (slackJ <= Real(0) || slackQ <= Real(0))
              {
                feasible = 0;
                continue;
              }
              if (incrementJ > Real(0))
                alpha = std::min(alpha, tau * slackJ / incrementJ);
              if (incrementQ > Real(0))
                alpha = std::min(alpha, tau * slackQ / incrementQ);
            }
          }
        }
        if (!feasible)
          return Real(0);
        return std::clamp(alpha, Real(0), Real(1));
      }

      /**
       * @brief Evaluates the mesh validity of a displacement.
       *
       * Sweeps every cell quadrature point for the extremes that decide whether
       * a trial step is acceptable. Only @f$j@f$ is needed at most points, so
       * the deformation is evaluated lazily, and the distortion is read only
       * where the cell is not inverted, since it is undefined otherwise.
       */
      template <class Mesh, class FES>
      AdmissibilityState getAdmissibilityState(const Mesh& mesh, const FES& fes,
        const std::vector<Index>& validationCells, const Displacement& u,
        std::size_t dimension) const
      {
        if constexpr (requires { u.acquire(); })
          u.acquire();
        Real minJ = std::numeric_limits<Real>::infinity();
        Real maxJ = -std::numeric_limits<Real>::infinity();
        Real maxQ = Real(0);
        std::size_t inadmissibleCount = 0;
#ifdef RODIN_USE_OPENMP
#pragma omp parallel reduction(min : minJ) reduction(max : maxJ, maxQ)                   \
  reduction(+ : inadmissibleCount)
#endif
        {
          CellDeformation deformation(dimension);
          auto displacementJacobian = Variational::Jacobian(u);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(validationCells.size()); ++i)
          {
            const Index cellIndex = validationCells[static_cast<std::size_t>(i)];
            const auto cell = mesh.getCell(cellIndex);
            const auto& qf = getQuadrature(*cell, fes);
            const auto& quad = cell->getQuadrature(qf);
            for (std::size_t q = 0; q < quad.getSize(); ++q)
            {
              const Variational::IntegrationPoint ip(quad.getPoint(q), &qf, q);
              deformation.setDisplacementGradient(displacementJacobian.getValue(ip));
              const Real j = deformation.getJacobian();
              minJ = std::min(minJ, j);
              maxJ = std::max(maxJ, j);
              if (j <= m_parameters.jMinRatio)
                ++inadmissibleCount;
              if (deformation.isAdmissible())
                maxQ = std::max(maxQ, deformation.getRelativeDistortion());
            }
          }
        }
        return {minJ, maxJ, maxQ, inadmissibleCount};
      }

      /**
       * @brief Evaluates the interface fit of a displacement.
       *
       * The residual is the level set read at the deformed image of each
       * interface quadrature point. Points whose robust weight has fallen below
       * @ref WNGIRParameters::omegaMin are excluded from the residual norms:
       * they are the ones the robust weighting has already rejected as
       * outliers, and including them would let a distant feature dominate an
       * otherwise converged fit.
       */
      template <class Mesh, class FES, class PhiType, class LocatorType>
      SurfaceState getSurfaceState(const Mesh& mesh, const FES& fes,
        const Displacement& u, const PhiType& phi,
        const std::vector<Index>& interfaceFacets, const WNGIRLoss& loss,
        Real normalization, std::size_t dimension, const LocatorType& locator) const
      {
        if constexpr (requires { u.acquire(); })
          u.acquire();
        if constexpr (requires { phi.acquire(); })
          phi.acquire();
        struct SurfaceAccumulation
        {
            SurfaceState state;
            Real squaredResidual = 0;
        };
        std::vector<SurfaceAccumulation> facetStates(interfaceFacets.size());
#ifdef RODIN_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (Index i = 0; i < static_cast<Index>(interfaceFacets.size()); ++i)
        {
          auto& accumulation = facetStates[static_cast<std::size_t>(i)];
          SurfaceState& facetState = accumulation.state;
          DeformationMap deformation(u, locator);
          const Index facetIdx = interfaceFacets[static_cast<std::size_t>(i)];
          const auto face = mesh.getFace(facetIdx);
          const auto& qf = getQuadrature(*face, fes);
          const auto& quad = face->getQuadrature(qf);
          for (std::size_t q = 0; q < quad.getSize(); ++q)
          {
            const auto& src = quad.getPoint(q);
            const Variational::IntegrationPoint ip(src, &qf, q);
            const Real w = qf.getWeight(q) * src.getDistortion();
            const Real r = phi.getValue(deformation.getMovedPoint(ip));
            const Real omega = loss.getWeight(r);
            facetState.totalLen += w;
            facetState.energy += w * normalization * loss.getValue(r);
            if (omega >= m_parameters.omegaMin)
            {
              facetState.activeLen += w;
              accumulation.squaredResidual += w * r * r;
              facetState.activeSup = std::max(facetState.activeSup, std::abs(r));
            }
          }
        }

        SurfaceState state;
        Real squared = 0;
        // Combine the per-facet positive quadrature sums before normalizing the
        // active residual over the complete interface.
        for (const auto& accumulation : facetStates)
        {
          const SurfaceState& facetState = accumulation.state;
          state.energy += facetState.energy;
          state.activeLen += facetState.activeLen;
          state.totalLen += facetState.totalLen;
          state.activeSup = std::max(state.activeSup, facetState.activeSup);
          squared += accumulation.squaredResidual;
        }
        state.activeRMS = state.activeLen > Real(0)
          ? std::sqrt(std::max(Real(0), squared) / state.activeLen)
          : std::numeric_limits<Real>::infinity();
        if (!(state.activeLen > Real(0)))
          state.activeSup = std::numeric_limits<Real>::infinity();
        return state;
      }

      /**
       * @brief Normal-jump statistics of the discrete interface.
       *
       * Facets are adjacent when they share a vertex in two dimensions or an
       * edge in three, and each adjacent pair contributes the angle between its
       * normals weighted by the mean of the two facet measures. Normals are
       * compared up to sign, since facet orientation is not consistent across
       * the interface.
       */
      template <class Mesh, class FES>
      NormalJump getNormalJump(const Mesh& mesh, const FES& fes,
        const std::vector<Index>& interfaceFacets, std::size_t dimension) const
      {
        // Averaged unit normal and measure of each facet. The normal is the
        // quadrature average of the pointwise unit normals, so a curved facet
        // contributes its mean orientation; a facet degenerate throughout
        // falls back to a fixed direction rather than a null vector.
        constexpr Real degenerate = Real(1e-30);
        std::vector<Math::SpatialVector<Real>> normals(interfaceFacets.size());
        std::vector<Real> measures(interfaceFacets.size(), Real(0));

        for (std::size_t i = 0; i < interfaceFacets.size(); ++i)
        {
          const auto face = mesh.getFace(interfaceFacets[i]);
          const auto& qf = getQuadrature(*face, fes);
          const auto& quad = face->getQuadrature(qf);

          Math::SpatialVector<Real> n = Math::SpatialVector<Real>::Zero(dimension);
          Real measure = 0;
          for (std::size_t q = 0; q < quad.getSize(); ++q)
          {
            const auto& pt = quad.getPoint(q);
            const Real w = qf.getWeight(q) * pt.getDistortion();
            const auto& J = pt.getJacobian();
            // The facet normal: in 2D the tangent rotated a quarter turn, in
            // 3D the cross product of the two tangents.
            Math::SpatialVector<Real> nq;
            if (dimension == 2)
            {
              nq.resize(2);
              nq(0) = J(1, 0);
              nq(1) = -J(0, 0);
            }
            else
            {
              Math::SpatialVector<Real> a = Math::SpatialVector<Real>::Zero(dimension);
              Math::SpatialVector<Real> b = Math::SpatialVector<Real>::Zero(dimension);
              for (std::size_t r = 0; r < dimension; ++r)
              {
                a(static_cast<Eigen::Index>(r)) = J(r, 0);
                b(static_cast<Eigen::Index>(r)) = J(r, 1);
              }
              nq = a.cross(b);
            }
            const Real norm = nq.norm();
            if (norm > degenerate)
            {
              n += w * (nq / norm);
              measure += w;
            }
          }

          if (n.norm() <= degenerate)
          {
            n.setZero();
            n(0) = Real(1);
          }
          else
          {
            n /= n.norm();
          }
          measures[i] = measure;
          normals[i] = n;
        }

        const auto edgeKey = [](Index a, Index b) -> std::uint64_t {
          const auto lo = static_cast<std::uint64_t>(std::min(a, b));
          const auto hi = static_cast<std::uint64_t>(std::max(a, b));
          return (lo << 32) ^ hi;
        };

        // Facets sharing an incidence, keyed by the shared entity.
        std::unordered_map<std::uint64_t, std::vector<std::size_t>> incident;
        if (dimension == 2)
        {
          incident.reserve(2 * interfaceFacets.size());
          for (std::size_t i = 0; i < interfaceFacets.size(); ++i)
            for (const Index v : mesh.getFace(interfaceFacets[i])->getVertices())
              incident[static_cast<std::uint64_t>(v)].push_back(i);
        }
        else if (dimension == 3)
        {
          incident.reserve(3 * interfaceFacets.size());
          for (std::size_t i = 0; i < interfaceFacets.size(); ++i)
          {
            const auto& vv = mesh.getFace(interfaceFacets[i])->getVertices();
            incident[edgeKey(vv[0], vv[1])].push_back(i);
            incident[edgeKey(vv[1], vv[2])].push_back(i);
            incident[edgeKey(vv[2], vv[0])].push_back(i);
          }
        }

        NormalJump jump;
        Real squared = 0;
        Real weight = 0;
        for (const auto& [key, faces] : incident)
        {
          (void)key;
          for (std::size_t a = 0; a < faces.size(); ++a)
          {
            for (std::size_t b = a + 1; b < faces.size(); ++b)
            {
              const auto ia = faces[a], ib = faces[b];
              const Real dot = std::abs(normals[ia].dot(normals[ib]));
              const Real theta = std::acos(std::max(Real(-1), std::min(Real(1), dot)));
              const Real w = Real(0.5) * (measures[ia] + measures[ib]);
              squared += w * theta * theta;
              weight += w;
              jump.max = std::max(jump.max, theta);
            }
          }
        }
        jump.rms = weight > Real(0) ? std::sqrt(squared / weight) : Real(0);
        return jump;
      }

      /**
       * @brief Robust scale @f$\sigma@f$ of the interface residual.
       *
       * The loss weight requires a scale separating the residuals to be fitted
       * from those treated as outliers. It is the 90th percentile of the
       * undeformed residual, floored at @f$3hG_\phi@f$, where @f$G_\phi@f$ is
       * the maximum sampled target-gradient norm. The percentile adapts to the
       * initial mismatch, while the floor converts the geometric resolution to
       * level-set units.
       *
       * Fixed once per frame rather than re-estimated per iteration, which
       * would let the scale chase its own progress and never reject anything.
       * A positive @ref WNGIRParameters::robustScale overrides the automatic
       * mesh-dependent selection.
       */
      struct RobustScale
      {
          /// @brief Fixed robust-loss scale in level-set units.
          Real sigma = 0;
          /// @brief Maximum sampled target-gradient norm.
          Real gradientScale = 0;
      };

      template <class Mesh, class FES, class PhiType, class GradType>
      RobustScale getRobustScale(const Mesh& mesh, const FES& fes, const PhiType& phi,
        const GradType& grad, const std::vector<Index>& interfaceFacets, Real h) const
      {
        std::vector<Real> residuals;
        Real gradientScale = 0;
        for (const Index facetIdx : interfaceFacets)
        {
          const auto face = mesh.getFace(facetIdx);
          const auto& qf = getQuadrature(*face, fes);
          const auto& quad = face->getQuadrature(qf);
          for (std::size_t q = 0; q < quad.getSize(); ++q)
          {
            const auto& point = quad.getPoint(q);
            residuals.push_back(std::abs(phi.getValue(point)));
            gradientScale = std::max(gradientScale, grad.getValue(point).norm());
          }
        }

        Real sigma = m_parameters.robustScale;
        if (!(sigma > Real(0)))
        {
          sigma = Real(3) * h * gradientScale;
          const std::size_t k90 =
            static_cast<std::size_t>(Real(0.9) * static_cast<Real>(residuals.size() - 1));
          std::nth_element(residuals.begin(), residuals.begin() + k90, residuals.end());
          sigma = std::max(sigma, residuals[k90]);
        }
        sigma = std::max(sigma, std::sqrt(std::numeric_limits<Real>::min()));
        return {sigma, gradientScale};
      }

      // Solves the currently-assembled step problem with CG and copies the
      // backend-matched solution GridFunction into @p out.
      bool solveStep(Displacement& out, std::size_t& iterations, Real& error)
      {
        m_stepProblem.solve(m_stepSolver);
        out = m_duStep.getSolution();
        if constexpr (requires(
                        const StepSolverType& solver) { solver.getIterationNumber(); })
          iterations = m_stepSolver.getIterationNumber();
        else
          iterations = 0;
        if constexpr (requires(const StepSolverType& solver) { solver.getError(); })
          error = m_stepSolver.getError();
        else
          error = Real(0);
        const bool success = [&]() {
          if constexpr (requires(const StepSolverType& solver) { solver.success(); })
            return static_cast<bool>(m_stepSolver.success());
          else
            return true;
        }();
        return success && std::isfinite(out.max()) && std::isfinite(out.min());
      }

      Displacement* m_u;
      TrialFunctionType m_duStep;
      TestFunctionType m_vStep;
      ProblemType m_stepProblem;
      StepSolverType m_stepSolver;
      BilinearFormType m_bulkForm;
      bool m_bulkFormAssembled = false;
      WNGIRParameters m_parameters;
      WNGIRReport m_report;
  };
}

#endif
