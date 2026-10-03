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
#include <utility>
#include <vector>

#include <Eigen/Eigenvalues>
#include <Eigen/IterativeLinearSolvers>
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Solver/MUMPS.h"
#include "Rodin/Variational/Integrator.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/Trace.h"
#include "Rodin/Variational/IdentityMatrix.h"
#include "Rodin/Variational/VectorFunction.h"

#include "Rodin/Location.h"

#include "CellDeformation.h"
#include "WNGIRLoss.h"
#include "WNGIRParameters.h"
#include "WNGIRDirectionalNewton.h"
#include "WNGIRRegularityMetric.h"
#include "WNGIRPrimalBarrierForce.h"
#include "WNGIRPrimalBarrierMetric.h"
#include "WNGIRReport.h"
#include "WNGIRObservationCoefficient.h"
#include "WNGIRSurfaceForceCoefficient.h"

namespace Rodin::Adaptation
{
  namespace Detail
  {
    /// @brief Sampled componentwise physical displacement norm, independent of FE coefficients.
    template <class Mesh, class FES, class Displacement>
    Real wngirPhysicalDisplacementNorm(const Mesh& mesh, const FES& fes,
      const std::vector<Index>& cells, const Displacement& step,
      std::size_t validationOrder = 0)
    {
      if constexpr (requires { step.acquire(); })
        step.acquire();
      Real maximum = 0;
      for (const Index index : cells)
      {
        const auto cell = mesh.getCell(index);
        const Geometry::Polytope::Traits traits(cell->getGeometry());
        const auto evaluate = [&](const Geometry::Point& point) {
          const auto value = step.getValue(point);
          Real norm = 0;
          for (std::size_t component = 0; component < mesh.getDimension(); ++component)
          {
            if (!std::isfinite(value(component)))
              return std::numeric_limits<Real>::infinity();
            norm = std::max(norm, std::abs(value(component)));
          }
          return norm;
        };
        for (std::size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
          maximum = std::max(maximum, evaluate(Geometry::Point(*cell, traits.getVertex(vertex))));
        const auto& fe = fes.getFiniteElement(cell->getDimension(), index);
        // This is a polynomial FE field, not the composed level-set residual.
        const auto& qf = QF::PolytopeQuadratureFormula::get(validationOrder > 0
            ? validationOrder
            : std::max<size_t>(6, 2 * fe.getOrder() + 4),
          cell->getGeometry());
        const auto& quadrature = cell->getQuadrature(qf);
        for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          maximum = std::max(maximum, evaluate(quadrature.getPoint(q)));
      }
      return maximum;
    }

  }

  /**
   * @brief Fitting--distribution WNGIR mesh-fitting solver for the local backend.
   *
   * Fits the interface skeleton of a mesh to the zero level set of @f$\phi@f$
   * with M = F + D (Fitting, Distribution) and affine quadratic hinges.
   * The independent weights are kappaF and kappaD. Each outer iteration
   * freezes the metric, force and constraint rows; inner Newton uses merit
   * backtracking. Directional Newton scales the physical predictor and metric
   * before constructing the hinges. Armijo then backtracks from the full increment.
   *
   * The form-language assembly uses the local Eigen backend, with CG, SparseLU
   * or MUMPS solving the same pointwise deviatoric current-strain operator. Linear systems select zero
   * along unresolved similarity modes; this does not modify the inner objective.
   * There is no full-space completion or inertia gate. Linear
   * residual, direction and actual-geometry line-search checks remain active.
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
      using BilinearFormType = std::decay_t<decltype(Variational::BilinearForm(
        std::declval<TrialFunctionType&>(), std::declval<TestFunctionType&>()))>;
      using LinearFormType = std::decay_t<decltype(Variational::LinearForm(
        std::declval<TestFunctionType&>()))>;

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
          Real affineMinJ = std::numeric_limits<Real>::infinity();
          Real affineMaxQ = 0;
          Index qualityCell = 0;
          Real qualityCurrent = 0;
          Real qualityLinearChange = 0;
          Real qualityQuadraticChange = 0;
      };

      /// @brief Interface-fit summary of a displacement.
      struct SurfaceState
      {
          /// @brief Fixed-scale energy used by the current Armijo search.
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
       * These statistics describe the classified interface only; they do not
       * alter the geometric target or the solver tolerances.
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
          m_regularityForm(m_duStep, m_vStep),
          m_localMetricForm(m_duStep, m_vStep),
          m_surfaceForm(m_vStep)
      {}

      /// @brief Sets WNGIR runtime parameters.
      WNGIR& setParameters(const WNGIRParameters& parameters)
      {
        using OperatorType =
          typename FormLanguage::Traits<LinearSystemType>::OperatorType;
        using VectorType = typename FormLanguage::Traits<LinearSystemType>::VectorType;
        if constexpr (!std::is_same_v<OperatorType, Math::SparseMatrix<Real>> ||
          !std::is_same_v<VectorType, Math::Vector<Real>>)
          Alert::Exception() << "WNGIR requires the local Eigen backend." << Alert::Raise;
#ifndef RODIN_USE_MUMPS
        if (parameters.directSolver == WNGIRParameters::DirectSolver::MUMPS)
          Alert::Exception() << "WNGIR MUMPS solves require RODIN_USE_MUMPS."
                             << Alert::Raise;
#endif
        if (!std::isfinite(parameters.kappaD) || !(parameters.kappaD > Real(0)) ||
          !std::isfinite(parameters.kappaF) || !(parameters.kappaF > Real(0)))
          Alert::Exception() << "WNGIR requires positive "
                                "fitting/distribution weights."
                             << Alert::Raise;
        for (const Real value : {parameters.kappaJ, parameters.kappaQ,
               parameters.primalBarrierAbsoluteTolerance, parameters.energyStagTol,
               parameters.stepTol, parameters.acceptedStepOverHTol, parameters.robustScale})
          if (!std::isfinite(value) || value < Real(0))
            Alert::Exception() << "WNGIR weights and tolerances must be finite and nonnegative."
                               << Alert::Raise;
        if (!std::isfinite(parameters.primalBarrierRelativeTolerance) ||
          !(parameters.primalBarrierRelativeTolerance > Real(0)) ||
          !std::isfinite(parameters.cgRelativeTolerance) ||
          !(parameters.cgRelativeTolerance > Real(0)) ||
          !(parameters.armijoCoefficient > Real(0) &&
            parameters.armijoCoefficient < Real(1)) ||
          !(parameters.omegaMin > Real(0) && parameters.omegaMin <= Real(1)) ||
          !(parameters.jMinRatio > Real(0) && parameters.jMinRatio < Real(1)) ||
          !(parameters.jLineSearchRatio > Real(0) &&
            parameters.jLineSearchRatio < Real(1)) ||
          parameters.primalBarrierIterations == 0 || parameters.stagnationIterations == 0)
          Alert::Exception() << "WNGIR requires valid residual tolerances, budgets and line search."
                             << Alert::Raise;
        if (!std::isfinite(parameters.directionalNewtonMaxStepOverH) ||
          !(parameters.directionalNewtonMaxStepOverH > Real(0)))
          Alert::Exception()
            << "Directional Newton requires a finite positive step/h bound."
            << Alert::Raise;
        if (!(parameters.qualityGuard > Real(0) && parameters.qualityGuard < Real(1)) ||
          !(parameters.jSafe > Real(0) && parameters.jSafe < Real(1)) ||
          !(parameters.qMax > Real(1)) || !std::isfinite(parameters.qMax) ||
          !std::isfinite(parameters.muHat) || parameters.muHat < Real(0) ||
          !std::isfinite(parameters.geometricSupTolerance) ||
          parameters.geometricSupTolerance < Real(0))
          Alert::Exception() << "WNGIR requires 0 < guard,jSafe < 1, finite qMax > 1 and "
                                "nonnegative penalty/fit tolerance."
                             << Alert::Raise;
        m_parameters = parameters;
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

        WNGIRReport rep;
        const Real h = p.h;
        if (!std::isfinite(h) || !(h > Real(0)))
          Alert::Exception() << "WNGIR requires a finite positive reference mesh size."
                             << Alert::Raise;
        if (&mesh != &fes.getMesh())
          Alert::Exception() << "WNGIR mesh must match the displacement space."
                             << Alert::Raise;
        const Real acceptedJacobianFloor = std::max(p.jLineSearchRatio, p.jSafe);
        const Real stepTol = p.stepTol;
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
        std::vector<Index> markedFacets;
        for (auto face = mesh.getFace(); face; ++face)
          if (face->getAttribute() == p.interfaceAttribute)
            markedFacets.push_back(face->getIndex());
        auto suppliedFacets = interfaceFacets;
        std::sort(markedFacets.begin(), markedFacets.end());
        std::sort(suppliedFacets.begin(), suppliedFacets.end());
        if (markedFacets != suppliedFacets)
          Alert::Exception() << "WNGIR facet list must exactly match the interface attribute."
                             << Alert::Raise;
        const auto normalJump = getNormalJump(mesh, fes, interfaceFacets, meshDim);
        rep.normalJumpRMS = normalJump.rms;
        rep.normalJumpMax = normalJump.max;
        const Real geometricTarget = p.geometricSupTolerance > Real(0)
          ? p.geometricSupTolerance
          : std::pow(h, fes.getFiniteElement(meshDim - 1, interfaceFacets.front()).getOrder() + 1);
        if (!std::isfinite(geometricTarget) || !(geometricTarget > Real(0)))
          Alert::Exception() << "WNGIR requires a finite positive geometric target."
                             << Alert::Raise;
        rep.geometricSupTarget = geometricTarget;
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
        const WNGIRLoss loss(sigma);
        rep.sigma = sigma;
        rep.levelSetGradientScale = robustScale.gradientScale;
        const Real domainMeasure = getDomainMeasure(mesh, fes, validationCells);

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
          return wngirInterfaceQuadratureOrder(fe.getOrder());
        };

        FastAdm initialAdm{};
        {
          initialAdm = fastAdmissibility(u);
          rep.minJ = initialAdm.minJ;
          rep.maxJ = initialAdm.maxJ;
          rep.maxQRel = initialAdm.maxQ;
        }
        if (initialAdm.inadmissibleCount > 0 || initialAdm.minJ <= acceptedJacobianFloor ||
          initialAdm.maxQ >= p.qMax)
        {
          rep.exitReason = "initial-state-not-strictly-feasible";
          m_report = rep;
          return rep;
        }

        SurfaceState currentSurface = surfaceState(u);
        recordSurfaceState(currentSurface);
        if (p.rigidDiagnostics)
        {
          const RigidModeState initialRigid = getRigidModeState(mesh, fes, u, phi, grad,
            interfaceFacets, sigma2, dataNormalization, meshDim, locator);
          rep.rigidModeCoercivity = initialRigid.minimum;
          rep.rigidModeCoercivityRatio = initialRigid.ratio;
          rep.rigidModeDimension = initialRigid.dimension;
        }
        if (!(currentSurface.activeLen > Real(0)))
        {
          rep.exitReason = "observation-degenerate-active-set";
          m_report = rep;
          return rep;
        }

        Real ePrev = currentSurface.energy;
        const auto traceStart = setupTic;
        // Validation concerns the accepted outer geometry, not the inner QP iterate.
        const auto traceGeometry = [&](const auto& geometry, std::size_t accepted,
                                     const char* phase) {
          const auto flags = std::cout.flags();
          const auto precision = std::cout.precision();
          std::cout << "      wngir geometry: outer=" << accepted << "  phase=" << phase
                    << std::scientific
                    << std::setprecision(std::numeric_limits<Real>::max_digits10)
                    << "  geom_rms=" << geometry.rms << "  geom_sup=" << geometry.sup
                    << "  normal_rms=" << geometry.normalRMS << "  h=" << h
                    << "  min_j=" << rep.minJ << "  max_qrel=" << rep.maxQRel
                    << "  inner_total=" << rep.primalBarrierIterations
                    << "  inner_last=" << rep.lastPrimalBarrierIterations
                    << "  inner_converged=" << rep.primalBarrierConverged
                    << "  inner_residual=" << rep.primalBarrierResidual
                    << "  inner_relative_residual=" << rep.primalBarrierRelativeResidual
                    << "  inner_residual_tolerance=" << rep.primalBarrierResidualTolerance
                    << "  geom_sup_target=" << geometricTarget
                    << "  inactive_hinge_skips=" << rep.inactiveHingeSkips
                    << "  direct_analyses=" << rep.directAnalyses
                    << "  direct_factorizations=" << rep.directFactorizations
                    << "  fit_energy=" << rep.energy
                    << "  inner_assembly_seconds=" << rep.tPrimalBarrierAssembly
                    << "  inner_solve_seconds=" << rep.tPrimalBarrierSolve
                    << "  inner_ls_seconds=" << rep.tPrimalBarrierLineSearch
                    << "  seconds=" << secondsSince(traceStart) << std::endl;
          std::cout.flags(flags);
          std::cout.precision(precision);
        };
        {
          const auto geometry = getInterfaceGeometryState(
            mesh, fes, u, phi, grad, interfaceFacets, meshDim, locator);
          if (p.trace)
            traceGeometry(geometry, 0, "initial");
          rep.geometricRMS = geometry.rms;
          rep.geometricSup = geometry.sup;
          rep.normalRMS = geometry.normalRMS;
          rep.qualityBudgetSatisfied = true;
          if (!std::isfinite(geometry.sup))
          {
            rep.exitReason = "geometric-validation-failed";
            m_report = rep;
            return rep;
          }
          if (geometry.sup <= geometricTarget)
          {
            rep.geometricRMS = geometry.rms;
            rep.geometricSup = geometry.sup;
            rep.normalRMS = geometry.normalRMS;
            rep.geometricTargetReached = true;
            rep.qualityBudgetSatisfied = true;
            rep.exitReason = "full-interface-geometric-sup-converged";
            m_report = rep;
            return rep;
          }
        }

        // ============================================================
        // Nonlinear iteration.
        // ============================================================
        std::size_t consecutiveSmallAcceptedSteps = 0;
        std::size_t consecutiveSmallEnergyChanges = 0;
        const auto recordLinearSolve = [&rep](std::size_t iterations, Real error) {
          rep.linearIterations += iterations;
          ++rep.linearSolveCount;
          rep.maxLinearIterations = std::max(rep.maxLinearIterations, iterations);
          rep.linearError = error;
        };
        for (; rep.iterations < p.maxIterations; ++rep.iterations)
        {
          auto tic = Clock::now();
          Detail::WNGIRObservationCoefficient obsCoeff(
            grad, u, locator, p, dataNormalization, meshDim);
          auto obsMetric =
            Variational::FaceIntegral(Variational::Dot(obsCoeff * m_duStep, m_vStep));
          obsMetric.setOrder(surfaceOrder);
          obsMetric.over(p.interfaceAttribute);
          Detail::WNGIRSurfaceForceCoefficient forceCoeff(
            phi, grad, u, locator, loss, dataNormalization, meshDim);
          auto surfaceForce = Variational::FaceIntegral(forceCoeff, m_vStep);
          surfaceForce.setOrder(surfaceOrder);
          surfaceForce.over(p.interfaceAttribute);
          // The observation metric and the fitting force depend on the outer
          // displacement, not on the barrier increment, so they are assembled
          // here and reused by every correction below.
          m_localMetricForm = obsMetric;
          m_localMetricForm.assemble();
          const Real coefficient = p.h * p.kappaD;
          size_t order = 2;
          for (auto cell = mesh.getCell(); cell; ++cell)
            order = std::max(
              order, 2 * fes.getFiniteElement(meshDim, cell->getIndex()).getOrder());
          if (p.quadratureOrder > 0)
            order = p.quadratureOrder;
          m_regularityForm = Detail::WNGIRCurrentStrainMetric(
            m_duStep, m_vStep, u, coefficient, order);
          m_regularityForm.assemble();
          m_surfaceForm = surfaceForce;
          m_surfaceForm.assemble();
          bool solveOk = true;
          Real predictorAction = Real(0);
          typename ProblemType::ProblemBodyType predictorBody(m_regularityForm);
          predictorBody = predictorBody + m_localMetricForm - m_surfaceForm;
          m_stepProblem = predictorBody;
          assembleStep();
          rep.tAssembly += secondsSince(tic);

          Math::SparseMatrix<Real> fixedMetric;
          Math::Vector<Real> fixedForce;
          RegularityProjection fixedProjection;
          if constexpr (std::is_same_v<
                          typename FormLanguage::Traits<LinearSystemType>::OperatorType,
                          Math::SparseMatrix<Real>>)
          {
            {
              const auto& system = m_stepProblem.getLinearSystem();
              fixedMetric = system.getOperator();
              fixedForce = system.getVector();
              m_similarityModes = getSimilarityModes(u);
            }
          }
          tic = Clock::now();
          std::size_t predictorIterations = 0;
          Real predictorError = std::numeric_limits<Real>::infinity();
          solveOk = solveStep(vK, predictorIterations, predictorError, &fixedProjection, &rep);
          recordLinearSolve(predictorIterations, predictorError);
          rep.tSolve += secondsSince(tic);
          if (!solveOk)
          {
            rep.exitReason = "solve-predictor-failed";
            break;
          }

          predictor = vK;
          predictorAction = getSurfaceForceAction(mesh, fes, u, predictor, phi, grad,
            interfaceFacets, loss, dataNormalization, meshDim, locator);
          if (!(predictorAction > Real(0)) || !std::isfinite(predictorAction))
          {
            rep.exitReason = "no-descent-predictor";
            break;
          }
          rep.predictorScale = Real(1);
          if (p.directionalNewton)
          {
            const auto curvatures =
              getSurfaceDirectionalCurvature(mesh, fes, u, predictor, phi, grad,
                interfaceFacets, loss, dataNormalization, meshDim, locator);
            const Real norm = Detail::wngirPhysicalDisplacementNorm(
              mesh, fes, validationCells, predictor, p.geometricValidationOrder);
            rep.predictorScale =
              Detail::wngirDirectionalNewtonStep(predictorAction, curvatures.first,
                curvatures.second, norm, h * p.directionalNewtonMaxStepOverH);
            if (!(rep.predictorScale > Real(0)) || !std::isfinite(rep.predictorScale))
            {
              rep.exitReason = "invalid-directional-scaling";
              break;
            }
            // The predictor, frozen quadratic form and hinge scale describe the same motion.
            predictor *= rep.predictorScale;
            vK = predictor;
            predictorAction *= rep.predictorScale;
            fixedMetric *= Real(1) / rep.predictorScale;
            m_regularityForm.getOperator() *= Real(1) / rep.predictorScale;
            m_localMetricForm.getOperator() *= Real(1) / rep.predictorScale;
            for (Real& weight : fixedProjection.weights)
              weight /= rep.predictorScale;
            if (p.trace)
            {
              const auto precision = std::cout.precision();
              std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
                        << "      wngir directional: outer=" << rep.iterations
                        << "  curvature=" << curvatures.first
                        << "  fitting_curvature=" << curvatures.second
                        << "  scale=" << rep.predictorScale
                        << "  step_over_h=" << norm * rep.predictorScale / h
                        << "  unresolved_similarity_modes="
                        << rep.unresolvedSimilarityModes << '\n';
              std::cout.precision(precision);
            }
          }
          rep.predictorAction = predictorAction;

          {
            const Real modelDecrease = Real(0.5) * std::max(Real(0), predictorAction);
            const Real barrierCoefficient =
              domainMeasure > Real(0) ? p.muHat * modelDecrease / domainMeasure : Real(0);
            rep.primalBarrierCoefficient = barrierCoefficient;

            {
              // The affine hinges permit any finite predictor increment.
              const std::size_t innerIterations =
                std::max<std::size_t>(1, p.primalBarrierIterations);
              const Real predictorNorm =
                std::max(std::abs(predictor.max()), std::abs(predictor.min()));
              bool innerConverged = false;
              rep.lastPrimalBarrierIterations = 0;
              rep.lastPrimalBarrierAlpha = 0;
              rep.minPrimalBarrierAlpha = 1;
              rep.fullPrimalBarrierSteps = 0;
              tic = Clock::now();
              const Real residualScale = fixedForce.norm();
              for (std::size_t inner = 0; inner <= innerIterations; ++inner)
              {
                if (inner == 0 &&
                  !hasActivePrimalHinges(mesh, fes, validationCells, u, vK,
                    meshDim, barrierCoefficient))
                {
                  Math::Vector<Real> image;
                  applyMetric(fixedMetric, fixedProjection, vK.getData(), image);
                  rep.primalBarrierResidual = (image - fixedForce).norm();
                  rep.primalBarrierResidualTolerance = p.primalBarrierAbsoluteTolerance +
                    p.primalBarrierRelativeTolerance * residualScale;
                  rep.primalBarrierRelativeResidual = residualScale > Real(0)
                    ? rep.primalBarrierResidual / residualScale : rep.primalBarrierResidual;
                  innerConverged = rep.primalBarrierResidual <= rep.primalBarrierResidualTolerance;
                  rep.primalBarrierConverged = innerConverged;
                  rep.primalBarrierRelativeCorrection = Real(0);
                  if (p.trace)
                    std::cout << "        barrier skip: outer=" << rep.iterations
                              << "  reason=inactive-hinges  residual=" << rep.primalBarrierResidual
                              << "  rel=" << rep.primalBarrierRelativeResidual
                              << "  converged=" << innerConverged << '\n';
                  if (innerConverged)
                  {
                    ++rep.inactiveHingeSkips;
                    break;
                  }
                }
                Detail::WNGIRPrimalBarrierMetric barrierMetric(
                  m_duStep, m_vStep, u, vK, p, barrierCoefficient);
                Detail::WNGIRPrimalBarrierForce barrierForce(
                  m_vStep, u, vK, p, barrierCoefficient);
                typename ProblemType::ProblemBodyType body(m_regularityForm);
                body =
                  body + m_localMetricForm + barrierMetric - m_surfaceForm - barrierForce;
                m_stepProblem = body;
                assembleStep();
                const Real innerAssemblySeconds = secondsSince(tic);
                rep.tAssembly += innerAssemblySeconds;
                rep.tPrimalBarrierAssembly += innerAssemblySeconds;

                Math::Vector<Real> residualImage;
                const auto& innerSystem = m_stepProblem.getLinearSystem();
                applyMetric(innerSystem.getOperator(), fixedProjection, vK.getData(), residualImage);
                rep.primalBarrierResidual =
                  (residualImage - innerSystem.getVector()).norm();
                rep.primalBarrierResidualTolerance = p.primalBarrierAbsoluteTolerance +
                  p.primalBarrierRelativeTolerance * residualScale;
                rep.primalBarrierRelativeResidual = residualScale > Real(0)
                  ? rep.primalBarrierResidual / residualScale : rep.primalBarrierResidual;
                innerConverged = std::isfinite(rep.primalBarrierResidual) &&
                  rep.primalBarrierResidual <= rep.primalBarrierResidualTolerance;
                if (p.trace)
                {
                  const auto precision = std::cout.precision();
                  std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
                            << "        barrier residual: outer=" << rep.iterations
                            << "  inner=" << inner
                            << "  residual=" << rep.primalBarrierResidual
                            << "  force_norm=" << fixedForce.norm()
                            << "  rel=" << rep.primalBarrierRelativeResidual
                            << "  tolerance=" << rep.primalBarrierResidualTolerance
                            << "  converged=" << innerConverged << '\n';
                  std::cout.precision(precision);
                }
                if (innerConverged || inner == innerIterations)
                  break;

                tic = Clock::now();
                std::size_t barrierIterations = 0;
                Real barrierError = std::numeric_limits<Real>::infinity();
                solveOk =
                  solveStep(uTrial, barrierIterations, barrierError, &fixedProjection, &rep);
                recordLinearSolve(barrierIterations, barrierError);
                const Real innerSolveSeconds = secondsSince(tic);
                rep.tSolve += innerSolveSeconds;
                rep.tPrimalBarrierSolve += innerSolveSeconds;
                if (!solveOk)
                {
                  if (p.trace)
                    std::cout << "        barrier inner=" << (inner + 1)
                              << "  outer=" << rep.iterations << "  linear_ok=0"
                              << "  cg_it=" << barrierIterations
                              << "  cg_err=" << barrierError << std::endl;
                  rep.exitReason = "solve-linear-failed";
                  break;
                }

                scratch = uTrial;
                scratch -= vK;
                const Real correctionNorm =
                  std::max(std::abs(scratch.max()), std::abs(scratch.min()));
                const Real iterateNorm = std::max(std::abs(vK.max()), std::abs(vK.min()));
                const Real correctionScale = std::max(iterateNorm, predictorNorm);
                const Real relativeCorrection = correctionScale > Real(0)
                  ? correctionNorm / correctionScale
                  : (correctionNorm == Real(0) ? Real(0)
                                               : std::numeric_limits<Real>::infinity());
                rep.primalBarrierRelativeCorrection = relativeCorrection;
                ++rep.primalBarrierIterations;
                ++rep.lastPrimalBarrierIterations;
                rep.maxPrimalBarrierIterations =
                  std::max(rep.maxPrimalBarrierIterations, rep.lastPrimalBarrierIterations);
                Real innerAlpha = Real(1);
                if constexpr (std::is_same_v<typename FormLanguage::Traits<
                                               LinearSystemType>::OperatorType,
                                Math::SparseMatrix<Real>>)
                {
                  {
                    const auto meritStart = Clock::now();
                    const auto merit = [&](const Displacement& increment) {
                      Math::Vector<Real> image;
                      applyMetric(
                        fixedMetric, fixedProjection, increment.getData(), image);
                      const Real quadratic = Real(0.5) * increment.getData().dot(image);
                      const Real force = fixedForce.dot(increment.getData());
                      const Real quality = getPrimalBarrierEnergy(mesh, fes,
                        validationCells, u, increment, meshDim, barrierCoefficient);
                      return std::pair<Real, Real>{quadratic - force + quality,
                        std::abs(quadratic) + std::abs(force) + std::abs(quality)};
                    };
                    const auto before = merit(vK);
                    Math::Vector<Real> image;
                    const auto& system = m_stepProblem.getLinearSystem();
                    applyMetric(
                      system.getOperator(), fixedProjection, vK.getData(), image);
                    const Real slope =
                      (image - system.getVector()).dot(scratch.getData());
                    bool meritAccepted = false;
                    std::size_t backtracks = 0;
                    Real after = std::numeric_limits<Real>::infinity();
                    for (; backtracks < 32; ++backtracks)
                    {
                      uTrial = scratch;
                      uTrial *= innerAlpha;
                      uTrial += vK;
                      const auto trial = merit(uTrial);
                      after = trial.first;
                      const Real rounding = Real(64) *
                        std::numeric_limits<Real>::epsilon() *
                        std::max(before.second, trial.second);
                      if (std::isfinite(before.first) && std::isfinite(after) &&
                        std::isfinite(slope) &&
                        slope < Real(0) &&
                        after <= before.first +
                            Real(1e-4) * innerAlpha * std::min(Real(0), slope) + rounding)
                      {
                        meritAccepted = true;
                        break;
                      }
                      innerAlpha *= Real(0.5);
                    }
                    rep.primalBarrierBacktracks += backtracks;
                    rep.tPrimalBarrierLineSearch += secondsSince(meritStart);
                    if (p.trace)
                    {
                      const auto precision = std::cout.precision();
                      std::cout
                        << std::setprecision(std::numeric_limits<Real>::max_digits10)
                        << "        barrier merit: outer=" << rep.iterations
                        << "  inner=" << (inner + 1) << "  before=" << before.first
                        << "  after=" << after << "  slope=" << slope
                        << "  alpha=" << innerAlpha << "  backtracks=" << backtracks
                        << "  accepted=" << meritAccepted << std::endl;
                      std::cout.precision(precision);
                    }
                    if (!meritAccepted)
                    {
                      if (p.trace)
                        std::cout
                          << "        barrier inner=" << (inner + 1)
                          << "  outer=" << rep.iterations << "  linear_ok=1"
                          << "  merit_ok=0  converged=0  cg_it=" << barrierIterations
                          << "  cg_err=" << barrierError << std::endl;
                      rep.exitReason = "primal-barrier-line-search-failure";
                      solveOk = false;
                      break;
                    }
                  }
                }
                rep.lastPrimalBarrierAlpha = innerAlpha;
                rep.minPrimalBarrierAlpha =
                  std::min(rep.minPrimalBarrierAlpha, innerAlpha);
                if (innerAlpha >= Real(1) - Real(1e-12))
                  ++rep.fullPrimalBarrierSteps;
                if (p.trace)
                {
                  const auto precision = std::cout.precision();
                  std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
                            << "        barrier inner=" << (inner + 1)
                            << "  outer=" << rep.iterations << "  corr=" << correctionNorm
                            << "  iterate=" << iterateNorm
                            << "  rel=" << relativeCorrection << "  alpha=" << innerAlpha
                            << "  linear_ok=1  cg_it=" << barrierIterations
                            << "  cg_err=" << barrierError
                            << "  mu_eff=" << barrierCoefficient << std::endl;
                  std::cout.precision(precision);
                }
                scratch *= innerAlpha;
                vK += scratch;
                // Reassemble at the accepted inner iterate to certify its stationarity residual.
                tic = Clock::now();
              }
              rep.primalBarrierConverged = innerConverged;
              if (solveOk && !innerConverged)
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

          if (!(std::isfinite(vK.max()) && std::isfinite(vK.min())))
          {
            rep.exitReason = "solve-nonfinite";
            break;
          }

          Real directionAction = getSurfaceForceAction(mesh, fes, u, vK, phi, grad,
            interfaceFacets, loss, dataNormalization, meshDim, locator);
          const Real predictorNorm =
            std::max(std::abs(predictor.max()), std::abs(predictor.min()));
          Real directionNorm = std::max(std::abs(vK.max()), std::abs(vK.min()));
          rep.directionAction = directionAction;
          rep.descentRatio =
            predictorAction > Real(0) ? directionAction / predictorAction : Real(0);
          rep.directionNormRatio =
            predictorNorm > Real(0) ? directionNorm / predictorNorm : Real(0);
          if (!(std::isfinite(directionAction) && directionAction > Real(0)))
          {
            rep.exitReason = "no-descent-direction";
            break;
          }

          // ---- Nonlinear line search on TRUE geometry ----
          tic = Clock::now();
          Real alpha = Real(1);
          bool accepted = false;
          std::size_t backtracks = 0;
          FastAdm adm{};
          Real eTrial = std::numeric_limits<Real>::infinity();
          SurfaceState trialSurface{};
          {
            previousU = u;
            for (; backtracks <= p.maxBacktracks; ++backtracks)
            {
              // uTrial = previousU + alpha * vK
              uTrial = vK;
              uTrial *= alpha;
              uTrial += previousU;
              if ((uTrial.getData().array() == previousU.getData().array()).all())
                break;
              bool jOK = true;
              bool qOK = true;
              {
                adm = p.trace ? getAdmissibilityState(
                                  mesh, fes, validationCells, uTrial, meshDim, &previousU)
                              : fastAdmissibility(uTrial);
                jOK = adm.inadmissibleCount == 0 && adm.minJ > acceptedJacobianFloor;
                qOK = adm.maxQ < p.qMax;
              }
              bool eOK = true;
              if (jOK && qOK)
              {
                trialSurface = surfaceState(uTrial);
                eTrial = trialSurface.energy;
                const Real sufficientDecrease =
                  p.armijoCoefficient * alpha * directionAction;
                eOK = std::isfinite(eTrial) && eTrial <= ePrev - sufficientDecrease;
              }
              if (p.trace)
              {
                const auto precision = std::cout.precision();
                std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
                          << "      wngir trial: outer=" << rep.iterations
                          << "  alpha=" << alpha << "  affine_min_j=" << adm.affineMinJ
                          << "  affine_max_qrel=" << adm.affineMaxQ
                          << "  min_j=" << adm.minJ << "  max_qrel=" << adm.maxQ
                          << "  j_ok=" << jOK << "  q_ok=" << qOK
                          << "  energy_checked=" << (jOK && qOK)
                          << "  e_ok=" << (jOK && qOK && eOK) << std::endl;
                if (p.traceQualityWitness)
                  std::cout << "      wngir quality witness: outer=" << rep.iterations
                            << "  alpha=" << alpha << "  cell=" << adm.qualityCell
                            << "  current_q=" << adm.qualityCurrent
                            << "  actual_q=" << adm.maxQ
                            << "  linear_change=" << adm.qualityLinearChange
                            << "  quadratic_change=" << adm.qualityQuadraticChange
                            << "  remainder="
                            << adm.maxQ - adm.qualityCurrent - adm.qualityLinearChange
                            << "  current_margin=" << p.qMax - adm.qualityCurrent
                            << std::endl;
                std::cout.precision(precision);
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
            }
          }
          rep.tLineSearch += secondsSince(tic);
          if (!accepted)
          {
            u = previousU;
            rep.exitReason = "line-search-failure";
            break;
          }

          rep.lastAlpha = alpha;
          {
            // Measure the accepted FE field, not its coefficient vector.
            scratch = u;
            scratch -= previousU;
            rep.acceptedStep = Detail::wngirPhysicalDisplacementNorm(
              mesh, fes, validationCells, scratch, p.geometricValidationOrder);
          }
          rep.minJ = adm.minJ;
          rep.maxJ = adm.maxJ;
          rep.maxQRel = adm.maxQ;

          const auto& surf = trialSurface;
          currentSurface = surf;
          const Real eNow = eTrial;
          rep.backtracks += backtracks;
          rep.actualPredictedDecrease = alpha * directionAction > Real(0)
            ? (ePrev - eNow) / (alpha * directionAction)
            : Real(0);
          recordSurfaceState(surf);
          bool geometricTargetReached = false;
          {
            const auto geometry = getInterfaceGeometryState(
              mesh, fes, u, phi, grad, interfaceFacets, meshDim, locator);
            if (p.trace)
              traceGeometry(geometry, rep.iterations + 1, "accepted");
            if (!std::isfinite(geometry.sup))
            {
              rep.exitReason = "geometric-validation-failed";
              ++rep.iterations;
              break;
            }
            geometricTargetReached = geometry.sup <= geometricTarget;
          }
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
                      << "  predictor_action=" << rep.predictorAction
                      << "  inner_residual=" << rep.primalBarrierResidual
                      << "  inner_relative_residual=" << rep.primalBarrierRelativeResidual
                      << "  desc=" << rep.descentRatio
                      << "  dirNorm=" << rep.directionNormRatio
                      << "  ared/pred=" << rep.actualPredictedDecrease
                      << "  bt=" << backtracks << "  rejJ=" << rep.jacobianRejections
                      << "  rejQ=" << rep.distortionRejections
                      << "  rejE=" << rep.energyRejections << "  min_j=" << rep.minJ
                      << "  max_j=" << rep.maxJ << "  max_Q=" << rep.maxQRel << '\n';

          if (geometricTargetReached)
          {
            rep.exitReason = "full-interface-geometric-sup-converged";
            ++rep.iterations;
            break;
          }
          const Real stepThreshold =
            stepTol + h * p.acceptedStepOverHTol;
          consecutiveSmallAcceptedSteps = rep.acceptedStep <= stepThreshold
            ? consecutiveSmallAcceptedSteps + 1 : 0;
          const Real energyScale =
            std::max({std::abs(ePrev), std::abs(eNow), std::numeric_limits<Real>::min()});
          const Real energyChange = std::abs(ePrev - eNow) / energyScale;
          consecutiveSmallEnergyChanges = p.energyStagTol > Real(0) &&
              energyChange <= p.energyStagTol
            ? consecutiveSmallEnergyChanges + 1 : 0;
          if (consecutiveSmallAcceptedSteps >= p.stagnationIterations ||
            consecutiveSmallEnergyChanges >= p.stagnationIterations)
          {
            rep.exitReason = consecutiveSmallAcceptedSteps >= p.stagnationIterations
              ? "best-effort-step-stagnation" : "best-effort-energy-stagnation";
            ++rep.iterations;
            break;
          }
          ePrev = eNow;
        }

        if (p.rigidDiagnostics)
        {
          const RigidModeState finalRigid = getRigidModeState(mesh, fes, u, phi, grad,
            interfaceFacets, sigma2, dataNormalization, meshDim, locator);
          rep.rigidModeCoercivity = finalRigid.minimum;
          rep.rigidModeCoercivityRatio = finalRigid.ratio;
          rep.rigidModeDimension = finalRigid.dimension;
        }

        const InterfaceGeometryState geometry = getInterfaceGeometryState(
          mesh, fes, u, phi, grad, interfaceFacets, meshDim, locator);
        rep.geometricRMS = geometry.rms;
        rep.geometricSup = geometry.sup;
        rep.normalRMS = geometry.normalRMS;
        rep.qualityBudgetSatisfied = rep.minJ > acceptedJacobianFloor && rep.maxQRel < p.qMax;
        rep.geometricTargetReached = rep.qualityBudgetSatisfied && geometry.sup <= geometricTarget;
        if (p.trace)
          traceGeometry(geometry, rep.iterations, "final");

        m_report = rep;
        return rep;
      }

    private:
      /// @brief Geometric discrepancy of the complete fitted interface.
      struct InterfaceGeometryState
      {
          Real rms = std::numeric_limits<Real>::infinity();
          Real sup = std::numeric_limits<Real>::infinity();
          Real normalRMS = std::numeric_limits<Real>::infinity();
      };

      /// Robust fitting Hessian action, omitting D2 phi; no metric assembly/solve.
      template <class Mesh, class FES, class PhiType, class GradType, class LocatorType>
      std::pair<Real, Real> getSurfaceDirectionalCurvature(const Mesh& mesh,
        const FES& fes, const Displacement& current, const Displacement& direction,
        const PhiType& phi, const GradType& grad,
        const std::vector<Index>& interfaceFacets, const WNGIRLoss& loss,
        Real normalization, size_t dimension, const LocatorType& locator) const
      {
        if constexpr (requires { current.acquire(); })
          current.acquire();
        if constexpr (requires { direction.acquire(); })
          direction.acquire();
        if constexpr (requires { phi.acquire(); })
          phi.acquire();
        if constexpr (requires { grad.acquire(); })
          grad.acquire();
        std::vector<std::pair<Real, Real>> facetCurvatures(interfaceFacets.size());
#ifdef RODIN_USE_OPENMP
#pragma omp parallel
#endif
        {
          DeformationMap<Displacement, LocatorType> deformation(current, locator);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(interfaceFacets.size()); ++i)
          {
            const auto face = mesh.getFace(interfaceFacets[static_cast<size_t>(i)]);
            const auto& qf = getQuadrature(*face, fes);
            const auto& quadrature = face->getQuadrature(qf);
            for (size_t q = 0; q < quadrature.getSize(); ++q)
            {
              const auto& point = quadrature.getPoint(q);
              const Variational::IntegrationPoint ip(point, &qf, q);
              const auto value = direction.getValue(point);
              const Detail::WNGIRResidualState state(
                phi, grad, deformation, ip, loss, true);
              const Real residual = state.getResidual();
              const Real action = state.getGradient().dot(value);
              const Real weight =
                qf.getWeight(q) * point.getDistortion() * normalization * action * action;
              facetCurvatures[i].first += weight * loss.getCurvature(residual);
              facetCurvatures[i].second += weight * state.getWeight();
            }
          }
        }
        std::pair<Real, Real> curvature{};
        for (const auto& value : facetCurvatures)
        {
          curvature.first += value.first;
          curvature.second += value.second;
        }
        return curvature;
      }

      /**
       * @brief Quadrature formula WNGIR uses on a polytope.
       *
         * Cell sampling integrates finite-element products at order @f$2k@f$.
         * Interface sampling uses a higher minimum because its composed
         * level-set coefficients are non-polynomial. A pinned parameter overrides both.
       */
      template <class FES>
      const QF::QuadratureFormulaBase& getQuadrature(
        const Geometry::Polytope& polytope, const FES& fes) const
      {
        const auto& fe =
          fes.getFiniteElement(polytope.getDimension(), polytope.getIndex());
        const bool isInterface = polytope.getDimension() < fes.getMesh().getDimension();
        const std::size_t automaticOrder = isInterface
          ? wngirInterfaceQuadratureOrder(fe.getOrder())
          : std::max<std::size_t>(2, 2 * fe.getOrder());
        const std::size_t order = m_parameters.quadratureOrder > 0
          ? m_parameters.quadratureOrder
          : automaticOrder;
        return QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
      }

      /// @brief Higher-order quadrature used only for reported geometric responses.
      template <class FES>
      const QF::QuadratureFormulaBase& getGeometricValidationQuadrature(
        const Geometry::Polytope& polytope, const FES& fes) const
      {
        const auto& fe =
          fes.getFiniteElement(polytope.getDimension(), polytope.getIndex());
        const std::size_t order = m_parameters.geometricValidationOrder > 0
          ? m_parameters.geometricValidationOrder
          : wngirGeometricValidationOrder(fe.getOrder());
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
        const GradType& grad, const std::vector<Index>& interfaceFacets,
        const WNGIRLoss& loss, Real normalization, std::size_t dimension,
        const LocatorType& locator) const
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
            phi, grad, current, locator, loss, normalization, dimension);
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
       * from the analytical model. No essential boundary condition is imposed,
       * so the rigid kernel is always present and the diagnostic is never
       * vacuous.
       */
      template <class Mesh, class FES, class PhiType, class GradType, class LocatorType>
      RigidModeState getRigidModeState(const Mesh& mesh, const FES& fes,
        const Displacement& current, const PhiType& phi, const GradType& grad,
        const std::vector<Index>& interfaceFacets, Real sigma2, Real normalization,
        std::size_t dimension, const LocatorType& locator) const
      {
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
            gradients[mode](static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) =
              Real(-1);
            gradients[mode](static_cast<Eigen::Index>(b), static_cast<Eigen::Index>(a)) =
              Real(1);
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
                  weight *
                  (values[i].dot(values[j]) +
                    (gradients[i].array() * gradients[j].array()).sum());
              }
            }
          }
        }

        Math::Matrix<Real> observationGram = Math::Matrix<Real>::Zero(count, count);
        Detail::WNGIRObservationCoefficient observation(
          grad, current, locator, m_parameters, normalization, dimension);
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

      /// @brief Integrated affine quality energy on the same quadrature as its force and Hessian.
      template <class Mesh, class FES>
      bool hasActivePrimalHinges(const Mesh& mesh, const FES& fes,
        const std::vector<Index>& cells, const Displacement& current,
        const Displacement& inner, std::size_t dimension, Real coefficient) const
      {
        bool active = false;
#ifdef RODIN_USE_OPENMP
#pragma omp parallel reduction(|| : active)
#endif
        {
          auto currentJacobian = Variational::Jacobian(current);
          auto innerJacobian = Variational::Jacobian(inner);
          CellDeformation deformation(dimension);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(cells.size()); ++i)
          {
            const auto cell = mesh.getCell(cells[static_cast<std::size_t>(i)]);
            const auto& qf = getQuadrature(*cell, fes);
            const auto& quadrature = cell->getQuadrature(qf);
            for (std::size_t q = 0; q < quadrature.getSize(); ++q)
            {
              const Variational::IntegrationPoint ip(quadrature.getPoint(q), &qf, q);
              deformation.setDisplacementGradient(currentJacobian.getValue(ip));
              const Detail::WNGIRPrimalBarrierState state(
                deformation, innerJacobian.getValue(ip), m_parameters, coefficient);
              active = active || !state.isFeasible() ||
                state.getJacobianHessian() != Real(0) ||
                state.getDistortionHessian() != Real(0);
            }
          }
        }
        return active;
      }

      template <class Mesh, class FES>
      Real getPrimalBarrierEnergy(const Mesh& mesh, const FES& fes,
        const std::vector<Index>& validationCells, const Displacement& current,
        const Displacement& inner, std::size_t dimension, Real coefficient) const
      {
        Real energy = 0;
#ifdef RODIN_USE_OPENMP
#pragma omp parallel reduction(+ : energy)
#endif
        {
          auto currentJacobian = Variational::Jacobian(current);
          auto innerJacobian = Variational::Jacobian(inner);
          CellDeformation deformation(dimension);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(validationCells.size()); ++i)
          {
            const auto cell = mesh.getCell(validationCells[static_cast<std::size_t>(i)]);
            const auto& qf = getQuadrature(*cell, fes);
            const auto& quadrature = cell->getQuadrature(qf);
            for (std::size_t q = 0; q < quadrature.getSize(); ++q)
            {
              const auto& point = quadrature.getPoint(q);
              const Variational::IntegrationPoint ip(point, &qf, q);
              deformation.setDisplacementGradient(currentJacobian.getValue(ip));
              const Detail::WNGIRPrimalBarrierState state(
                deformation, innerJacobian.getValue(ip), m_parameters, coefficient);
              energy += qf.getWeight(q) * point.getDistortion() *
                state.getEnergy(m_parameters, coefficient);
            }
          }
        }
        return energy;
      }

      template <class Mesh, class FES>
      AdmissibilityState getAdmissibilityState(const Mesh& mesh, const FES& fes,
        const std::vector<Index>& validationCells, const Displacement& u,
        std::size_t dimension, const Displacement* current = nullptr) const
      {
        if constexpr (requires { u.acquire(); })
          u.acquire();
        Real minJ = std::numeric_limits<Real>::infinity();
        Real maxJ = -std::numeric_limits<Real>::infinity();
        Real maxQ = Real(0);
        Real affineMinJ = std::numeric_limits<Real>::infinity();
        Real affineMaxQ = Real(0);
        if (current)
        {
          if constexpr (requires { current->acquire(); })
            current->acquire();
        }
        std::size_t inadmissibleCount = 0;
        struct Witness
        {
            Real actual = 0, current = 0, linear = 0, quadratic = 0;
            Index cell = 0;
        };
        const bool recordWitness = current && m_parameters.traceQualityWitness;
        std::vector<Witness> witnesses(recordWitness ? validationCells.size() : 0);
#ifdef RODIN_USE_OPENMP
#pragma omp parallel reduction(min : minJ, affineMinJ)                                   \
  reduction(max : maxJ, maxQ, affineMaxQ) reduction(+ : inadmissibleCount)
#endif
        {
          CellDeformation deformation(dimension);
          CellDeformation currentDeformation(dimension);
          auto displacementJacobian = Variational::Jacobian(u);
          auto currentJacobian = Variational::Jacobian(current ? *current : u);
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
              if (!std::isfinite(j) || j <= m_parameters.jMinRatio || !deformation.isAdmissible())
                ++inadmissibleCount;
              if (deformation.isAdmissible())
              {
                const Real quality = deformation.getRelativeDistortion();
                if (!std::isfinite(quality))
                  ++inadmissibleCount;
                else
                  maxQ = std::max(maxQ, quality);
              }
              if (current)
              {
                const auto gradient = currentJacobian.getValue(ip);
                const SpatialMat increment(displacementJacobian.getValue(ip) - gradient);
                currentDeformation.setDisplacementGradient(gradient);
                affineMinJ = std::min(affineMinJ,
                  currentDeformation.getJacobian() +
                    currentDeformation.getJacobianAction(increment));
                affineMaxQ = std::max(affineMaxQ,
                  currentDeformation.getRelativeDistortion() +
                    currentDeformation.getRelativeDistortionAction(increment));
                if (recordWitness && deformation.isAdmissible())
                {
                  auto& witness = witnesses[static_cast<size_t>(i)];
                  const Real actual = deformation.getRelativeDistortion();
                  if (actual > witness.actual)
                    witness = {actual, currentDeformation.getRelativeDistortion(),
                      currentDeformation.getRelativeDistortionAction(increment),
                      Real(0.5) *
                        currentDeformation.getRelativeDistortionSecondAction(
                          increment, increment),
                      cellIndex};
                }
              }
            }
          }
        }
        AdmissibilityState result{
          minJ, maxJ, maxQ, inadmissibleCount, affineMinJ, affineMaxQ};
        Witness limiting;
        for (const auto& witness : witnesses)
          if (witness.actual > limiting.actual)
            limiting = witness;
        result.qualityCell = limiting.cell;
        result.qualityCurrent = limiting.current;
        result.qualityLinearChange = limiting.linear;
        result.qualityQuadraticChange = limiting.quadratic;
        return result;
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
       * @brief Evaluates position and normal errors on the fitted interface.
       *
       * At a mapped interface point @f$y=T_h(x)@f$, the normalized residual
       * @f$\phi(y)/|\nabla\phi(y)|@f$ is the first-order signed distance to the
       * target zero set. The fitted normal is obtained from the tangents mapped
       * by @f$I+\nabla u_h@f$; consequently, curved and higher-order displacement
       * fields are measured without replacing their image by an affine facet.
       * Normal orientation is ignored because interior-facet orientation is not
       * intrinsic to the interface skeleton.
       */
      template <class Mesh, class FES, class PhiType, class GradType, class LocatorType>
      InterfaceGeometryState getInterfaceGeometryState(const Mesh& mesh, const FES& fes,
        const Displacement& u, const PhiType& phi, const GradType& grad,
        const std::vector<Index>& interfaceFacets, std::size_t dimension,
        const LocatorType& locator) const
      {
        if constexpr (requires { u.acquire(); })
          u.acquire();
        if constexpr (requires { phi.acquire(); })
          phi.acquire();
        if constexpr (requires { grad.acquire(); })
          grad.acquire();

        struct Accumulation
        {
            Real measure = 0;
            Real squaredDistance = 0;
            Real maximumDistance = 0;
            Real squaredNormal = 0;
            std::size_t invalidCount = 0;
        };
        std::vector<Accumulation> facetStates(interfaceFacets.size());
        constexpr Real gradientFloor = Real(1e-14);

#ifdef RODIN_USE_OPENMP
#pragma omp parallel
#endif
        {
          DeformationMap deformation(u, locator);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(interfaceFacets.size()); ++i)
          {
            Accumulation& state = facetStates[static_cast<std::size_t>(i)];
            const Index facetIndex = interfaceFacets[static_cast<std::size_t>(i)];
            const auto face = mesh.getFace(facetIndex);
            const auto& fe = fes.getFiniteElement(face->getDimension(), facetIndex);
            const auto& qf = getGeometricValidationQuadrature(*face, fes);
            const auto& quadrature = face->getQuadrature(qf);
            // Endpoints/vertices contribute to the maximum, not the surface quadrature integral.
            const Geometry::Polytope::Traits traits(face->getGeometry());
            for (std::size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
            {
              const Geometry::Point point(*face, traits.getVertex(vertex));
              const auto& moved = deformation.getMovedPoint(Variational::IntegrationPoint(point));
              const Real norm = grad.getValue(moved).norm();
              const Real distance = std::abs(phi.getValue(moved)) / norm;
              if (!(norm > gradientFloor) || !std::isfinite(norm) || !std::isfinite(distance))
                ++state.invalidCount;
              else
                state.maximumDistance = std::max(state.maximumDistance, distance);
            }
            for (std::size_t q = 0; q < quadrature.getSize(); ++q)
            {
              const auto& point = quadrature.getPoint(q);
              const Variational::IntegrationPoint ip(point, &qf, q);
              const auto& J = point.getJacobian();
              SpatialMat mappedJacobian = J;
              const auto& rc = point.getReferenceCoordinates();
              for (std::size_t local = 0; local < fe.getCount(); ++local)
              {
                const Index dof =
                  fes.getGlobalIndex({face->getDimension(), facetIndex}, local);
                for (std::size_t axis = 0; axis < face->getDimension(); ++axis)
                {
                  for (std::size_t component = 0; component < dimension; ++component)
                  {
                    mappedJacobian(static_cast<Eigen::Index>(component),
                      static_cast<Eigen::Index>(axis)) += u[dof] *
                      fe.getBasis(local).template getDerivative<1>(component, axis)(rc);
                  }
                }
              }

              SpatialVec fittedNormal = SpatialVec::Zero(dimension);
              Real mappedMeasure = 0;
              if (dimension == 1)
              {
                mappedMeasure = Real(1);
                fittedNormal(0) = Real(1);
              }
              else if (dimension == 2)
              {
                SpatialVec tangent = SpatialVec::Zero(2);
                tangent(0) = mappedJacobian(0, 0);
                tangent(1) = mappedJacobian(1, 0);
                mappedMeasure = tangent.norm();
                fittedNormal(0) = tangent(1);
                fittedNormal(1) = -tangent(0);
              }
              else
              {
                SpatialVec a = SpatialVec::Zero(3);
                SpatialVec b = SpatialVec::Zero(3);
                for (std::size_t r = 0; r < 3; ++r)
                {
                  a(static_cast<Eigen::Index>(r)) = mappedJacobian(r, 0);
                  b(static_cast<Eigen::Index>(r)) = mappedJacobian(r, 1);
                }
                fittedNormal = a.cross(b);
                mappedMeasure = fittedNormal.norm();
              }

              const Geometry::Point movedPoint = deformation.getMovedPoint(ip);
              const SpatialVec targetGradient = grad.getValue(movedPoint);
              const Real gradientNorm = targetGradient.norm();
              const Real fittedNormalNorm = fittedNormal.norm();
              if (!(mappedMeasure > gradientFloor && gradientNorm > gradientFloor &&
                    fittedNormalNorm > gradientFloor) || !std::isfinite(mappedMeasure) ||
                  !std::isfinite(gradientNorm) || !std::isfinite(fittedNormalNorm))
              {
                ++state.invalidCount;
                continue;
              }

              fittedNormal /= fittedNormalNorm;
              const SpatialVec targetNormal = targetGradient / gradientNorm;
              const Real distance = std::abs(phi.getValue(movedPoint)) / gradientNorm;
              if (!std::isfinite(distance))
              {
                ++state.invalidCount;
                continue;
              }
              const Real normalErrorSquared = Real(2) -
                Real(2) *
                  std::clamp(std::abs(fittedNormal.dot(targetNormal)), Real(0), Real(1));
              const Real weight = qf.getWeight(q) * mappedMeasure;
              state.measure += weight;
              state.squaredDistance += weight * distance * distance;
              state.maximumDistance = std::max(state.maximumDistance, distance);
              state.squaredNormal += weight * normalErrorSquared;
            }
          }
        }

        Accumulation total;
        for (const Accumulation& state : facetStates)
        {
          total.measure += state.measure;
          total.squaredDistance += state.squaredDistance;
          total.maximumDistance = std::max(total.maximumDistance, state.maximumDistance);
          total.squaredNormal += state.squaredNormal;
          total.invalidCount += state.invalidCount;
        }
        if (!(total.measure > Real(0)) ||
          total.invalidCount > 0)
          return {};
        return {std::sqrt(std::max(Real(0), total.squaredDistance) / total.measure),
          total.maximumDistance,
          std::sqrt(std::max(Real(0), total.squaredNormal) / total.measure)};
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

      struct RegularityProjection
      {
          std::vector<Math::Vector<Real>> modes;
          std::vector<Real> weights;
      };

      /// @brief Assemble one symmetric operator for the model, factors and residuals.
      void assembleStep()
      {
        m_stepProblem.assemble();
        using OperatorType =
          typename FormLanguage::Traits<LinearSystemType>::OperatorType;
        if constexpr (std::is_same_v<OperatorType, Math::SparseMatrix<Real>>)
        {
          auto& matrix = m_stepProblem.getLinearSystem().getOperator();
          matrix =
            (Real(0.5) * (matrix + Math::SparseMatrix<Real>(matrix.transpose()))).eval();
          matrix.makeCompressed();
        }
      }

      /**
       * @brief Coefficient-orthonormal current translations, rotations and dilation.
       * Field interpolation respects the finite element's DOF functionals.
       */
      Math::Matrix<Real> getSimilarityModes(const Displacement& current) const
      {
        const auto& fes = current.getFiniteElementSpace();
        const size_t d = fes.getMesh().getDimension();
        const size_t count = d * (d + 1) / 2 + 1;
        Math::Matrix<Real> modes(current.getData().size(), count);
        Displacement mode(fes);
        size_t k = 0;
        for (size_t axis = 0; axis < d; ++axis)
        {
          mode = Variational::VectorFunction(d, [=](const Geometry::Point&) {
            SpatialVec value = SpatialVec::Zero(d);
            value(axis) = Real(1);
            return value;
          });
          modes.col(k++) = mode.getData();
        }
        for (size_t a = 0; a < d; ++a)
          for (size_t b = a + 1; b < d; ++b)
          {
            mode =
              Variational::VectorFunction(d, [&, a, b](const Geometry::Point& point) {
                const SpatialVec position(
                  point.getCoordinates() + current.getValue(point));
                SpatialVec value = SpatialVec::Zero(d);
                value(a) = -position(b);
                value(b) = position(a);
                return value;
              });
            modes.col(k++) = mode.getData();
          }
        mode = Variational::VectorFunction(d, [&](const Geometry::Point& point) {
          return SpatialVec(point.getCoordinates() + current.getValue(point));
        });
        modes.col(k) = mode.getData();
        for (size_t column = 0; column < count; ++column)
        {
          for (size_t pass = 0; pass < 2; ++pass)
            for (size_t previous = 0; previous < column; ++previous)
              modes.col(column) -=
                modes.col(previous).dot(modes.col(column)) * modes.col(previous);
          const Real norm = modes.col(column).norm();
          if (!(norm > Real(0)) || !std::isfinite(norm))
            return {};
          modes.col(column) /= norm;
        }
        return modes;
      }

      /**
       * @brief Selects zero in unresolved similarity modes of one linear system.
       *
       * Gauges are confined to the solve, not the metric or hinge objective.
       * They are recomputed for each active-set matrix, so a hinge that resolves
       * a previously free dilation releases that mode. True residuals are checked
       * against the ungauged operator. No large sparse eigenproblem is required.
       */
      size_t gaugeSimilarityModes(const Math::SparseMatrix<Real>& matrix,
        const Math::Vector<Real>& force, RegularityProjection& projection) const
      {
        const auto& modes = m_similarityModes;
        const size_t count = modes.cols();
        if (count == 0)
          return 0;
        Math::Matrix<Real> images(matrix.rows(), count);
        for (size_t column = 0; column < count; ++column)
        {
          Math::Vector<Real> image;
          applyMetric(matrix, projection, modes.col(column), image);
          images.col(column) = image;
        }
        Math::Matrix<Real> gram = modes.transpose() * images;
        gram = (Real(0.5) * (gram + gram.transpose())).eval();
        Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eigen(gram);
        if (eigen.info() != Eigen::Success)
          return 0;
        const Real scale = matrix.diagonal().cwiseAbs().maxCoeff();
        const Real tolerance = Real(256) * std::numeric_limits<Real>::epsilon();
        size_t unresolved = 0;
        for (size_t column = 0; column < count; ++column)
        {
          const Math::Vector<Real> candidate = modes * eigen.eigenvectors().col(column);
          if ((images * eigen.eigenvectors().col(column)).norm() <= tolerance * scale &&
            std::abs(force.dot(candidate)) <= tolerance * force.norm())
          {
            projection.modes.push_back(candidate);
            projection.weights.push_back(scale);
            ++unresolved;
          }
        }
        return unresolved;
      }

      void applyMetric(const Math::SparseMatrix<Real>& A,
        const RegularityProjection& projection, const Math::Vector<Real>& x,
        Math::Vector<Real>& y) const
      {
        y = A * x;
        for (std::size_t k = 0; k < projection.weights.size(); ++k)
          y += projection.weights[k] * projection.modes[k].dot(x) * projection.modes[k];
      }

      bool metricConjugateGradient(const Math::SparseMatrix<Real>& A,
        const RegularityProjection& projection, const Math::Vector<Real>& b,
        Math::Vector<Real>& x, std::size_t maxIterations, Real relativeTolerance,
        std::size_t& iterations, Real& error) const
      {
        Math::Vector<Real> jacobi(A.rows());
        for (Eigen::Index i = 0; i < A.rows(); ++i)
        {
          const Real d = A.coeff(i, i);
          jacobi(i) = (std::abs(d) > Real(0)) ? Real(1) / d : Real(1);
        }
        auto applyPreconditioner = [&](const Math::Vector<Real>& residual) {
          return Math::Vector<Real>(jacobi.cwiseProduct(residual));
        };
        const Real rhsNorm = b.norm();
        iterations = 0;
        if (!(rhsNorm > Real(0)))
        {
          x.setZero();
          error = Real(0);
          return true;
        }
        Math::Vector<Real> Ax;
        applyMetric(A, projection, x, Ax);
        Math::Vector<Real> r = b - Ax;
        Math::Vector<Real> z = applyPreconditioner(r);
        Math::Vector<Real> p = z;
        Real rz = r.dot(z);
        const Real threshold = relativeTolerance * rhsNorm;
        Math::Vector<Real> Ap;
        for (std::size_t it = 0; it < maxIterations; ++it)
        {
          if (r.norm() <= threshold)
            break;
          applyMetric(A, projection, p, Ap);
          const Real pAp = p.dot(Ap);
          if (!(pAp > Real(0)))
            break;
          const Real alpha = rz / pAp;
          x += alpha * p;
          r -= alpha * Ap;
          z = applyPreconditioner(r);
          const Real rzNext = r.dot(z);
          p = z + (rzNext / rz) * p;
          rz = rzNext;
          ++iterations;
        }
        error = r.norm() / rhsNorm;
        return std::isfinite(error) && error <= relativeTolerance;
      }

      // Solves the currently-assembled step problem and copies the
      // backend-matched solution GridFunction into @p out.
      bool solveStep(Displacement& out, std::size_t& iterations, Real& error,
        const RegularityProjection* fixedProjection = nullptr,
        WNGIRReport* report = nullptr)
      {
        auto& axb = m_stepProblem.getLinearSystem();
        using OperatorType =
          typename FormLanguage::Traits<LinearSystemType>::OperatorType;
        using VectorType = typename FormLanguage::Traits<LinearSystemType>::VectorType;
        if constexpr (std::is_same_v<OperatorType, Math::SparseMatrix<Real>> &&
          std::is_same_v<VectorType, Math::Vector<Real>>)
        {
          const auto& rhs = axb.getVector();
          const auto projection =
            fixedProjection ? *fixedProjection : RegularityProjection{};
          auto gauged = projection;
          const size_t unresolved = gaugeSimilarityModes(axb.getOperator(), rhs, gauged);
          if (report)
            report->unresolvedSimilarityModes = unresolved;
          if (m_parameters.directSolver != WNGIRParameters::DirectSolver::CG)
          {
            const auto solveDirect = [&](auto& direct, std::string_view backend) {
              const auto solveSystem = [&](LinearSystemType& system) {
                if constexpr (requires { direct.factorize(system); })
                {
                  auto& next = system.getOperator();
                  next.makeCompressed();
                  const auto& previous = m_directSystem.getOperator();
                  const bool samePattern = next.rows() == previous.rows() &&
                    next.cols() == previous.cols() &&
                    next.nonZeros() == previous.nonZeros() &&
                    std::equal(next.outerIndexPtr(),
                      next.outerIndexPtr() + next.outerSize() + 1,
                      previous.outerIndexPtr()) &&
                    std::equal(next.innerIndexPtr(),
                      next.innerIndexPtr() + next.nonZeros(), previous.innerIndexPtr());
                  const bool sameValues = samePattern &&
                    std::equal(next.valuePtr(), next.valuePtr() + next.nonZeros(),
                      previous.valuePtr());
                  if (!samePattern)
                    direct.clear(Solver::Factorization::Symbolic);
                  const bool numeric = direct.success() &&
                    direct.getInfo().factorization == Solver::Factorization::Numeric;
                  m_directSystem.getOperator() = next;
                  m_directSystem.getVector() = system.getVector();
                  if (!sameValues || !numeric)
                  {
                    if (report)
                    {
                      ++report->directFactorizations;
                      if (!direct.getInfo().factorization)
                        ++report->directAnalyses;
                    }
                    direct.factorize(m_directSystem);
                  }
                  if (direct.success())
                    direct.solve(m_directSystem);
                  if (direct.success())
                    system.getSolution() = m_directSystem.getSolution();
                }
                else
                  direct.solve(system);
              };
              const auto& matrix = axb.getOperator();
              const auto& rigid = gauged;
              if (rigid.weights.empty())
                solveSystem(axb);
              else
              {
                // Auxiliary diagonal signs also support negative low-rank regularity corrections.
                const auto n = matrix.rows();
                const auto rank = static_cast<Eigen::Index>(rigid.weights.size());
                std::vector<Math::SparseTriplet<Real>> entries;
                entries.reserve(matrix.nonZeros() + 2 * n * rank + rank);
                for (Eigen::Index column = 0; column < matrix.outerSize(); ++column)
                  for (Math::SparseMatrix<Real>::InnerIterator entry(matrix, column);
                    entry; ++entry)
                    entries.emplace_back(entry.row(), entry.col(), entry.value());
                for (Eigen::Index k = 0; k < rank; ++k)
                {
                  const Real scale = std::sqrt(std::abs(rigid.weights[k]));
                  for (Eigen::Index i = 0; i < n; ++i)
                  {
                    const Real value = scale * rigid.modes[k](i);
                    entries.emplace_back(i, n + k, value);
                    entries.emplace_back(n + k, i, value);
                  }
                  entries.emplace_back(
                    n + k, n + k, rigid.weights[k] > Real(0) ? Real(-1) : Real(1));
                }
                LinearSystemType augmented;
                augmented.getOperator().resize(n + rank, n + rank);
                augmented.getOperator().setFromTriplets(entries.begin(), entries.end());
                augmented.getVector() = Math::Vector<Real>::Zero(n + rank);
                augmented.getVector().head(n) = rhs;
                solveSystem(augmented);
                if (direct.success())
                  axb.getSolution() = augmented.getSolution().head(n);
              }
              iterations = 0;
              Math::Vector<Real> image;
              if (direct.success())
                applyMetric(matrix, projection, axb.getSolution(), image);
              error = direct.success() ? (image - rhs).norm() /
                  std::max(rhs.norm(), std::numeric_limits<Real>::min())
                                       : std::numeric_limits<Real>::infinity();
              const bool ok = direct.success() && axb.getSolution().allFinite() &&
                error <= m_parameters.cgRelativeTolerance;
              if (m_parameters.trace)
              {
                std::cout << "        metric " << backend << ": ok=" << ok
                          << "  rel=" << error << "  status=" << direct.getInfo().status;
                if constexpr (requires { direct.getResources().instance.infog[1]; })
                  std::cout << "  detail=" << direct.getResources().instance.infog[1];
                std::cout << '\n';
              }
              if (!ok)
                return false;
              m_duStep.getSolution().setData(axb.getSolution());
              out = m_duStep.getSolution();
              return true;
            };
#ifdef RODIN_USE_MUMPS
            if (m_parameters.directSolver == WNGIRParameters::DirectSolver::MUMPS)
            {
              if (!m_mumps)
                m_mumps = std::make_unique<Solver::MUMPS<LinearSystemType>>(m_stepProblem);
              m_mumps->setSymmetric(Solver::MUMPS<LinearSystemType>::Symmetry::General)
                .setMaxThreads(m_parameters.directSolverThreads);
              return solveDirect(*m_mumps, "MUMPS");
            }
#endif
            Solver::SparseLU direct(m_stepProblem);
            return solveDirect(direct, "LU");
          }
          const auto& guess = axb.getSolution();
          const std::size_t maxIterations = m_parameters.cgMaxIterations > 0
            ? m_parameters.cgMaxIterations
            : std::min<std::size_t>(2000,
                std::max<std::size_t>(
                  100, 2 * m_duStep.getFiniteElementSpace().getSize()));
          Math::Vector<Real> solution =
            (guess.size() == rhs.size()) ? guess : Math::Vector<Real>::Zero(rhs.size());
          const bool solved = metricConjugateGradient(axb.getOperator(), gauged, rhs,
            solution, maxIterations, m_parameters.cgRelativeTolerance, iterations, error);
          Math::Vector<Real> image;
          applyMetric(axb.getOperator(), projection, solution, image);
          error =
            (image - rhs).norm() / std::max(rhs.norm(), std::numeric_limits<Real>::min());
          const bool ok =
            solved && std::isfinite(error) && error <= m_parameters.cgRelativeTolerance;
          axb.getSolution() = solution;
          m_duStep.getSolution().setData(solution);
          out = m_duStep.getSolution();
          if (m_parameters.trace &&
            (!ok || !std::isfinite(out.max()) || !std::isfinite(out.min())))
            std::cout << "        cg failure: ok=" << ok << "  it=" << iterations
                      << "  max_it=" << maxIterations << "  err=" << error << "  finite="
                      << (std::isfinite(out.max()) && std::isfinite(out.min())) << '\n';
          return ok && std::isfinite(out.max()) && std::isfinite(out.min());
        }
        else
        {
          Alert::Exception() << "WNGIR requires the local Eigen backend." << Alert::Raise;
          return false;
        }
      }

      Displacement* m_u;
      TrialFunctionType m_duStep;
      TestFunctionType m_vStep;
      ProblemType m_stepProblem;
      LinearSystemType m_directSystem;
#ifdef RODIN_USE_MUMPS
      std::unique_ptr<Solver::MUMPS<LinearSystemType>> m_mumps;
#endif
      BilinearFormType m_regularityForm;
      Math::Matrix<Real> m_similarityModes;
      /// @brief Observation metric and fitting force at the outer displacement.
      ///
      /// Both depend on the outer displacement only, so they are assembled once
      /// per nonlinear iteration and reused by every barrier correction.
      BilinearFormType m_localMetricForm;
      LinearFormType m_surfaceForm;
      WNGIRParameters m_parameters;
      WNGIRReport m_report;
  };
}

#endif
