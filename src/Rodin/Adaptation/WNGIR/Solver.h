/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIR_SOLVER_H
#define RODIN_ADAPTATION_WNGIR_SOLVER_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "Rodin/Solver/NewtonSolver.h"
#include "HingeProblem.h"
#include "LinearSolver.h"
#include "Rodin/Variational/Integrator.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/Trace.h"
#include "Rodin/Variational/IdentityMatrix.h"
#include "Rodin/Variational/VectorFunction.h"

#include "Rodin/Location.h"

#include "../CellDeformation.h"
#include "Loss.h"
#include "Parameters.h"
#include "DirectionalNewton.h"
#include "Distribution.h"
#include "HingeForce.h"
#include "HingeMetric.h"
#include "Report.h"
#include "FittingCoefficient.h"
#include "FittingForce.h"

namespace Rodin::Adaptation
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
        maximum =
          std::max(maximum, evaluate(Geometry::Point(*cell, traits.getVertex(vertex))));
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
   *
   * @par Architecture
   * The robust fitting energy supplies the force and outer Armijo merit.
   * Fitting and current-strain forms define the frozen metric. Optional user
   * bilinear integrators augment that metric, not the fitting energy.
   * Homogeneous Dirichlet conditions constrain every predictor and inner
   * increment through native DOF assembly and linear-system elimination.
   * The affine quadratic hinges correct the predictor before actual quality
   * and energy acceptance. Solve-only similarity gauges resolve compatible
   * null modes without adding a mass penalty.
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
      using Monitor = std::function<void(const WNGIRReport&)>;

      /// @brief Observes accepted steps and the final report, without changing the solve.
      WNGIR& setMonitor(Optional<Monitor> monitor)
      {
        m_monitor = std::move(monitor);
        return *this;
      }

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
          m_trialUUID(du.getUUID()),
          m_duStep(du.getFiniteElementSpace()),
          m_vStep(v.getFiniteElementSpace()),
          m_metricProblem(m_duStep, m_vStep),
          m_hingeProblem(m_duStep, m_vStep),
          m_linearSolver(m_hingeProblem),
          m_distributionForm(m_duStep, m_vStep),
          m_fittingMetric(m_duStep, m_vStep),
          m_fittingForce(m_vStep),
          m_additionalMetric(du, v)
      {}

      WNGIR(const WNGIR&) = delete;
      WNGIR& operator=(const WNGIR&) = delete;

      /// @brief Additional metric terms, assembled afresh at each outer iteration.
      /// These augment F + D; they do not change the fitting energy or force.
      BilinearFormType& getMetric()
      {
        return m_additionalMetric;
      }

      const BilinearFormType& getMetric() const
      {
        return m_additionalMetric;
      }

      /// @brief Selects the marked interface to fit on the displacement mesh.
      WNGIR& setInterfaceAttribute(Geometry::Attribute attribute)
      {
        m_parameters.interfaceAttribute = attribute;
        return *this;
      }

      /**
       * @brief Holds the current boundary displacement fixed during fitting.
       *
       * Only value-prescribing, homogeneous conditions on the constructor's
       * trial function are supported. They constrain increments, not accumulated
       * displacement. Nonzero values and identification constraints are rejected.
       */
      WNGIR& operator+=(const Variational::DirichletBCBase<Real>& condition)
      {
        if (condition.getOperand().getUUID() != m_trialUUID)
          Alert::Exception() << "WNGIR boundary conditions must use its trial function."
                             << Alert::Raise;
        m_boundaryConditions.add(condition);
        return *this;
      }

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
        if (parameters.linearSolver == WNGIRParameters::LinearSolver::MUMPS)
          Alert::Exception() << "WNGIR MUMPS solves require RODIN_USE_MUMPS."
                             << Alert::Raise;
#endif
        if (!std::isfinite(parameters.kappaD) || !(parameters.kappaD > Real(0)) ||
          !std::isfinite(parameters.kappaF) || !(parameters.kappaF > Real(0)))
          Alert::Exception() << "WNGIR requires positive "
                                "fitting/distribution weights."
                             << Alert::Raise;
        for (const Real value :
          {parameters.kappaJ, parameters.kappaQ, parameters.innerAbsoluteTolerance,
            parameters.energyStagTol, parameters.stepTol, parameters.acceptedStepOverHTol,
            parameters.robustScale})
          if (!std::isfinite(value) || value < Real(0))
            Alert::Exception()
              << "WNGIR weights and tolerances must be finite and nonnegative."
              << Alert::Raise;
        if (!std::isfinite(parameters.innerRelativeTolerance) ||
          !(parameters.innerRelativeTolerance > Real(0)) ||
          !std::isfinite(parameters.linearRelativeTolerance) ||
          !(parameters.linearRelativeTolerance > Real(0)) ||
          !(parameters.armijoCoefficient > Real(0) &&
            parameters.armijoCoefficient < Real(1)) ||
          !(parameters.omegaMin > Real(0) && parameters.omegaMin <= Real(1)) ||
          !(parameters.jMinRatio > Real(0) && parameters.jMinRatio < Real(1)) ||
          !(parameters.jLineSearchRatio > Real(0) &&
            parameters.jLineSearchRatio < Real(1)) ||
          parameters.innerIterations == 0 || parameters.stagnationIterations == 0)
          Alert::Exception()
            << "WNGIR requires valid residual tolerances, budgets and line search."
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
      template <class PhiDerived, class GradDerived>
      WNGIRReport solve(const Variational::RealFunctionBase<PhiDerived>& phi,
        const Variational::VectorFunctionBase<Real, GradDerived>& grad)
      {
        using Rodin::Index;
        const WNGIRParameters& p = m_parameters;
        Displacement& u = *m_u;
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        using Mesh = std::remove_cvref_t<decltype(mesh)>;
        const std::size_t meshDim = mesh.getDimension();

        WNGIRReport rep;
        const Real h = p.h;
        if (!std::isfinite(h) || !(h > Real(0)))
          Alert::Exception() << "WNGIR requires a finite positive reference mesh size."
                             << Alert::Raise;
        assembleBoundaryConditions();
        buildBoundaryConstraints(mesh);
        for (const Index dof : m_pinnedDofs)
          m_fixedIncrementDOFs[dof] = Real(0);
        m_hingeProblem.setSystemTransform(
          [this](LinearSystemType& system) { applySlipConstraint(system); });
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
          validationCells.push_back(cellIndex);
        }

        if (!p.interfaceAttribute)
        {
          rep.reason = WNGIRReport::Reason::MissingInterface;
          return finish(rep);
        }

        // One bounding-volume locator per solve: the background mesh is fixed
        // for the frame, and the index builds lazily on first query.
        const Location::AABB<Mesh> locator(mesh);
        std::vector<Index> interfaceFacets;
        for (auto face = mesh.getFace(); face; ++face)
          if (face->getAttribute() == *p.interfaceAttribute)
            interfaceFacets.push_back(face->getIndex());
        if (interfaceFacets.empty())
        {
          rep.reason = WNGIRReport::Reason::EmptyInterface;
          return finish(rep);
        }
        const auto normalJump = getNormalJump(mesh, fes, interfaceFacets, meshDim);
        rep.normalJumpRMS = normalJump.rms;
        rep.normalJumpMax = normalJump.max;
        const Real geometricTarget = p.geometricSupTolerance > Real(0)
          ? p.geometricSupTolerance
          : std::pow(h,
              fes.getFiniteElement(meshDim - 1, interfaceFacets.front()).getOrder() + 1);
        if (!std::isfinite(geometricTarget) || !(geometricTarget > Real(0)))
          Alert::Exception() << "WNGIR requires a finite positive geometric target."
                             << Alert::Raise;
        rep.geometricSupTarget = geometricTarget;
        const auto robustScale = getRobustScale(mesh, fes, phi, grad, interfaceFacets, h);
        if (!(robustScale.gradientScale > Real(0)) ||
          !std::isfinite(robustScale.gradientScale))
        {
          rep.reason = WNGIRReport::Reason::DegenerateGradient;
          return finish(rep);
        }
        const Real sigma = robustScale.sigma;
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
        if (initialAdm.inadmissibleCount > 0 ||
          initialAdm.minJ <= acceptedJacobianFloor || initialAdm.maxQ >= p.qMax)
        {
          rep.reason = WNGIRReport::Reason::InvalidInitialGeometry;
          return finish(rep);
        }

        SurfaceState currentSurface = surfaceState(u);
        recordSurfaceState(currentSurface);
        if (!(currentSurface.activeLen > Real(0)))
        {
          rep.reason = WNGIRReport::Reason::EmptyActiveSet;
          return finish(rep);
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
                    << "  inner_total=" << rep.innerIterations
                    << "  inner_last=" << rep.lastInnerIterations
                    << "  inner_converged=" << rep.innerConverged
                    << "  inner_residual=" << rep.innerResidual
                    << "  inner_relative_residual=" << rep.innerRelativeResidual
                    << "  inner_residual_tolerance=" << rep.innerResidualTolerance
                    << "  geom_sup_target=" << geometricTarget
                    << "  inactive_hinge_skips=" << rep.inactiveHingeSkips
                    << "  direct_analyses=" << rep.directAnalyses
                    << "  direct_factorizations=" << rep.directFactorizations
                    << "  fit_energy=" << rep.energy
                    << "  interface_measure=" << rep.interfaceMeasure
                    << "  inner_assembly_seconds=" << rep.tInnerAssembly
                    << "  inner_solve_seconds=" << rep.tInnerSolve
                    << "  inner_ls_seconds=" << rep.tInnerLineSearch
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
            rep.reason = WNGIRReport::Reason::InvalidGeometry;
            return finish(rep);
          }
          if (geometry.sup <= geometricTarget)
          {
            rep.geometricRMS = geometry.rms;
            rep.geometricSup = geometry.sup;
            rep.normalRMS = geometry.normalRMS;
            rep.geometricTargetReached = true;
            rep.qualityBudgetSatisfied = true;
            rep.reason = WNGIRReport::Reason::GeometricTarget;
            return finish(rep);
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
          WNGIRFittingCoefficient obsCoeff(
            grad, u, locator, p, dataNormalization, meshDim);
          auto obsMetric =
            Variational::FaceIntegral(Variational::Dot(obsCoeff * m_duStep, m_vStep));
          obsMetric.setOrder(surfaceOrder);
          obsMetric.over(*p.interfaceAttribute);
          WNGIRFittingForce forceCoeff(
            phi, grad, u, locator, loss, dataNormalization, meshDim);
          auto surfaceForce = Variational::FaceIntegral(forceCoeff, m_vStep);
          surfaceForce.setOrder(surfaceOrder);
          surfaceForce.over(*p.interfaceAttribute);
          // The observation metric and the fitting force depend on the outer
          // displacement, not on the barrier increment, so they are assembled
          // here and reused by every correction below.
          m_fittingMetric = obsMetric;
          m_fittingMetric.assemble();
          if (!m_additionalMetric.getLocalIntegrators().empty() ||
            !m_additionalMetric.getGlobalIntegrators().empty())
          {
            m_additionalMetric.assemble();
            m_fittingMetric.getOperator() += m_additionalMetric.getOperator();
          }
          const Real coefficient = p.h * p.kappaD;
          size_t order = 2;
          for (auto cell = mesh.getCell(); cell; ++cell)
            order = std::max(
              order, 2 * fes.getFiniteElement(meshDim, cell->getIndex()).getOrder());
          if (p.quadratureOrder > 0)
            order = p.quadratureOrder;
          m_distributionForm =
            WNGIRDistribution(m_duStep, m_vStep, u, coefficient, order);
          m_distributionForm.assemble();
          m_fittingForce = surfaceForce;
          m_fittingForce.assemble();
          bool solveOk = true;
          Real predictorAction = Real(0);
          typename ProblemType::ProblemBodyType predictorBody(m_distributionForm);
          predictorBody = predictorBody + m_fittingMetric - m_fittingForce;
          m_metricProblem = predictorBody;
          assembleStep();
          rep.tAssembly += secondsSince(tic);

          Math::SparseMatrix<Real> fixedMetric;
          Math::Vector<Real> fixedForce;
          if constexpr (std::is_same_v<
                          typename FormLanguage::Traits<LinearSystemType>::OperatorType,
                          Math::SparseMatrix<Real>>)
          {
            {
              const auto& system = m_metricProblem.getLinearSystem();
              fixedMetric = system.getOperator();
              fixedForce = system.getVector();
              m_similarityModes = getSimilarityModes(u);
            }
          }
          tic = Clock::now();
          std::size_t predictorIterations = 0;
          Real predictorError = std::numeric_limits<Real>::infinity();
          solveOk =
            solveStep(m_metricProblem, vK, predictorIterations, predictorError, &rep);
          recordLinearSolve(predictorIterations, predictorError);
          rep.tSolve += secondsSince(tic);
          if (!solveOk)
          {
            rep.reason = WNGIRReport::Reason::PredictorFailure;
            break;
          }

          predictor = vK;
          predictorAction = getSurfaceForceAction(mesh, fes, u, predictor, phi, grad,
            interfaceFacets, loss, dataNormalization, meshDim, locator);
          if (!(predictorAction > Real(0)) || !std::isfinite(predictorAction))
          {
            rep.reason = WNGIRReport::Reason::NonDescentPredictor;
            break;
          }
          rep.predictorScale = Real(1);
          if (p.directionalNewton)
          {
            const auto curvatures =
              getSurfaceDirectionalCurvature(mesh, fes, u, predictor, phi, grad,
                interfaceFacets, loss, dataNormalization, meshDim, locator);
            const Real norm = wngirPhysicalDisplacementNorm(
              mesh, fes, validationCells, predictor, p.geometricValidationOrder);
            rep.predictorScale =
              wngirDirectionalNewtonStep(predictorAction, curvatures.first,
                curvatures.second, norm, h * p.directionalNewtonMaxStepOverH);
            if (!(rep.predictorScale > Real(0)) || !std::isfinite(rep.predictorScale))
            {
              rep.reason = WNGIRReport::Reason::InvalidScaling;
              break;
            }
            // The predictor, frozen quadratic form and hinge scale describe the same motion.
            predictor *= rep.predictorScale;
            vK = predictor;
            predictorAction *= rep.predictorScale;
            fixedMetric *= Real(1) / rep.predictorScale;
            m_distributionForm.getOperator() *= Real(1) / rep.predictorScale;
            m_fittingMetric.getOperator() *= Real(1) / rep.predictorScale;
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
            const Real hingeCoefficient =
              domainMeasure > Real(0) ? p.muHat * modelDecrease / domainMeasure : Real(0);
            rep.hingeCoefficient = hingeCoefficient;
            const size_t innerIterations = p.innerIterations;
            const Real predictorNorm =
              std::max(std::abs(predictor.max()), std::abs(predictor.min()));
            const Real residualScale = fixedForce.norm();
            const Real residualTolerance =
              p.innerAbsoluteTolerance + p.innerRelativeTolerance * residualScale;
            bool innerConverged = false;
            rep.lastInnerIterations = 0;
            rep.lastInnerAlpha = 0;
            rep.minInnerAlpha = 1;
            rep.fullInnerSteps = 0;
            if (!hasActiveHinges(
                  mesh, fes, validationCells, u, vK, meshDim, hingeCoefficient))
            {
              const Math::Vector<Real> image = fixedMetric * vK.getData();
              rep.innerResidual = (image - fixedForce).norm();
              rep.innerResidualTolerance = residualTolerance;
              rep.innerRelativeResidual = residualScale > Real(0)
                ? rep.innerResidual / residualScale
                : rep.innerResidual;
              innerConverged = rep.innerResidual <= residualTolerance;
              rep.innerConverged = innerConverged;
              rep.innerRelativeCorrection = Real(0);
              if (p.trace)
                std::cout << "        barrier skip: outer=" << rep.iterations
                          << "  reason=inactive-hinges  residual=" << rep.innerResidual
                          << "  rel=" << rep.innerRelativeResidual
                          << "  converged=" << innerConverged << '\n';
              if (innerConverged)
                ++rep.inactiveHingeSkips;
            }
            if (!innerConverged)
            {
              WNGIRHingeMetric hingeMetric(m_duStep, m_vStep, u, vK, p, hingeCoefficient);
              WNGIRHingeForce hingeForce(m_vStep, u, vK, p, hingeCoefficient);
              typename ProblemType::ProblemBodyType body(m_distributionForm);
              body = body + m_fittingMetric + hingeMetric - m_fittingForce - hingeForce;
              m_hingeProblem = body;
              m_hingeProblem.setState(vK).setBoundaryDOFs(m_fixedIncrementDOFs);
              using Newton = Solver::NewtonSolver<WNGIRLinearSolver<LinearSystemType>>;
              Newton newton(m_linearSolver);
              newton.setMaxIterations(innerIterations)
                .setAbsoluteTolerance(residualTolerance)
                .setRelativeTolerance(Real(0))
                .setStepTolerance(Real(0));
              const auto recordResidual = [&](Real residual, size_t inner) {
                const Real assemblySeconds = m_hingeProblem.getAssemblySeconds();
                rep.tAssembly += assemblySeconds;
                rep.tInnerAssembly += assemblySeconds;
                rep.innerResidual = residual;
                rep.innerResidualTolerance = residualTolerance;
                rep.innerRelativeResidual =
                  residualScale > Real(0) ? residual / residualScale : residual;
                innerConverged = std::isfinite(residual) && residual <= residualTolerance;
                if (p.trace)
                {
                  const auto precision = std::cout.precision();
                  std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
                            << "        barrier residual: outer=" << rep.iterations
                            << "  inner=" << inner << "  residual=" << residual
                            << "  force_norm=" << residualScale
                            << "  rel=" << rep.innerRelativeResidual
                            << "  tolerance=" << residualTolerance
                            << "  converged=" << innerConverged << '\n';
                  std::cout.precision(precision);
                }
              };
              newton.setMonitor([&](const typename Newton::Report& report) {
                if (report.converged || !std::isfinite(report.finalResidual))
                  recordResidual(report.finalResidual, report.iterations);
              });
              newton.setStepPolicy([&](auto&, auto& system, auto& newtonReport) {
                recordResidual(newtonReport.finalResidual, newtonReport.iterations);
                const size_t linearIterations = m_linearSolver.getIterations();
                const Real linearError = m_linearSolver.getError();
                recordLinearSolve(linearIterations, linearError);
                const Real innerSolveSeconds = m_linearSolver.getSolveSeconds();
                rep.tSolve += innerSolveSeconds;
                rep.tInnerSolve += innerSolveSeconds;
                solveOk = m_linearSolver.success();
                if (!solveOk)
                {
                  if (p.trace)
                    std::cout << "        barrier inner=" << (newtonReport.iterations + 1)
                              << "  outer=" << rep.iterations << "  linear_ok=0"
                              << "  cg_it=" << linearIterations
                              << "  cg_err=" << linearError << std::endl;
                  rep.reason = WNGIRReport::Reason::LinearFailure;
                  return typename Newton::StepResult{false, false, Real(0)};
                }
                scratch.setData(system.getSolution());
                const Real correctionNorm =
                  std::max(std::abs(scratch.max()), std::abs(scratch.min()));
                const Real iterateNorm = std::max(std::abs(vK.max()), std::abs(vK.min()));
                const Real correctionScale = std::max(iterateNorm, predictorNorm);
                const Real relativeCorrection = correctionScale > Real(0)
                  ? correctionNorm / correctionScale
                  : (correctionNorm == Real(0) ? Real(0)
                                               : std::numeric_limits<Real>::infinity());
                rep.innerRelativeCorrection = relativeCorrection;
                ++rep.innerIterations;
                ++rep.lastInnerIterations;
                rep.maxInnerIterations =
                  std::max(rep.maxInnerIterations, rep.lastInnerIterations);
                Real innerAlpha = Real(1);
                {
                  const auto meritStart = Clock::now();
                  const auto merit = [&](const Displacement& increment) {
                    Math::Vector<Real> image;
                    image = fixedMetric * increment.getData();
                    const Real quadratic = Real(0.5) * increment.getData().dot(image);
                    const Real force = fixedForce.dot(increment.getData());
                    const Real quality = getHingeEnergy(mesh, fes, validationCells, u,
                      increment, meshDim, hingeCoefficient);
                    return std::pair<Real, Real>{quadratic - force + quality,
                      std::abs(quadratic) + std::abs(force) + std::abs(quality)};
                  };
                  const auto before = merit(vK);
                  const Real slope = -system.getVector().dot(scratch.getData());
                  bool meritAccepted = false;
                  std::size_t backtracks = 0;
                  Real after = std::numeric_limits<Real>::infinity();
                  for (; backtracks < InnerMaxBacktracks; ++backtracks)
                  {
                    uTrial = scratch;
                    uTrial *= innerAlpha;
                    uTrial += vK;
                    const auto trial = merit(uTrial);
                    after = trial.first;
                    const Real rounding = MeritRoundoffFactor *
                      std::numeric_limits<Real>::epsilon() *
                      std::max(before.second, trial.second);
                    if (std::isfinite(before.first) && std::isfinite(after) &&
                      std::isfinite(slope) && slope < Real(0) &&
                      after <= before.first +
                          InnerArmijoCoefficient * innerAlpha * std::min(Real(0), slope) +
                          rounding)
                    {
                      meritAccepted = true;
                      break;
                    }
                    innerAlpha *= BacktrackingReduction;
                  }
                  rep.innerBacktracks += backtracks;
                  rep.tInnerLineSearch += secondsSince(meritStart);
                  if (p.trace)
                  {
                    const auto precision = std::cout.precision();
                    std::cout << std::setprecision(
                                   std::numeric_limits<Real>::max_digits10)
                              << "        barrier merit: outer=" << rep.iterations
                              << "  inner=" << (newtonReport.iterations + 1)
                              << "  before=" << before.first << "  after=" << after
                              << "  slope=" << slope << "  alpha=" << innerAlpha
                              << "  backtracks=" << backtracks
                              << "  accepted=" << meritAccepted << std::endl;
                    std::cout.precision(precision);
                  }
                  if (!meritAccepted)
                  {
                    if (p.trace)
                      std::cout
                        << "        barrier inner=" << (newtonReport.iterations + 1)
                        << "  outer=" << rep.iterations << "  linear_ok=1"
                        << "  merit_ok=0  converged=0  cg_it=" << linearIterations
                        << "  cg_err=" << linearError << std::endl;
                    rep.reason = WNGIRReport::Reason::InnerLineSearchFailure;
                    solveOk = false;
                    return typename Newton::StepResult{false, false, Real(0)};
                  }
                }
                rep.lastInnerAlpha = innerAlpha;
                rep.minInnerAlpha = std::min(rep.minInnerAlpha, innerAlpha);
                if (innerAlpha >= Real(1) - FullStepTolerance)
                  ++rep.fullInnerSteps;
                if (p.trace)
                {
                  const auto precision = std::cout.precision();
                  std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
                            << "        barrier inner=" << (newtonReport.iterations + 1)
                            << "  outer=" << rep.iterations << "  corr=" << correctionNorm
                            << "  iterate=" << iterateNorm
                            << "  rel=" << relativeCorrection << "  alpha=" << innerAlpha
                            << "  linear_ok=1  cg_it=" << linearIterations
                            << "  cg_err=" << linearError
                            << "  mu_eff=" << hingeCoefficient << std::endl;
                  std::cout.precision(precision);
                }
                scratch *= innerAlpha;
                vK += scratch;

                return typename Newton::StepResult{true, false, scratch.getData().norm()};
              });
              newton.solve(vK);
              innerConverged = newton.converged();
              // Certify the last accepted iterate without spending another correction.
              if (solveOk &&
                newton.getReport().reason == Newton::ConvergedReason::MaxIterations)
              {
                m_hingeProblem.assemble();
                recordResidual(
                  m_hingeProblem.getLinearSystem().getVector().norm(), innerIterations);
              }
            }
            rep.innerConverged = innerConverged;
            if (solveOk && !innerConverged)
            {
              rep.reason = WNGIRReport::Reason::InnerIterationLimit;
              solveOk = false;
            }
          }
          if (!solveOk)
            break;

          if (!(std::isfinite(vK.max()) && std::isfinite(vK.min())))
          {
            rep.reason = WNGIRReport::Reason::NonfiniteDirection;
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
            rep.reason = WNGIRReport::Reason::NonDescentDirection;
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
              alpha *= BacktrackingReduction;
            }
          }
          rep.tLineSearch += secondsSince(tic);
          if (!accepted)
          {
            u = previousU;
            rep.reason = WNGIRReport::Reason::LineSearchFailure;
            break;
          }

          rep.lastAlpha = alpha;
          {
            // Measure the accepted FE field, not its coefficient vector.
            scratch = u;
            scratch -= previousU;
            rep.acceptedStep = wngirPhysicalDisplacementNorm(
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
            rep.geometricRMS = geometry.rms;
            rep.geometricSup = geometry.sup;
            rep.normalRMS = geometry.normalRMS;
            rep.qualityBudgetSatisfied =
              rep.minJ > acceptedJacobianFloor && rep.maxQRel < p.qMax;
            rep.geometricTargetReached =
              rep.qualityBudgetSatisfied && geometry.sup <= geometricTarget;
            if (m_monitor)
            {
              auto snapshot = rep;
              ++snapshot.iterations;
              (*m_monitor)(snapshot);
            }
            if (p.trace)
              traceGeometry(geometry, rep.iterations + 1, "accepted");
            if (!std::isfinite(geometry.sup))
            {
              rep.reason = WNGIRReport::Reason::InvalidGeometry;
              ++rep.iterations;
              break;
            }
            geometricTargetReached = geometry.sup <= geometricTarget;
          }
          if (!(surf.activeLen > Real(0)))
          {
            rep.reason = WNGIRReport::Reason::EmptyActiveSet;
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
                      << "  muEff=" << rep.hingeCoefficient
                      << "  pbIt=" << rep.lastInnerIterations
                      << "  pbRel=" << rep.innerRelativeCorrection
                      << "  predictor_action=" << rep.predictorAction
                      << "  inner_residual=" << rep.innerResidual
                      << "  inner_relative_residual=" << rep.innerRelativeResidual
                      << "  desc=" << rep.descentRatio
                      << "  dirNorm=" << rep.directionNormRatio
                      << "  ared/pred=" << rep.actualPredictedDecrease
                      << "  bt=" << backtracks << "  rejJ=" << rep.jacobianRejections
                      << "  rejQ=" << rep.distortionRejections
                      << "  rejE=" << rep.energyRejections << "  min_j=" << rep.minJ
                      << "  max_j=" << rep.maxJ << "  max_Q=" << rep.maxQRel << '\n';

          if (geometricTargetReached)
          {
            rep.reason = WNGIRReport::Reason::GeometricTarget;
            ++rep.iterations;
            break;
          }
          const Real stepThreshold = stepTol + h * p.acceptedStepOverHTol;
          consecutiveSmallAcceptedSteps =
            rep.acceptedStep <= stepThreshold ? consecutiveSmallAcceptedSteps + 1 : 0;
          const Real energyScale =
            std::max({std::abs(ePrev), std::abs(eNow), std::numeric_limits<Real>::min()});
          const Real energyChange = std::abs(ePrev - eNow) / energyScale;
          consecutiveSmallEnergyChanges =
            p.energyStagTol > Real(0) && energyChange <= p.energyStagTol
            ? consecutiveSmallEnergyChanges + 1
            : 0;
          if (consecutiveSmallAcceptedSteps >= p.stagnationIterations ||
            consecutiveSmallEnergyChanges >= p.stagnationIterations)
          {
            rep.reason = consecutiveSmallAcceptedSteps >= p.stagnationIterations
              ? WNGIRReport::Reason::SmallAcceptedSteps
              : WNGIRReport::Reason::SmallEnergyChanges;
            ++rep.iterations;
            break;
          }
          ePrev = eNow;
        }

        const InterfaceGeometryState geometry = getInterfaceGeometryState(
          mesh, fes, u, phi, grad, interfaceFacets, meshDim, locator);
        rep.geometricRMS = geometry.rms;
        rep.geometricSup = geometry.sup;
        rep.normalRMS = geometry.normalRMS;
        rep.qualityBudgetSatisfied =
          rep.minJ > acceptedJacobianFloor && rep.maxQRel < p.qMax;
        rep.geometricTargetReached =
          rep.qualityBudgetSatisfied && geometry.sup <= geometricTarget;
        if (p.trace)
          traceGeometry(geometry, rep.iterations, "final");

        return finish(rep);
      }

    private:
      WNGIRReport finish(const WNGIRReport& report)
      {
        m_report = report;
        if (m_monitor)
          (*m_monitor)(m_report);
        return m_report;
      }

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
              const WNGIRResidualState state(phi, grad, deformation, ip, loss, true);
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
          WNGIRFittingForce force(
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

      /// @brief Integrated affine quality energy on the same quadrature as its force and Hessian.
      template <class Mesh, class FES>
      bool hasActiveHinges(const Mesh& mesh, const FES& fes,
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
              const WNGIRHingeState state(
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
      Real getHingeEnergy(const Mesh& mesh, const FES& fes,
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
              const WNGIRHingeState state(
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
              if (!std::isfinite(j) || j <= m_parameters.jMinRatio ||
                !deformation.isAdmissible())
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
              const auto& moved =
                deformation.getMovedPoint(Variational::IntegrationPoint(point));
              const Real norm = grad.getValue(moved).norm();
              const Real distance = std::abs(phi.getValue(moved)) / norm;
              if (!(norm > gradientFloor) || !std::isfinite(norm) ||
                !std::isfinite(distance))
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
                    fittedNormalNorm > gradientFloor) ||
                !std::isfinite(mappedMeasure) || !std::isfinite(gradientNorm) ||
                !std::isfinite(fittedNormalNorm))
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
        if (!(total.measure > Real(0)) || total.invalidCount > 0)
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
        UnorderedMap<std::uint64_t, std::vector<std::size_t>> incident;
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
          sigma = RobustScaleMeshFloor * h * gradientScale;
          const std::size_t k90 = static_cast<std::size_t>(
            RobustScaleQuantile * static_cast<Real>(residuals.size() - 1));
          std::nth_element(residuals.begin(), residuals.begin() + k90, residuals.end());
          sigma = std::max(sigma, residuals[k90]);
        }
        sigma = std::max(sigma, std::sqrt(std::numeric_limits<Real>::min()));
        return {sigma, gradientScale};
      }

      struct SlipNode
      {
          std::vector<Index> dofs;
          Math::Matrix<Real> complement;
      };

      template <class MeshType>
      void buildBoundaryConstraints(const MeshType& mesh)
      {
        m_slipNodes.clear();
        m_pinnedDofs.clear();
        const auto& fes = m_duStep.getFiniteElementSpace();
        const auto size = static_cast<Eigen::Index>(fes.getSize());
        const std::size_t meshDim = mesh.getDimension();
        const auto& parameters = m_parameters;

        std::map<Index, std::vector<Math::SpatialVector<Real>>> normals;
        for (auto face = mesh.getPolytope(meshDim - 1); face; ++face)
        {
          const auto attribute = face->getAttribute();
          if (!attribute)
            continue;
          if (parameters.fixedBoundaryAttributes.contains(*attribute))
          {
            for (const Index vertex : face->getVertices())
            {
              const auto& dofs = fes.getDOFs(0, vertex);
              for (Eigen::Index k = 0; k < dofs.size(); ++k)
                m_pinnedDofs.insert(static_cast<Index>(dofs(k)));
            }
          }
          if (!parameters.slipBoundaryAttributes.contains(*attribute))
            continue;
          const auto& vertices = face->getVertices();
          if (vertices.size() < meshDim)
            continue;
          Math::SpatialVector<Real> normal(meshDim);
          normal.setZero();
          const auto x0 = mesh.getVertexCoordinates(vertices[0]);
          if (meshDim == 3)
          {
            const Math::SpatialVector<Real> a =
              mesh.getVertexCoordinates(vertices[1]) - x0;
            const Math::SpatialVector<Real> b =
              mesh.getVertexCoordinates(vertices[2]) - x0;
            normal(0) = a(1) * b(2) - a(2) * b(1);
            normal(1) = a(2) * b(0) - a(0) * b(2);
            normal(2) = a(0) * b(1) - a(1) * b(0);
          }
          else
          {
            const Math::SpatialVector<Real> a =
              mesh.getVertexCoordinates(vertices[1]) - x0;
            normal(0) = -a(1);
            normal(1) = a(0);
          }
          const Real length = normal.norm();
          if (!(length > Real(0)))
            continue;
          normal /= length;
          for (const Index vertex : vertices)
            normals[vertex].push_back(normal);
        }

        if (normals.empty())
        {
          m_slipProjector.resize(0, 0);
          return;
        }

        std::vector<Math::SparseTriplet<Real>> entries;
        entries.reserve(static_cast<std::size_t>(size) + 9 * normals.size());
        std::vector<char> constrained(static_cast<std::size_t>(size), 0);
        for (const auto& [vertex, facetNormals] : normals)
        {
          const auto& dofs = fes.getDOFs(0, vertex);
          const auto components = static_cast<Eigen::Index>(dofs.size());
          Math::Matrix<Real> complement(components, components);
          complement.setZero();
          // Gram--Schmidt over the facet normals gives an orthonormal basis of
          // the forbidden subspace, whose rank is one, two, or three.
          std::vector<Math::SpatialVector<Real>> basis;
          for (const auto& candidate : facetNormals)
          {
            Math::SpatialVector<Real> q = candidate;
            for (const auto& previous : basis)
              q -= previous.dot(q) * previous;
            const Real norm = q.norm();
            if (norm > Real(1e-8))
              basis.push_back(q / norm);
          }
          for (const auto& q : basis)
          {
            for (Eigen::Index i = 0; i < components; ++i)
              for (Eigen::Index j = 0; j < components; ++j)
                complement(i, j) += q(i) * q(j);
          }

          SlipNode node;
          node.complement = complement;
          node.dofs.reserve(static_cast<std::size_t>(components));
          for (Eigen::Index k = 0; k < components; ++k)
          {
            const auto dof = static_cast<Index>(dofs(k));
            node.dofs.push_back(dof);
            constrained[static_cast<std::size_t>(dof)] = 1;
          }
          for (Eigen::Index i = 0; i < components; ++i)
          {
            for (Eigen::Index j = 0; j < components; ++j)
            {
              const Real value = (i == j ? Real(1) : Real(0)) - complement(i, j);
              if (value != Real(0))
                entries.emplace_back(
                  static_cast<Math::SparseIndex>(node.dofs[static_cast<std::size_t>(i)]),
                  static_cast<Math::SparseIndex>(node.dofs[static_cast<std::size_t>(j)]),
                  value);
            }
          }
          m_slipNodes.push_back(std::move(node));
        }
        for (Eigen::Index dof = 0; dof < size; ++dof)
        {
          if (!constrained[static_cast<std::size_t>(dof)])
            entries.emplace_back(static_cast<Math::SparseIndex>(dof),
              static_cast<Math::SparseIndex>(dof), Real(1));
        }
        m_slipProjector.resize(size, size);
        m_slipProjector.setFromTriplets(entries.begin(), entries.end());
        m_slipProjector.makeCompressed();
      }

      /**
       * @brief Restricts an assembled step to the slip-admissible subspace.
       *
       * With @f$ P @f$ the projector onto that subspace and @f$ Q=I-P @f$, the
       * step solves @f$ (PAP + \alpha Q)u = Pb @f$. Its solution satisfies
       * @f$ Qu=0 @f$ exactly and @f$ P(Au-b)=0 @f$, so the facets keep their
       * surface to roundoff rather than to a penalty. The diagonal of the
       * eliminated block carries the scale of the row it replaces.
       */
      void applySlipConstraint(LinearSystemType& axb) const
      {
        if (m_slipNodes.empty())
          return;
        if constexpr (std::is_same_v<
                        typename FormLanguage::Traits<LinearSystemType>::OperatorType,
                        Math::SparseMatrix<Real>>)
        {
          auto& A = axb.getOperator();
          auto& b = axb.getVector();
          if (A.rows() != m_slipProjector.rows())
            return;
          A.makeCompressed();

          std::vector<Math::SparseTriplet<Real>> eliminated;
          eliminated.reserve(9 * m_slipNodes.size());
          for (const auto& node : m_slipNodes)
          {
            Real scale = 0;
            for (const Index dof : node.dofs)
              scale += std::abs(
                A.coeff(static_cast<Eigen::Index>(dof), static_cast<Eigen::Index>(dof)));
            scale /= static_cast<Real>(node.dofs.size());
            if (!(scale > Real(0)))
              scale = Real(1);
            const auto components = static_cast<Eigen::Index>(node.dofs.size());
            for (Eigen::Index i = 0; i < components; ++i)
            {
              for (Eigen::Index j = 0; j < components; ++j)
              {
                const Real value = scale * node.complement(i, j);
                if (value != Real(0))
                  eliminated.emplace_back(static_cast<Math::SparseIndex>(
                                            node.dofs[static_cast<std::size_t>(i)]),
                    static_cast<Math::SparseIndex>(
                      node.dofs[static_cast<std::size_t>(j)]),
                    value);
              }
            }
          }
          Math::SparseMatrix<Real> complement(A.rows(), A.cols());
          complement.setFromTriplets(eliminated.begin(), eliminated.end());

          Math::SparseMatrix<Real> restricted =
            (m_slipProjector.transpose() * A * m_slipProjector).eval();
          restricted = (restricted + complement).eval();
          restricted.makeCompressed();
          A = std::move(restricted);
          b = m_slipProjector.transpose() * b;
        }
      }

      void assembleBoundaryConditions()
      {
        m_fixedIncrementDOFs.clear();
        for (auto& condition : m_boundaryConditions)
        {
          condition.assemble();
          const auto* values =
            std::get_if<typename Variational::DirichletBCBase<Real>::ValueDOFs>(
              &condition.getDOFs());
          if (!values)
            Alert::Exception() << "WNGIR supports homogeneous value conditions only."
                               << Alert::Raise;
          for (const auto& [dof, value] : *values)
          {
            if (value != Real(0))
              Alert::Exception()
                << "WNGIR boundary increments must be zero." << Alert::Raise;
            m_fixedIncrementDOFs[dof] = Real(0);
          }
        }
      }

      /// @brief Assemble one symmetric operator for the model, factors and residuals.
      void assembleStep()
      {
        m_metricProblem.assemble();
        using OperatorType =
          typename FormLanguage::Traits<LinearSystemType>::OperatorType;
        if constexpr (std::is_same_v<OperatorType, Math::SparseMatrix<Real>>)
        {
          auto& matrix = m_metricProblem.getLinearSystem().getOperator();
          matrix =
            (Real(0.5) * (matrix + Math::SparseMatrix<Real>(matrix.transpose()))).eval();
          if (!m_fixedIncrementDOFs.empty())
            m_metricProblem.getLinearSystem().eliminate(m_fixedIncrementDOFs);
          applySlipConstraint(m_metricProblem.getLinearSystem());
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
       * @brief Solves the predictor and transfers it to a displacement field.
       * The retained linear adapter also serves the native inner Newton solve.
       */
      bool solveStep(ProblemType& problem, Displacement& out, size_t& iterations,
        Real& error, WNGIRReport* report = nullptr)
      {
        m_linearSolver.setParameters(m_parameters)
          .setSimilarityModes(m_similarityModes)
          .setReport(report);
        auto& system = problem.getLinearSystem();
        m_linearSolver.solve(system);
        iterations = m_linearSolver.getIterations();
        error = m_linearSolver.getError();
        if (!m_linearSolver.success())
          return false;
        m_duStep.getSolution().setData(system.getSolution());
        out = m_duStep.getSolution();
        return true;
      }

      /// Dimensionless sufficient-decrease coefficient for the frozen inner merit.
      static constexpr Real InnerArmijoCoefficient = Real(1e-4);
      /// Safeguard on trials per inner correction, independent of the Newton cap.
      static constexpr size_t InnerMaxBacktracks = 32;
      /// Dimensionless reduction shared by the inner and outer Armijo searches.
      static constexpr Real BacktrackingReduction = Real(0.5);
      /// Relative floating-point allowance in the quadratic inner merit.
      static constexpr Real MeritRoundoffFactor = Real(64);
      /// Dimensionless tolerance for reporting an undamped inner correction.
      static constexpr Real FullStepTolerance = Real(1e-12);
      /// Minimum robust scale in units of the reference mesh spacing.
      static constexpr Real RobustScaleMeshFloor = Real(3);
      /// Residual quantile used by automatic robust-scale selection.
      static constexpr Real RobustScaleQuantile = Real(0.9);

      Displacement* m_u;
      Identifiable::UUID m_trialUUID;
      TrialFunctionType m_duStep;
      TestFunctionType m_vStep;
      ProblemType m_metricProblem;
      WNGIRHingeProblem<TrialFunctionType, TestFunctionType> m_hingeProblem;
      WNGIRLinearSolver<LinearSystemType> m_linearSolver;
      BilinearFormType m_distributionForm;
      Math::Matrix<Real> m_similarityModes;
      /// @brief Observation metric and fitting force at the outer displacement.
      ///
      /// Both depend on the outer displacement only, so they are assembled once
      /// per nonlinear iteration and reused by every barrier correction.
      BilinearFormType m_fittingMetric;
      LinearFormType m_fittingForce;
      BilinearFormType m_additionalMetric;
      Variational::EssentialBoundary<Real> m_boundaryConditions;
      IndexMap<Real> m_fixedIncrementDOFs;
      std::vector<SlipNode> m_slipNodes;
      FlatSet<Index> m_pinnedDofs;
      Math::SparseMatrix<Real> m_slipProjector;
      WNGIRParameters m_parameters;
      WNGIRReport m_report;
      Optional<Monitor> m_monitor;
  };
}

#endif
