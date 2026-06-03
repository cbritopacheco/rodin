/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_LSR_H
#define RODIN_ADAPTATION_LSR_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Rodin/Assembly.h"
#include "Rodin/Solid/Linear/LinearElasticityIntegral.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Variational.h"

#include "CellGeomCache.h"
#include "JacobianAdmissibilityBarrierSampled.h"
#include "LSRAdmissibility.h"
#include "LSRParameters.h"
#include "LSRReport.h"
#include "LSRRegistration.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Level-set registration solver.
   *
   * `LSR` owns no mesh data. It is bound to the displacement field `u`,
   * infers the carrier finite-element space and mesh from that field, and
   * overwrites `u` with the accepted registration displacement.
   */
  template <class Displacement>
  class LSR
  {
    public:
      explicit LSR(Displacement& u)
        : m_u(u)
      {}

      LSR& setParameters(LSRParameters params)
      {
        m_params = std::move(params);
        return *this;
      }

      const LSRParameters& getParameters() const noexcept
      {
        return m_params;
      }

      const LSRReport& getReport() const noexcept
      {
        return m_report;
      }

      template <class SLF, class PhiDerived, class GradDerived, class HessDerived>
      LSRReport solve(
          const SLF& sLF,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& gradPhi,
          const Variational::MatrixFunctionBase<Real, HessDerived>& hessPhi)
      {
        using Variational::DirichletBC;
        using Variational::Integral;
        using Variational::Jacobian;
        using Variational::Problem;
        using Variational::RealFunction;
        using Variational::TestFunction;
        using Variational::TrialFunction;
        using Variational::VectorFunction;
        using Variational::Zero;

        m_report = {};
        m_prevAcceptedAlpha = Real(0);

        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        LSRParameters params = m_params;
        completeParameters(params, sLF);
        const bool barrierEnabled =
             params.shapeWeight != Real(0)
          || params.qBarrierWeight != Real(0);

        TrialFunction du(fes);
        TestFunction  v(fes);
        auto zero = VectorFunction{ Zero(), Zero() };

	        RealFunction<Real> shapeWeight(params.shapeWeight);
	        RealFunction<Real> h1RegularizationWeight(
	            params.h1RegularizationWeight);

        BarrierParameters barrierParams;
        barrierParams.jMin = params.jMinRatio;
        barrierParams.domainMeasure = computeDomainMeasure();
        barrierParams.qBarrierWeight = params.qBarrierWeight;
        barrierParams.qBarrierAct    = params.qBarrierAct;
        barrierParams.qBarrierMax    = params.qBarrierMax;

        LSRRegistration lsrTerm(du, v, u);
        JacobianAdmissibilityBarrierSampled barrier(
            du, v, u, params.quadratureOrder);

        if (params.initialGuess == LSRInitialGuess::Zero)
          u.getData().setZero();
        else if (params.initialGuess == LSRInitialGuess::Hilbert)
          applyHilbertInitialGuess(
              sLF, phi, gradPhi, shapeWeight,
              barrierParams, params);

        Problem newton(du, v);
        switch (params.tangent)
        {
          case LSRTangent::GaussNewton:
            if (barrierEnabled)
	              newton =
	                  lsrTerm.Tangent(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
	                + barrier.Tangent(shapeWeight, barrierParams)
	                + barrier.Residual(shapeWeight, barrierParams)
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
	                + DirichletBC(du, zero);
	            else
	              newton =
	                  lsrTerm.Tangent(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
	                + DirichletBC(du, zero);
            break;
          case LSRTangent::Newton:
            if (barrierEnabled)
	              newton =
	                  lsrTerm.Tangent(phi, gradPhi, hessPhi, sLF, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
	                + barrier.Tangent(shapeWeight, barrierParams)
	                + barrier.Residual(shapeWeight, barrierParams)
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
	                + DirichletBC(du, zero);
	            else
	              newton =
	                  lsrTerm.Tangent(phi, gradPhi, hessPhi, sLF, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
	                + DirichletBC(du, zero);
            break;
          case LSRTangent::PSDProjectedNewton:
            if (barrierEnabled)
	              newton =
	                  lsrTerm.TangentPSDProjected(
	                      phi, gradPhi, hessPhi, sLF, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
	                + barrier.Tangent(shapeWeight, barrierParams)
	                + barrier.Residual(shapeWeight, barrierParams)
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
	                + DirichletBC(du, zero);
	            else
	              newton =
	                  lsrTerm.TangentPSDProjected(
	                      phi, gradPhi, hessPhi, sLF, makeLSRIntegratorParameters(params))
	                + lsrTerm.Residual(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
	                + Integral(h1RegularizationWeight * Jacobian(du), Jacobian(v))
	                + Integral(h1RegularizationWeight * Jacobian(u), Jacobian(v))
	                + DirichletBC(du, zero);
            break;
        }

        Solver::SparseLU linearSolver(newton);
        Solver::NewtonSolver newtonSolver(linearSolver);

          auto evaluator = [&](const Math::Vector<Real>& uTry)
          {
            return evaluateLSRAdmissibilitySampled(
                u, uTry, params.jMinRatio, params.quadratureOrder);
          };

        const Real jLineSearchRatio =
          std::max(params.jMinRatio,
                   params.lineSearchSafetyMargin * params.jSafeRatio);
        m_report.jLineSearchRatio = jLineSearchRatio;

          auto objective = [&](const Math::Vector<Real>& uTry) -> Real
          {
            u.getData() = uTry;
            return computeObjective(
                sLF, phi, params, barrierParams,
                barrierEnabled);
          };

        auto residualNorm = [&](const Math::Vector<Real>& uTry) -> Real
        {
          u.getData() = uTry;
          newton.assemble();
          return newton.getLinearSystem().getVector().norm();
        };

        Math::Vector<Real> uBest = u.getData();
        Real energyBest = std::numeric_limits<Real>::infinity();
        Real residualBest = std::numeric_limits<Real>::infinity();
        std::size_t stallCount = 0;
        bool initialized = false;

        newtonSolver
          .setMaxIterations(params.maxNewtonIterations)
          .setAbsoluteTolerance(params.absoluteTolerance)
          .setRelativeTolerance(params.relativeTolerance)
          .setStepTolerance(Real(0))
          .setDampingFactor(Real(1));

        newtonSolver.setStepPolicy(
          typename decltype(newtonSolver)::StepPolicy(
          [&](auto& x, auto& linearSystem, auto& solverReport)
            -> typename decltype(newtonSolver)::StepResult
        {
          const Math::Vector<Real> uOld = x;
          const Real residual = solverReport.final_residual;
          const Real energy = objective(uOld);

          if (!std::isfinite(residual) || !std::isfinite(energy))
          {
            m_report.converged = false;
            x = uBest;
            return {false, false, Real(0)};
          }

          if (!initialized)
          {
            m_report.initialResidual = residual;
            m_report.initialEnergy = energy;
            uBest = uOld;
            energyBest = energy;
            residualBest = residual;
            initialized = true;
          }

          Math::Vector<Real> newtonDirection = linearSystem.getSolution();
          Math::Vector<Real> residualDirection = linearSystem.getVector();
          const Real newtonDirectionNorm = newtonDirection.norm();
          if (!std::isfinite(newtonDirectionNorm))
          {
            m_report.converged = false;
            x = uBest;
            return {false, false, Real(0)};
          }

          LSRLineSearchResult ls;
          Real acceptedResidual = std::numeric_limits<Real>::quiet_NaN();
          Real acceptedEnergy = std::numeric_limits<Real>::quiet_NaN();
          Math::Vector<Real> acceptedU = uOld;
          Real acceptedStepNorm = Real(0);
          LSRAdmissibilityReport lastAdm;
          const Real energyTol =
            params.energyDecreaseTolerance
            * std::max<Real>(Real(1), std::abs(energy));

          auto tryLineSearchDirection =
            [&](const Math::Vector<Real>& direction) -> bool
          {
            const Real directionNorm = direction.norm();
            if (!std::isfinite(directionNorm) || directionNorm == Real(0))
              return false;

            // Warm-start the initial step from the previously accepted
            // alpha (no-op on the first iteration and whenever the
            // previous iteration accepted alpha = 1).
            Real alphaStart = params.alphaInit;
            if (params.useWarmStartAlpha
                && m_prevAcceptedAlpha > Real(0))
            {
              alphaStart = std::min<Real>(
                  Real(1),
                  params.alphaWarmStartGrowth * m_prevAcceptedAlpha);
              alphaStart = std::max<Real>(alphaStart, params.alphaMin);
            }

            Real alpha = alphaStart;
            while (alpha >= params.alphaMin)
            {
              const Math::Vector<Real> uTrial = uOld + alpha * direction;
              const LSRAdmissibilityReport adm = evaluator(uTrial);

              lastAdm = adm;

              const bool jOK =
                   adm.minJRatio > jLineSearchRatio
                && adm.inadmissibleCount == 0;
              const bool qOK = adm.maxQShape <= params.qShapeMax;

              bool energyOK = false;
              Real trialEnergy = std::numeric_limits<Real>::infinity();
              if (jOK && qOK)
              {
                trialEnergy = objective(uTrial);
                energyOK =
                     std::isfinite(trialEnergy)
                  && trialEnergy <= energy + energyTol;
              }

              if (jOK && qOK && energyOK)
              {
                acceptedResidual = residualNorm(uTrial);
                if (std::isfinite(acceptedResidual))
                {
                  acceptedEnergy = trialEnergy;
                  acceptedU = uTrial;
                  acceptedStepNorm = alpha * directionNorm;
                  ls.alphaAccepted = alpha;
                  ls.minJRatioAccepted = adm.minJRatio;
                  ls.inadmissibleCountAccepted = adm.inadmissibleCount;
                  ls.maxQShapeAccepted = adm.maxQShape;
                  ls.succeeded = true;
                  return true;
                }
              }

              alpha *= params.alphaReduction;
              ++ls.backtracks;
            }
            return false;
          };

          const bool acceptedNewton =
            tryLineSearchDirection(newtonDirection);
          if (!acceptedNewton)
          {
            const bool acceptedResidualDirection =
              tryLineSearchDirection(residualDirection);
            if (!acceptedResidualDirection)
              tryLineSearchDirection(-residualDirection);
          }

          if (!ls.succeeded)
          {
            x = uOld;
            ls.minJRatioAccepted = lastAdm.minJRatio;
            ls.inadmissibleCountAccepted = lastAdm.inadmissibleCount;
            ls.maxQShapeAccepted = lastAdm.maxQShape;
          }
          else
          {
            x = acceptedU;
          }

          m_report.iterations = solverReport.iterations + 1;
          m_report.finalResidual =
            ls.succeeded ? acceptedResidual : residualBest;
          m_report.finalEnergy =
            ls.succeeded ? acceptedEnergy : energyBest;
          m_report.finalStepNorm = acceptedStepNorm;
          m_report.lastAcceptedAlpha = ls.alphaAccepted;
          m_report.totalBacktracks += ls.backtracks;
          if (ls.succeeded)
            m_prevAcceptedAlpha = ls.alphaAccepted;
          m_report.minJRatio = ls.minJRatioAccepted;
          m_report.inadmissibleCount = ls.inadmissibleCountAccepted;

          if (!ls.succeeded)
          {
            m_report.lineSearchFailed = true;
            m_report.converged = false;
            x = uBest;
            solverReport.final_residual = m_report.finalResidual;
            solverReport.final_step_norm = m_report.finalStepNorm;
            return {false, false, m_report.finalStepNorm};
          }

          const bool energyImproved =
            acceptedEnergy < energyBest
              - params.energyDecreaseTolerance
                * std::max<Real>(Real(1), std::abs(energyBest));
          if (energyImproved)
          {
            uBest = u.getData();
            energyBest = acceptedEnergy;
            residualBest = acceptedResidual;
            stallCount = 0;
          }
          else
          {
            ++stallCount;
          }

          if (params.acceptedStateConvergenceTest
              && params.acceptedStateConvergenceTest(m_report))
          {
            m_report.converged = true;
            solverReport.final_residual = m_report.finalResidual;
            solverReport.final_step_norm = m_report.finalStepNorm;
            return {true, true, m_report.finalStepNorm};
          }

          if (acceptedResidual <= params.absoluteTolerance
              || (m_report.initialResidual > 0
                  && acceptedResidual <= params.relativeTolerance
                    * m_report.initialResidual))
          {
            m_report.converged = true;
            solverReport.final_residual = m_report.finalResidual;
            solverReport.final_step_norm = m_report.finalStepNorm;
            return {true, true, m_report.finalStepNorm};
          }

          if (params.stepTolerance > Real(0)
              && m_report.finalStepNorm <= params.stepTolerance)
          {
            m_report.converged = true;
            solverReport.final_residual = m_report.finalResidual;
            solverReport.final_step_norm = m_report.finalStepNorm;
            return {true, true, m_report.finalStepNorm};
          }

          if (params.stallPatience > 0 && stallCount >= params.stallPatience)
          {
            m_report.converged = true;
            x = uBest;
            m_report.finalResidual = residualBest;
            m_report.finalEnergy = energyBest;
            solverReport.final_residual = m_report.finalResidual;
            solverReport.final_step_norm = m_report.finalStepNorm;
            return {true, true, m_report.finalStepNorm};
          }

          solverReport.final_residual = m_report.finalResidual;
          solverReport.final_step_norm = m_report.finalStepNorm;
          return {true, false, m_report.finalStepNorm};
        }));

        newtonSolver.solve(u.getData());
        const auto& solverReport = newtonSolver.getReport();

        if (!initialized)
        {
          m_report.iterations = solverReport.iterations;
          m_report.initialResidual = solverReport.initial_residual;
          m_report.finalResidual = solverReport.final_residual;
          m_report.converged = solverReport.converged;
          return m_report;
        }

        if (!solverReport.converged && !m_report.lineSearchFailed)
        {
          u.getData() = uBest;
          m_report.finalResidual = residualBest;
          m_report.finalEnergy = energyBest;
        }
        m_report.converged = solverReport.converged;
        return m_report;
      }

    private:
      template <class SLF, class PhiDerived>
      Real computeObjective(
          const SLF& sLF,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const LSRParameters& params,
          const BarrierParameters& barrierParams,
          bool barrierEnabled) const
      {
        const auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        Real energy = 0;
        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        {
          const auto& cell = *cellIt;
          const auto& fe = fes.getFiniteElement(
              cell.getDimension(), cell.getIndex());
	          const auto& qf =
	            QF::PolytopeQuadratureFormula::get(
	                lsrQuadOrderFor(fe.getOrder(), makeLSRIntegratorParameters(params)),
	                cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& pt = quadrature.getPoint(q);
            const Variational::IntegrationPoint ip(pt, &qf, q);
            const Real s = sLF.getValue(ip);
            const auto uq = u.getValue(ip);
            Math::SpatialVector<Real> y(cell.getDimension());
            for (std::size_t c = 0; c < fes.getVectorDimension(); ++c)
              y(c) = pt.getCoordinates()(c) + uq(c);

            const Geometry::Point movedPoint =
              detail::makeMovedPoint(pt, y);
            const Real r = phi.getValue(movedPoint) - s;
            const Real W =
              std::exp(-s * s / (2 * params.deltaW * params.deltaW));
            energy += Real(0.5)
              * qf.getWeight(q) * pt.getDistortion()
              * params.rhoS * W * params.normalizer * r * r;
          }
        }

        if (barrierEnabled)
        {
          for (auto cellIt2 = mesh.getCell(); cellIt2; ++cellIt2)
          {
            const auto& barrierCell = *cellIt2;
            const Real e = computeBarrierSampledCellEnergy(
                barrierCell, u, params.shapeWeight, barrierParams);
            if (!std::isfinite(e))
              return std::numeric_limits<Real>::infinity();
            energy += e;
          }
        }
	        if (params.h1RegularizationWeight != Real(0))
	          energy += computeH1RegularizationEnergy(params);
	        return energy;
	      }

	      Real computeH1RegularizationEnergy(const LSRParameters& params) const
	      {
	        using Variational::IntegrationPoint;
	        using Variational::Jacobian;

	        const auto& u = m_u.get();
	        const auto& fes = u.getFiniteElementSpace();
	        const auto& mesh = fes.getMesh();
	        auto gradU = Jacobian(u);

	        Real energy = 0;
	        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
	        {
	          const auto& cell = *cellIt;
	          const auto& fe =
	            fes.getFiniteElement(cell.getDimension(), cell.getIndex());
	          const auto& qf =
	            QF::PolytopeQuadratureFormula::get(
	                lsrQuadOrderFor(
	                    fe.getOrder(), makeLSRIntegratorParameters(params)),
	                cell.getGeometry());
	          const auto& quadrature = cell.getQuadrature(qf);
	          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
	          {
	            const auto& pt = quadrature.getPoint(q);
	            const IntegrationPoint ip(pt, &qf, q);
	            energy += Real(0.5)
	              * params.h1RegularizationWeight
	              * qf.getWeight(q)
	              * pt.getDistortion()
	              * gradU.getValue(ip).squaredNorm();
	          }
	        }
	        return energy;
	      }


      template <class SLF>
      Real computeWeightedBandMeasure(const SLF& sLF, Real deltaW) const
      {
        const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
        Real weightedBandMeasure = 0;

        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        {
          const auto& cell = *cellIt;
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(cell.getDimension(), cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& pt = quadrature.getPoint(q);
            Math::SpatialMatrix<Real> J;
            cell.getTransformation().jacobian(J, pt.getReferenceCoordinates());
            const Real detJ = std::abs(J.determinant());
            const Real wq = qf.getWeight(q) * detJ;
            const Real s = sLF.getValue(pt);
            const Real W = std::exp(-s * s / (2 * deltaW * deltaW));
            weightedBandMeasure += W * wq;
          }
        }
        return weightedBandMeasure;
      }

      Real computeDomainMeasure() const
      {
        const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
        Real measure = 0;
        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        {
          const auto& cell = *cellIt;
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(cell.getDimension(), cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          for (std::size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& pt = quadrature.getPoint(q);
            Math::SpatialMatrix<Real> J;
            cell.getTransformation().jacobian(J, pt.getReferenceCoordinates());
            measure += qf.getWeight(q) * std::abs(J.determinant());
          }
        }
        return measure;
      }

      template <class SLF>
      void completeParameters(LSRParameters& params, const SLF& sLF) const
      {
        const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
        const Real domainMeasure = computeDomainMeasure();

        if (params.hRef <= 0)
          params.hRef = std::sqrt(domainMeasure / static_cast<Real>(mesh.getCellCount()));
        if (params.deltaW <= 0)
          params.deltaW = params.hRef;
        if (params.normalizer <= 0)
        {
          const Real weightedBandMeasure = computeWeightedBandMeasure(sLF, params.deltaW);
          if (weightedBandMeasure <= 0)
            throw std::runtime_error("LSR: weighted level-set band is empty.");
          params.normalizer = Real(1) / (weightedBandMeasure * params.hRef * params.hRef);
        }
      }

      Real initialGuessScaleFactor(const LSRParameters& params) const
      {
        switch (params.initialGuessScaling)
        {
          case LSRInitialGuessScaling::Unnormalized:
            return Real(1);
          case LSRInitialGuessScaling::EnergyNormalized:
            return params.normalizer;
          case LSRInitialGuessScaling::BandNormalized:
            return params.normalizer * params.hRef * params.hRef;
        }
        return Real(1);
      }

      template <class SLF, class PhiDerived, class GradDerived,
                class ShapeWeightDerived>
      void applyHilbertInitialGuess(
          const SLF& sLF,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& gradPhi,
          const Variational::RealFunctionBase<ShapeWeightDerived>& shapeWeight,
          const BarrierParameters& barrierParams,
          const LSRParameters& params)
      {
        using Variational::DirichletBC;
        using Variational::Integral;
        using Variational::Jacobian;
        using Variational::LinearElasticityIntegral;
        using Variational::Problem;
        using Variational::RealFunction;
        using Variational::TestFunction;
        using Variational::TrialFunction;
        using Variational::VectorFunction;
        using Variational::Zero;

        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();

        TrialFunction duH(fes);
        TestFunction vH(fes);
        auto zero = VectorFunction{ Zero(), Zero() };
        const Real scale = initialGuessScaleFactor(params);

        RealFunction rhsScalar(
            [&](const Geometry::Point& p) -> Real
            {
              const Real phiVal = phi.getValue(p);
              const Real sVal = sLF.getValue(p);
              const Real W =
                std::exp(-sVal * sVal / (2 * params.deltaW * params.deltaW));
              return scale * params.rhoS * W * (phiVal - sVal);
            });

        u.getData().setZero();

        Problem hilbert(duH, vH);
        switch (params.initialGuessMetric)
        {
          case LSRHilbertMetric::Harmonic:
            hilbert =
                Integral(Jacobian(duH), Jacobian(vH))
              + Integral(rhsScalar * gradPhi, vH)
              + DirichletBC(duH, zero);
            break;
          case LSRHilbertMetric::Elasticity:
            hilbert =
                LinearElasticityIntegral(duH, vH)(
                    params.initialGuessElasticityLambda,
                    params.initialGuessElasticityMu)
              + Integral(rhsScalar * gradPhi, vH)
              + DirichletBC(duH, zero);
            break;
          case LSRHilbertMetric::ShapeHessian:
          {
            JacobianAdmissibilityBarrierSampled barrier(
                duH, vH, u, params.quadratureOrder);
            hilbert =
                barrier.Tangent(shapeWeight, barrierParams)
              + Integral(rhsScalar * gradPhi, vH)
              + DirichletBC(duH, zero);
            break;
          }
        }
        Solver::SparseLU(hilbert).solve();
        u.getData() = duH.getSolution().getData();
      }

      std::reference_wrapper<Displacement> m_u;
      LSRParameters m_params;
      LSRReport m_report;
      Real m_prevAcceptedAlpha = Real(0);
  };

  template <class Displacement>
  LSR(Displacement&) -> LSR<Displacement>;
}

#endif
