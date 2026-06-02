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
#include <unordered_map>
#include <utility>
#include <vector>

#include "Rodin/Assembly.h"
#include "Rodin/Solid/Linear/LinearElasticityIntegral.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Variational.h"

#include "CellGeomCache.h"
#include "JacobianAdmissibilityBarrier.h"
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

        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        auto cellGeom = precomputeCellGeometry(mesh);
        const auto& cellCache = cellGeom.first;
        const auto& cellToLocal = cellGeom.second;

        LSRParameters params = m_params;
        completeParameters(params, sLF);

        TrialFunction du(fes);
        TestFunction  v(fes);
        auto zero = VectorFunction{ Zero(), Zero() };

        RealFunction<Real> shapeWeight(params.shapeWeight);

        BarrierParameters barrierParams;
        barrierParams.jMin = params.jMinRatio;
        barrierParams.domainMeasure = computeDomainMeasure();
        barrierParams.qBarrierWeight = params.qBarrierWeight;
        barrierParams.qBarrierAct    = params.qBarrierAct;
        barrierParams.qBarrierMax    = params.qBarrierMax;

        LSRRegistration lsrTerm(du, v, u);
        JacobianAdmissibilityBarrier barrier(
            du, v, u, cellCache, cellToLocal);

        if (params.initialGuess == LSRInitialGuess::Zero)
          u.getData().setZero();
        else if (params.initialGuess == LSRInitialGuess::Hilbert)
          applyHilbertInitialGuess(
              sLF, phi, gradPhi, shapeWeight,
              barrierParams, cellCache, cellToLocal, params);

        Problem newton(du, v);
        switch (params.tangent)
        {
          case LSRTangent::GaussNewton:
            newton =
                lsrTerm.Tangent(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
              + lsrTerm.Residual(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
              + barrier.Tangent(shapeWeight, barrierParams)
              + barrier.Residual(shapeWeight, barrierParams)
              + DirichletBC(du, zero);
            break;
          case LSRTangent::Newton:
            newton =
                lsrTerm.Tangent(phi, gradPhi, hessPhi, sLF, makeLSRIntegratorParameters(params))
              + lsrTerm.Residual(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
              + barrier.Tangent(shapeWeight, barrierParams)
              + barrier.Residual(shapeWeight, barrierParams)
              + DirichletBC(du, zero);
            break;
          case LSRTangent::PSDProjectedNewton:
            newton =
                lsrTerm.TangentPSDProjected(
                    phi, gradPhi, hessPhi, sLF, makeLSRIntegratorParameters(params))
              + lsrTerm.Residual(phi, gradPhi, sLF, makeLSRIntegratorParameters(params))
              + barrier.Tangent(shapeWeight, barrierParams)
              + barrier.Residual(shapeWeight, barrierParams)
              + DirichletBC(du, zero);
            break;
        }

        Solver::SparseLU linearSolver(newton);
        Solver::NewtonSolver newtonSolver(linearSolver);
        newtonSolver
          .setMaxIterations(1)
          .setDampingFactor(Real(1))
          .setAbsoluteTolerance(Real(0))
          .setRelativeTolerance(Real(0))
          .setStepTolerance(params.stepTolerance);

        auto evaluator = [&](const Math::Vector<Real>& uTry)
        {
          return evaluateLSRAdmissibilityP1(
              uTry, cellCache, fes, params.jMinRatio);
        };

        const Real jLineSearchRatio =
          std::max(params.jMinRatio,
                   params.lineSearchSafetyMargin * params.jSafeRatio);
        m_report.jLineSearchRatio = jLineSearchRatio;

        Math::Vector<Real> uBest = u.getData();
        Real residualBest = std::numeric_limits<Real>::infinity();
        std::size_t stallCount = 0;
        Real r0 = 0;

        for (std::size_t it = 0; it < params.maxNewtonIterations; ++it)
        {
          const Math::Vector<Real> uOld = u.getData();

          newtonSolver.solve(u.getData());
          const auto step = newtonSolver.getReport();
          const Real residual = step.final_residual;
          if (it == 0)
          {
            r0 = residual;
            m_report.initialResidual = residual;
          }

          const bool residualImproved = (residual < residualBest);
          if (residualImproved)
          {
            residualBest = residual;
            stallCount = 0;
          }
          else
          {
            ++stallCount;
          }

          Math::Vector<Real> p = u.getData() - uOld;
          u.getData() = uOld;

          const auto ls = runLSRAdmissibilityLineSearch(
              u.getData(), p, evaluator, jLineSearchRatio,
              params.alphaInit, params.alphaReduction, params.alphaMin,
              params.qShapeMax);

          m_report.iterations = it + 1;
          m_report.finalResidual = residual;
          m_report.finalStepNorm = ls.alphaAccepted * p.norm();
          m_report.lastAcceptedAlpha = ls.alphaAccepted;
          m_report.totalBacktracks += ls.backtracks;
          m_report.minJRatio = ls.minJRatioAccepted;
          m_report.inadmissibleCount = ls.inadmissibleCountAccepted;

          if (!ls.succeeded)
          {
            m_report.lineSearchFailed = true;
            m_report.converged = false;
            u.getData() = uBest;
            return m_report;
          }

          if (residualImproved)
            uBest = u.getData();

          if (residual <= params.absoluteTolerance
              || (r0 > 0 && residual <= params.relativeTolerance * r0))
          {
            m_report.converged = true;
            return m_report;
          }

          if (m_report.finalStepNorm <= params.stepTolerance)
          {
            m_report.converged = true;
            return m_report;
          }

          if (params.stallPatience > 0 && stallCount >= params.stallPatience)
          {
            m_report.converged = true;
            u.getData() = uBest;
            return m_report;
          }
        }

        u.getData() = uBest;
        m_report.converged = false;
        m_report.finalResidual = residualBest;
        return m_report;
      }

    private:
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
          const std::vector<CellGeomCache>& cellCache,
          const std::unordered_map<Index, std::size_t>& cellToLocal,
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
            JacobianAdmissibilityBarrier barrier(
                duH, vH, u, cellCache, cellToLocal);
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
  };

  template <class Displacement>
  LSR(Displacement&) -> LSR<Displacement>;
}

#endif
