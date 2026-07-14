/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRSOLVER_H
#define RODIN_ADAPTATION_WNGIRSOLVER_H

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "Rodin/Solver/CG.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/Trace.h"
#include "Rodin/Variational/IdentityMatrix.h"

#include "WNGIRAdmissibilityMetric.h"
#include "WNGIRNormalOffsetInitializer.h"
#include "WNGIRParameters.h"
#include "WNGIRReport.h"
#include "WNGIRSurfaceForce.h"
#include "WNGIRSurfaceObservationMetric.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Backend-independent WNGIR mesh-fitting solver.
   *
   * Fits the interface skeleton of a mesh to the zero level set of @f$\phi@f$
   * by repeatedly solving a regularized first-derivative metric problem and
   * line-searching the resulting displacement on the true geometry.
   *
   * Constructed from backend-matched trial and test functions, the solver uses
   * the ordinary Rodin form-language deduction path for the step
   * @ref Rodin::Variational::Problem, the metric @ref Rodin::Variational::BilinearForm,
   * and @ref Rodin::Solver::CG. Every global quantity is a backend-matched
   * GridFunction; only element-local geometry uses @ref Math::SpatialVector /
   * @ref Math::SpatialMatrix. WNGIR contains no backend-specific code.
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

      template <class Mesh, class PhiDerived, class GradDerived>
      /// @brief Solves the WNGIR fitting problem on marked interface facets.
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
        const Real gammaM = p.gammaM >= Real(0) ? p.gammaM : Real(1) / h;
        const Real gammaH = p.gammaH >= Real(0) ? p.gammaH : Real(1) / h;
        const Real gammaDiv = p.gammaDiv >= Real(0) ? p.gammaDiv : gammaH;
        const Real ellM = p.ellM >= Real(0) ? p.ellM : Real(3) * h;
        const Real activeRMSTol =
          p.activeRMSTol > Real(0) ? p.activeRMSTol : Real(4) * h * h;
        const Real activeSupTol =
          p.activeSupTol > Real(0) ? p.activeSupTol : Real(10) * h * h;
        const Real stepTol = p.stepTol > Real(0) ? p.stepTol : Real(1e-4) * h;
        constexpr Real epsG = Real(1e-12);
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
        struct FacetQP
        {
            Real w = 0;
            SpatialVec X;
            Math::SpatialPoint rc;
        };
        struct FacetTable
        {
            Index facet = 0;
            std::vector<FacetQP> qps;
        };
        std::vector<FacetTable> facetTables;
        facetTables.reserve(interfaceFacets.size());
        for (const Index facetIdx : interfaceFacets)
        {
          const auto face = mesh.getFace(facetIdx);
          const auto& fe = fes.getFiniteElement(meshDim - 1, facetIdx);
          const std::size_t feOrder = fe.getOrder();
          const std::size_t qOrder = p.quadratureOrder > 0
            ? p.quadratureOrder
            : std::max<std::size_t>(2, 2 * feOrder);
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(qOrder, face->getGeometry());
          const auto& quad = face->getQuadrature(qf);

          FacetTable ft;
          ft.facet = facetIdx;
          ft.qps.resize(quad.getSize());
          for (std::size_t q = 0; q < quad.getSize(); ++q)
          {
            const auto& pt = quad.getPoint(q);
            auto& fq = ft.qps[q];
            fq.w = qf.getWeight(q) * pt.getDistortion();
            const auto& pc = pt.getPhysicalCoordinates();
            fq.X = SpatialVec::Zero(meshDim);
            for (std::size_t r = 0; r < meshDim; ++r)
              fq.X(static_cast<Eigen::Index>(r)) = pc(r);
            fq.rc = pt.getReferenceCoordinates();
          }
          facetTables.push_back(std::move(ft));
        }

        auto clampUnit = [](Real x) { return std::max(Real(-1), std::min(Real(1), x)); };
        auto edgeKey = [](Index a, Index b) -> std::uint64_t {
          const auto lo = static_cast<std::uint64_t>(std::min(a, b));
          const auto hi = static_cast<std::uint64_t>(std::max(a, b));
          return (lo << 32) ^ hi;
        };
        std::vector<SpatialVec> facetNormals(interfaceFacets.size());
        std::vector<Real> facetMeasures(interfaceFacets.size(), Real(0));
        if (meshDim == 2)
        {
          for (std::size_t i = 0; i < interfaceFacets.size(); ++i)
          {
            const auto face = mesh.getFace(interfaceFacets[i]);
            SpatialVec n = SpatialVec::Zero(meshDim);
            Real measure = 0;
            for (const auto& fq : facetTables[i].qps)
            {
              const Geometry::Point pt(*face, fq.rc, fq.X);
              const auto& J = pt.getJacobian();
              SpatialVec t = SpatialVec::Zero(meshDim);
              for (std::size_t r = 0; r < meshDim; ++r)
                t(static_cast<Eigen::Index>(r)) = J(r, 0);
              SpatialVec nq(2);
              nq(0) = t(1);
              nq(1) = -t(0);
              const Real nn = nq.norm();
              if (nn > Real(1e-30))
              {
                n += fq.w * (nq / nn);
                measure += fq.w;
              }
            }
            if (n.norm() <= Real(1e-30))
            {
              n.setZero();
              n(0) = Real(1);
            }
            else
              n /= n.norm();
            facetMeasures[i] = measure;
            facetNormals[i] = n;
          }
        }
        else if (meshDim == 3)
        {
          for (std::size_t i = 0; i < interfaceFacets.size(); ++i)
          {
            const auto face = mesh.getFace(interfaceFacets[i]);
            SpatialVec n = SpatialVec::Zero(meshDim);
            Real measure = 0;
            for (const auto& fq : facetTables[i].qps)
            {
              const Geometry::Point pt(*face, fq.rc, fq.X);
              const auto& J = pt.getJacobian();
              SpatialVec a = SpatialVec::Zero(meshDim);
              SpatialVec b = SpatialVec::Zero(meshDim);
              for (std::size_t r = 0; r < meshDim; ++r)
              {
                a(static_cast<Eigen::Index>(r)) = J(r, 0);
                b(static_cast<Eigen::Index>(r)) = J(r, 1);
              }
              SpatialVec nq = a.cross(b);
              const Real nn = nq.norm();
              if (nn > Real(1e-30))
              {
                n += fq.w * (nq / nn);
                measure += fq.w;
              }
            }
            if (n.norm() <= Real(1e-30))
            {
              n.setZero();
              n(0) = Real(1);
            }
            else
              n /= n.norm();
            facetMeasures[i] = measure;
            facetNormals[i] = n;
          }
        }
        Real normalJumpSq = 0;
        Real normalJumpWeight = 0;
        Real normalJumpMax = 0;
        if (meshDim == 2)
        {
          std::unordered_map<Index, std::vector<std::size_t>> incident;
          incident.reserve(2 * interfaceFacets.size());
          for (std::size_t i = 0; i < interfaceFacets.size(); ++i)
            for (const Index v : mesh.getFace(interfaceFacets[i])->getVertices())
              incident[v].push_back(i);
          for (const auto& [_, faces] : incident)
          {
            for (std::size_t a = 0; a < faces.size(); ++a)
              for (std::size_t b = a + 1; b < faces.size(); ++b)
              {
                const auto ia = faces[a], ib = faces[b];
                const Real dot = std::abs(facetNormals[ia].dot(facetNormals[ib]));
                const Real theta = std::acos(clampUnit(dot));
                const Real w = Real(0.5) * (facetMeasures[ia] + facetMeasures[ib]);
                normalJumpSq += w * theta * theta;
                normalJumpWeight += w;
                normalJumpMax = std::max(normalJumpMax, theta);
              }
          }
        }
        else if (meshDim == 3)
        {
          std::unordered_map<std::uint64_t, std::vector<std::size_t>> incident;
          incident.reserve(3 * interfaceFacets.size());
          for (std::size_t i = 0; i < interfaceFacets.size(); ++i)
          {
            const auto& vv = mesh.getFace(interfaceFacets[i])->getVertices();
            incident[edgeKey(vv[0], vv[1])].push_back(i);
            incident[edgeKey(vv[1], vv[2])].push_back(i);
            incident[edgeKey(vv[2], vv[0])].push_back(i);
          }
          for (const auto& [_, faces] : incident)
          {
            for (std::size_t a = 0; a < faces.size(); ++a)
              for (std::size_t b = a + 1; b < faces.size(); ++b)
              {
                const auto ia = faces[a], ib = faces[b];
                const Real dot = std::abs(facetNormals[ia].dot(facetNormals[ib]));
                const Real theta = std::acos(clampUnit(dot));
                const Real w = Real(0.5) * (facetMeasures[ia] + facetMeasures[ib]);
                normalJumpSq += w * theta * theta;
                normalJumpWeight += w;
                normalJumpMax = std::max(normalJumpMax, theta);
              }
          }
        }
        rep.normalJumpRMS = normalJumpWeight > Real(0)
          ? std::sqrt(normalJumpSq / normalJumpWeight)
          : Real(0);
        rep.normalJumpMax = normalJumpMax;
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

        auto facetPoint = [&](const FacetTable& ft, const FacetQP& fq) {
          const auto face = mesh.getFace(ft.facet);
          return Geometry::Point(*face, fq.rc, fq.X);
        };

        // ---- Fixed Welsch scale: σ = max(3h, q90(|φ| on Γ at u=0)) ----
        std::vector<Real> r0abs;
        for (const auto& ft : facetTables)
          for (const auto& fq : ft.qps)
          {
            SpatialVec zero = SpatialVec::Zero(meshDim);
            const auto src = facetPoint(ft, fq);
            r0abs.push_back(std::abs(
              Detail::evaluateTranslatedPoint(phi, src, zero, p.pointLocationTolerance)));
          }
        Real sigma = Real(3) * h;
        if (!r0abs.empty())
        {
          const std::size_t k90 =
            static_cast<std::size_t>(Real(0.9) * static_cast<Real>(r0abs.size() - 1));
          std::nth_element(r0abs.begin(), r0abs.begin() + k90, r0abs.end());
          sigma = std::max(sigma, r0abs[k90]);
        }
        const Real sigma2 = sigma * sigma;
        rep.sigma = sigma;
        if (!p.hasInterfaceAttribute)
        {
          rep.exitReason = "missing-interface-attribute";
          m_report = rep;
          return rep;
        }

        struct CellQP
        {
            Real w = 0;
            Math::SpatialPoint rc;
            SpatialVec X;
        };
        struct CellTable
        {
            Index cell = 0;
            std::vector<CellQP> qps;
        };
        std::vector<CellTable> cellTables;
        cellTables.reserve(mesh.getCellCount());
        for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        {
          const Index cellIdx = cellIt->getIndex();
          const auto& fe = fes.getFiniteElement(meshDim, cellIdx);
          const std::size_t feOrder = fe.getOrder();
          const std::size_t qOrder = p.quadratureOrder > 0
            ? p.quadratureOrder
            : std::max<std::size_t>(2, 2 * feOrder);
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(qOrder, cellIt->getGeometry());
          const auto& quad = cellIt->getQuadrature(qf);

          CellTable ct;
          ct.cell = cellIdx;
          ct.qps.resize(quad.getSize());
          for (std::size_t q = 0; q < quad.getSize(); ++q)
          {
            const auto& pt = quad.getPoint(q);
            auto& cq = ct.qps[q];
            cq.w = qf.getWeight(q) * pt.getDistortion();
            cq.rc = pt.getReferenceCoordinates();
            const auto& pc = pt.getPhysicalCoordinates();
            cq.X = SpatialVec::Zero(meshDim);
            for (std::size_t r = 0; r < meshDim; ++r)
              cq.X(static_cast<Eigen::Index>(r)) = pc(r);
          }
          cellTables.push_back(std::move(ct));
        }

        // ============================================================
        // Per-iteration field evaluations through the GridFunction.
        // ============================================================
        auto fieldAtFacetQP = [&](const Displacement& gf, const FacetTable& ft,
                                const FacetQP& fq) -> SpatialVec {
          return gf.getValue(facetPoint(ft, fq));
        };

        auto deformationAtCellQP = [&](const Displacement& gf, const CellTable& ct,
                                     const CellQP& cq) -> SpatialMat {
          const auto cell = mesh.getCell(ct.cell);
          const Geometry::Point pt(*cell, cq.rc, cq.X);
          SpatialMat F = SpatialMat::Identity(meshDim, meshDim);
          F += Variational::Jacobian(gf).getValue(pt);
          return F;
        };

        struct FastAdm
        {
            Real minJ = std::numeric_limits<Real>::infinity();
            Real maxQ = 0;
            std::size_t inadmissibleCount = 0;
        };
        auto fastAdmissibility = [&](const Displacement& gf) {
          FastAdm a;
          for (const auto& ct : cellTables)
            for (const auto& cq : ct.qps)
            {
              const SpatialMat F = deformationAtCellQP(gf, ct, cq);
              const Real j = F.determinant();
              if (j < a.minJ)
                a.minJ = j;
              if (j <= p.jMinRatio)
                ++a.inadmissibleCount;
              if (j > Real(0))
              {
                const Real q = F.squaredNorm() / (d * std::pow(j, Real(2) / d));
                if (q > a.maxQ)
                  a.maxQ = q;
              }
            }
          return a;
        };

        struct SurfaceState
        {
            Real energy = 0;
            Real activeLen = 0, totalLen = 0;
            Real activeRMS = 0, activeSup = 0;
        };
        auto surfaceState = [&](const Displacement& gf) {
          SurfaceState s;
          Real sq = 0;
          for (const auto& ft : facetTables)
            for (const auto& fq : ft.qps)
            {
              const SpatialVec y = fq.X + fieldAtFacetQP(gf, ft, fq);
              const auto src = facetPoint(ft, fq);
              const SpatialVec displacement = y - fq.X;
              const Real r = Detail::evaluateTranslatedPoint(
                phi, src, displacement, p.pointLocationTolerance);
              const Real omega = std::exp(-r * r / sigma2);
              s.totalLen += fq.w;
              s.energy += fq.w * Real(0.5) * sigma2 * (Real(1) - omega);
              if (omega >= p.omegaMin)
              {
                s.activeLen += fq.w;
                sq += fq.w * r * r;
                s.activeSup = std::max(s.activeSup, std::abs(r));
              }
            }
          s.activeRMS = s.activeLen > Real(0) ? std::sqrt(sq / s.activeLen) : Real(0);
          return s;
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
          m_bulkForm = Variational::Integral(gammaM * m_duStep, m_vStep) +
            Variational::Integral(gammaH * ellM * ellM * devTrial, devTest) +
            Variational::Integral(gammaDiv * ellM * ellM * divTrial, divTest);
          m_bulkForm.assemble();
          rep.tBulk = secondsSince(bulkTic);
          m_bulkFormAssembled = true;
        }

        Displacement vK(fes);
        Displacement scratch(fes);
        Displacement previousU(fes);
        Displacement uTrial(fes);

        if (p.initialGuessGamma > Real(0))
        {
          auto tic = Clock::now();
          Detail::WNGIRNormalOffsetMetric initMetric(
            phi, grad, m_duStep, m_vStep, p, sigma2);
          initMetric.over(p.interfaceAttribute);
          Detail::WNGIRNormalOffsetForce initForce(phi, grad, m_vStep, p, sigma2);
          initForce.over(p.interfaceAttribute);
          typename ProblemType::ProblemBodyType body(m_bulkForm);
          body = body + initMetric - initForce;
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
          previousU = u;
          const Real eBase = surfaceState(previousU).energy;
          Real alpha = Real(1);
          bool accepted = false;
          FastAdm adm{};
          Real eTrial = std::numeric_limits<Real>::infinity();
          while (alpha >= p.alphaMin)
          {
            uTrial = vK;
            uTrial *= alpha;
            uTrial += previousU;
            adm = fastAdmissibility(uTrial);
            const bool jOK = adm.inadmissibleCount == 0 && adm.minJ > p.jLineSearchRatio;
            const bool qOK = adm.maxQ < p.qMax;
            bool eOK = true;
            if (jOK && qOK && p.energyLineSearch)
            {
              eTrial = surfaceState(uTrial).energy;
              eOK = std::isfinite(eTrial) && eTrial <= eBase;
            }
            if (jOK && qOK && eOK)
            {
              u = uTrial;
              accepted = true;
              break;
            }
            alpha *= Real(0.5);
          }
          rep.tLineSearch += secondsSince(tic);
          if (p.trace)
          {
            const auto surf = surfaceState(u);
            std::cout << "      wngir init=normal-offset"
                      << "  accepted=" << (accepted ? 1 : 0) << "  alpha=" << alpha
                      << "  actRMS=" << std::scientific << surf.activeRMS
                      << "  actRMS/h=" << (h > Real(0) ? surf.activeRMS / h : Real(0))
                      << "  min_j=" << (accepted ? adm.minJ : rep.minJ)
                      << "  max_Q=" << (accepted ? adm.maxQ : rep.maxQRel) << '\n';
          }
        }

        Real ePrev = surfaceState(u).energy;

        // ============================================================
        // Nonlinear iteration.
        // ============================================================
        for (; rep.iterations < p.maxIterations; ++rep.iterations)
        {
          auto tic = Clock::now();
          Detail::WNGIRAdmissibilityMetric admMetric(m_duStep, m_vStep, u, p);
          Detail::WNGIRSurfaceObservationMetric obsMetric(
            phi, grad, m_duStep, m_vStep, u, p, sigma2);
          obsMetric.over(p.interfaceAttribute);
          Detail::WNGIRSurfaceForce surfaceForce(phi, grad, m_vStep, u, p, sigma2);
          surfaceForce.over(p.interfaceAttribute);

          std::size_t linearIterations = 0;
          Real linearError = std::numeric_limits<Real>::infinity();

          typename ProblemType::ProblemBodyType body(m_bulkForm);
          body = body + obsMetric + admMetric - surfaceForce;
          m_stepProblem = body;
          m_stepProblem.assemble();
          rep.tAssembly += secondsSince(tic);

          tic = Clock::now();
          const bool solveOk = solveStep(vK, linearIterations, linearError);
          rep.tSolve += secondsSince(tic);
          if (!solveOk)
          {
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

          // ---- Optimal step rescale β = ⟨d,v⟩_Γ / ⟨v,v⟩_Γ ----
          Real beta = Real(1);
          Real rawDemonsRMS = Real(0);
          Real liftedTraceRMS = Real(0);
          Real scaledTraceRMS = Real(0);
          if (p.betaMax > Real(1) || p.trace)
          {
            Real bNum = 0, bDen = 0, dDen = 0, gammaMeasure = 0;
            for (const auto& ft : facetTables)
              for (const auto& fq : ft.qps)
              {
                const SpatialVec y = fq.X + fieldAtFacetQP(u, ft, fq);
                const auto src = facetPoint(ft, fq);
                const SpatialVec disp = y - fq.X;
                const Real r = Detail::evaluateTranslatedPoint(
                  phi, src, disp, p.pointLocationTolerance);
                const SpatialVec g = Detail::evaluateTranslatedPoint(
                  grad, src, disp, p.pointLocationTolerance);
                const Real obsWeight = g.dot(g) + epsG +
                  (p.residualStabilizedObservationMetric ? (r * r) / sigma2 : Real(0));
                const Real omega = std::exp(-r * r / sigma2);
                const SpatialVec dVec = (-omega * r / obsWeight) * g;
                const SpatialVec v = fieldAtFacetQP(vK, ft, fq);
                bNum += fq.w * dVec.dot(v);
                bDen += fq.w * v.dot(v);
                dDen += fq.w * dVec.dot(dVec);
                gammaMeasure += fq.w;
              }
            if (gammaMeasure > Real(0))
            {
              rawDemonsRMS = std::sqrt(std::max<Real>(Real(0), dDen) / gammaMeasure);
              liftedTraceRMS = std::sqrt(std::max<Real>(Real(0), bDen) / gammaMeasure);
            }
            if (p.betaMax > Real(1) && bDen > Real(0) && std::isfinite(bNum))
            {
              beta = std::clamp(bNum / bDen, Real(1), p.betaMax);
              vK *= beta;
            }
            scaledTraceRMS = beta * liftedTraceRMS;
          }

          const Real maxStep = std::max(std::abs(vK.max()), std::abs(vK.min()));
          if (maxStep <= stepTol)
          {
            rep.exitReason = "best-effort-step-stagnation";
            break;
          }

          // ---- Nonlinear line search on TRUE geometry ----
          tic = Clock::now();
          previousU = u;
          Real alpha = Real(1);
          bool accepted = false;
          std::size_t backtracks = 0;
          FastAdm adm{};
          Real eTrial = std::numeric_limits<Real>::infinity();
          while (alpha >= p.alphaMin)
          {
            // uTrial = previousU + alpha * vK
            uTrial = vK;
            uTrial *= alpha;
            uTrial += previousU;
            adm = fastAdmissibility(uTrial);
            const bool jOK = adm.inadmissibleCount == 0 && adm.minJ > p.jLineSearchRatio;
            const bool qOK = adm.maxQ < p.qMax;
            bool eOK = true;
            if (jOK && qOK && p.energyLineSearch)
            {
              eTrial = surfaceState(uTrial).energy;
              eOK = std::isfinite(eTrial) && eTrial <= ePrev;
            }
            if (jOK && qOK && eOK)
            {
              u = uTrial;
              accepted = true;
              break;
            }
            alpha *= Real(0.5);
            ++backtracks;
          }
          rep.tLineSearch += secondsSince(tic);
          if (!accepted)
          {
            u = previousU;
            rep.exitReason = "line-search-failure";
            break;
          }

          FastAdm acceptedAdm = adm;
          SurfaceState acceptedSurf = surfaceState(u);
          Real acceptedEnergy =
            p.energyLineSearch && std::isfinite(eTrial) ? eTrial : acceptedSurf.energy;

          rep.lastAlpha = alpha;
          // acceptedStep = max_i |u_i - previousU_i|
          scratch = u;
          scratch -= previousU;
          rep.acceptedStep = std::max(std::abs(scratch.max()), std::abs(scratch.min()));
          rep.minJ = acceptedAdm.minJ;
          rep.maxQRel = acceptedAdm.maxQ;

          const auto surf = acceptedSurf;
          const Real eNow = acceptedEnergy;
          rep.activeRMS = surf.activeRMS;
          rep.activeSup = surf.activeSup;
          rep.activeFraction =
            surf.totalLen > Real(0) ? surf.activeLen / surf.totalLen : Real(0);
          rep.energy = eNow;
          if (p.trace)
            std::cout << "      wngir it=" << std::setw(3) << rep.iterations
                      << "  E=" << std::scientific << std::setprecision(3) << eNow
                      << "  actRMS=" << surf.activeRMS
                      << "  actRMS/h=" << (h > Real(0) ? surf.activeRMS / h : Real(0))
                      << "  actSup=" << surf.activeSup
                      << "  actFrac=" << rep.activeFraction
                      << "  dΓ/h=" << (h > Real(0) ? rawDemonsRMS / h : Real(0))
                      << "  vΓ/h=" << (h > Real(0) ? liftedTraceRMS / h : Real(0))
                      << "  β=" << beta
                      << "  step/h=" << (h > Real(0) ? rep.acceptedStep / h : Real(0))
                      << "  linIt=" << rep.linearIterations << "  α=" << alpha
                      << "  bt=" << backtracks << "  min_j=" << rep.minJ
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
          if (activeRMSOverHTol > Real(0) && h > Real(0) &&
            surf.activeRMS / h <= activeRMSOverHTol)
          {
            rep.exitReason = "geometric-rms-converged";
            ++rep.iterations;
            break;
          }
          if (activeSupOverHTol > Real(0) && h > Real(0) &&
            surf.activeSup / h <= activeSupOverHTol)
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
        m_report = rep;
        return rep;
      }

    private:
      // Solves the currently-assembled step problem with CG and copies the
      // backend-matched solution GridFunction into @p out.
      bool solveStep(Displacement& out, std::size_t& iterations, Real& error)
      {
        m_stepProblem.solve(m_stepSolver);
        out = m_duStep.getSolution();
        iterations = 0;
        error = Real(0);
        return std::isfinite(out.max()) && std::isfinite(out.min());
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
