/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRSOLVE_H
#define RODIN_ADAPTATION_WNGIRSOLVE_H

#include "WNGIRAdmissibilityGradient.h"
#include "WNGIRAdmissibilityMetric.h"
#include "WNGIRReport.h"
#include "WNGIRSurfaceForce.h"
#include "WNGIRSurfaceObservationMetric.h"

namespace Rodin::Adaptation
{
  /// Backend-neutral reusable-scratch overload. The backend-native
  /// trial/test functions, preassembled bulk form, step Problem, and linear
  /// solver are caller-owned and persist across iterations and frames.
  template <class Mesh, class Displacement, class PhiDerived, class GradDerived,
            class TrialT, class TestT, class BilinearFormT, class ProblemT,
            class LinearSolverT>
  WNGIRReport solveWNGIR(
      const Mesh& mesh,
      Displacement& u,
      const std::vector<Rodin::Index>& interfaceFacets,
      const Variational::RealFunctionBase<PhiDerived>& phi,
      const Variational::VectorFunctionBase<Real, GradDerived>& grad,
      const WNGIRParameters& p,
      TrialT& duStep,
      TestT& vStep,
      BilinearFormT& bulkForm,
      ProblemT& stepProblem,
      LinearSolverT& linearSolver,
      bool& bulkFormAssembled)
  {
    using Vec = Math::SpatialVector<Real>;
    using Mat = Math::SpatialMatrix<Real>;
    using Rodin::Index;
    auto zeroVec = [](std::size_t dim) {
      Vec out(dim);
      out.setZero();
      return out;
    };
    auto zeroMat = [](std::size_t dim) {
      Mat out(dim, dim);
      out.setZero();
      return out;
    };

    WNGIRReport rep;
    const Real h = p.h;
    const Real gammaM = p.gammaM > Real(0) ? p.gammaM : Real(1) / h;
    const Real gammaH = p.gammaH > Real(0) ? p.gammaH : Real(1) / h;
    const Real gammaDiv = p.gammaDiv > Real(0) ? p.gammaDiv : gammaH;
    const Real ellM = p.ellM > Real(0) ? p.ellM : Real(3) * h;
    const Real activeRMSTol =
      p.activeRMSTol > Real(0) ? p.activeRMSTol : Real(4) * h * h;
    const Real activeSupTol =
      p.activeSupTol > Real(0) ? p.activeSupTol : Real(10) * h * h;
    const Real activeRMSOverHTol = p.activeRMSOverHTol;
    const Real activeSupOverHTol = p.activeSupOverHTol;
    const Real stepTol =
      p.stepTol > Real(0) ? p.stepTol : Real(1e-4) * h;
    constexpr Real epsG = Real(1e-12);

    const auto& fes = u.getFiniteElementSpace();
    const std::size_t meshDim = mesh.getDimension();
    const Real d = static_cast<Real>(meshDim);

    // =================================================================
    // Per-frame tabulation (geometry is fixed; only u changes).
    // =================================================================

    // ---- Facet tables: trace basis at facet quadrature points ----
    struct FacetQP
    {
      Real w = 0;                 ///< quadrature weight × |J_facet|.
      Vec X;                      ///< physical coordinates.
      Math::SpatialPoint rc;      ///< reference coordinates on the face.
      std::vector<Vec> val;       ///< trace basis values (per local).
    };
    struct FacetTable
    {
      Index facet = 0;
      std::vector<Index> dofs;
      std::vector<FacetQP> qps;
    };
    std::vector<FacetTable> facetTables;
    facetTables.reserve(interfaceFacets.size());
    for (const Index facetIdx : interfaceFacets)
    {
      const auto face = mesh.getFace(facetIdx);
      const auto& fe = fes.getFiniteElement(meshDim - 1, facetIdx);
      const std::size_t nLocal = fe.getCount();
      const std::size_t feOrder = fe.getOrder();
      const std::size_t qOrder = p.quadratureOrder > 0
        ? p.quadratureOrder
        : std::max<std::size_t>(2, 2 * feOrder);
      const auto& qf = QF::PolytopeQuadratureFormula::get(
          qOrder, face->getGeometry());
      const auto& quad = face->getQuadrature(qf);

      FacetTable ft;
      ft.facet = facetIdx;
      ft.dofs.resize(nLocal);
      for (std::size_t l = 0; l < nLocal; ++l)
        ft.dofs[l] = fes.getGlobalIndex({meshDim - 1, facetIdx}, l);
      ft.qps.resize(quad.getSize());
      for (std::size_t q = 0; q < quad.getSize(); ++q)
      {
        const auto& pt = quad.getPoint(q);
        auto& fq = ft.qps[q];
        fq.w = qf.getWeight(q) * pt.getDistortion();
        const auto& pc = pt.getPhysicalCoordinates();
        fq.X = zeroVec(meshDim);
        for (std::size_t r = 0; r < meshDim; ++r)
          fq.X(static_cast<Eigen::Index>(r)) = pc(r);
        const auto& rc = pt.getReferenceCoordinates();
        fq.rc = rc;
        fq.val.resize(nLocal);
        for (std::size_t l = 0; l < nLocal; ++l)
        {
          const auto bv = fe.getBasis(l)(rc);
          fq.val[l] = zeroVec(meshDim);
          for (std::size_t r = 0; r < meshDim; ++r)
            fq.val[l](static_cast<Eigen::Index>(r)) = bv(r);
        }
      }
      facetTables.push_back(std::move(ft));
    }

    auto sourcePoint =
      [&](const FacetTable& ft, const FacetQP& fq)
    {
      const auto face = mesh.getFace(ft.facet);
      return Geometry::Point(*face, fq.rc, fq.X);
    };

    // ---- Fixed Welsch scale: σ = max(3h, q90(|φ| on Γ at u=0)) ----
    std::vector<Real> r0abs;
    for (const auto& ft : facetTables)
      for (const auto& fq : ft.qps)
      {
        Vec zero = zeroVec(meshDim);
        const auto src = sourcePoint(ft, fq);
        r0abs.push_back(std::abs(
              Detail::evaluateTranslatedPoint(
                phi, src, zero, p.pointLocationTolerance)));
      }
    Real sigma = Real(3) * h;
    if (!r0abs.empty())
    {
      const std::size_t k90 = static_cast<std::size_t>(
          Real(0.9) * static_cast<Real>(r0abs.size() - 1));
      std::nth_element(r0abs.begin(), r0abs.begin() + k90, r0abs.end());
      sigma = std::max(sigma, r0abs[k90]);
    }
    const Real sigma2 = sigma * sigma;
    rep.sigma = sigma;
    if (!p.hasInterfaceAttribute)
    {
      rep.exitReason = "missing-interface-attribute";
      return rep;
    }

    // ---- Cell tables: validation geometry for line search/diagnostics.
    struct CellQP
    {
      Real w = 0;                  ///< quadrature weight × |J_cell|.
      std::vector<Mat> jac;        ///< physical basis Jacobians.
    };
    struct CellTable
    {
      std::vector<Index> dofs;
      std::vector<CellQP> qps;
    };
    std::vector<CellTable> cellTables;
    cellTables.reserve(mesh.getCellCount());
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const Index cellIdx = cellIt->getIndex();
      const auto& fe = fes.getFiniteElement(meshDim, cellIdx);
      const std::size_t nLocal = fe.getCount();
      const std::size_t feOrder = fe.getOrder();
      const std::size_t qOrder = p.quadratureOrder > 0
        ? p.quadratureOrder
        : std::max<std::size_t>(2, 2 * feOrder);
      const auto& qf = QF::PolytopeQuadratureFormula::get(
          qOrder, cellIt->getGeometry());
      const auto& quad = cellIt->getQuadrature(qf);

      CellTable ct;
      ct.dofs.resize(nLocal);
      for (std::size_t l = 0; l < nLocal; ++l)
        ct.dofs[l] = fes.getGlobalIndex({meshDim, cellIdx}, l);

      ct.qps.resize(quad.getSize());
      for (std::size_t q = 0; q < quad.getSize(); ++q)
      {
        const auto& pt = quad.getPoint(q);
        auto& cq = ct.qps[q];
        cq.w = qf.getWeight(q) * pt.getDistortion();
        const auto& rc = pt.getReferenceCoordinates();
        // Physical Jacobians: phys = ref · Jinv, with the (affine)
        // cell Jacobian inverse from the point.
        const auto Jinv = pt.getJacobianInverse();
        cq.jac.resize(nLocal);
        for (std::size_t l = 0; l < nLocal; ++l)
        {
          const auto jref = fe.getBasis(l).getJacobian()(rc);
          Mat jp = zeroMat(meshDim);
          for (std::size_t r = 0; r < meshDim; ++r)
            for (std::size_t c = 0; c < meshDim; ++c)
              for (std::size_t a = 0; a < meshDim; ++a)
                jp(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c))
                  += jref(r, a) * Jinv(a, c);
          cq.jac[l] = jp;
        }
      }
      cellTables.push_back(std::move(ct));
    }

    // =================================================================
    // Per-iteration helpers (pure contractions of the tables).
    // =================================================================

    // u (or any coefficient vector) evaluated at a facet QP.
    auto fieldAtFacetQP =
      [&](const Math::Vector<Real>& coef,
          const FacetTable& ft, const FacetQP& fq) -> Vec
    {
      Vec out = zeroVec(meshDim);
      for (std::size_t l = 0; l < ft.dofs.size(); ++l)
        out += coef(ft.dofs[l]) * fq.val[l];
      return out;
    };

    // F = I + ∇u at a cell QP.
    auto deformationAtCellQP =
      [&](const Math::Vector<Real>& coef,
          const CellTable& ct, const CellQP& cq) -> Mat
    {
      Mat F = Mat::Identity(meshDim, meshDim);
      for (std::size_t l = 0; l < ct.dofs.size(); ++l)
      {
        const Real c = coef(ct.dofs[l]);
        F += c * cq.jac[l];
      }
      return F;
    };

    // Closed-form admissibility over all validation points.
    struct FastAdm
    {
      Real minJ = std::numeric_limits<Real>::infinity();
      Real maxQ = 0;
      std::size_t inadmissibleCount = 0;
    };
    auto fastAdmissibility = [&](const Math::Vector<Real>& uData)
    {
      FastAdm a;
      for (const auto& ct : cellTables)
        for (const auto& cq : ct.qps)
        {
          const Mat F = deformationAtCellQP(uData, ct, cq);
          const Real j = F.determinant();
          if (j < a.minJ) a.minJ = j;
          if (j <= p.jMinRatio) ++a.inadmissibleCount;
          if (j > Real(0))
          {
            const Real q =
              F.squaredNorm() / (d * std::pow(j, Real(2) / d));
            if (q > a.maxQ) a.maxQ = q;
          }
        }
      return a;
    };

    // Welsch surface energy + active-set diagnostics on Γ_ψ,h.
    struct SurfaceState
    {
      Real energy = 0;
      Real activeLen = 0, totalLen = 0;
      Real activeRMS = 0, activeSup = 0;
    };
    auto surfaceState = [&](const Math::Vector<Real>& uData)
    {
      SurfaceState s;
      Real sq = 0;
      for (const auto& ft : facetTables)
        for (const auto& fq : ft.qps)
        {
          const Vec y = fq.X + fieldAtFacetQP(uData, ft, fq);
          const auto src = sourcePoint(ft, fq);
          const Vec displacement = y - fq.X;
          const Real r =
            Detail::evaluateTranslatedPoint(
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
      s.activeRMS = s.activeLen > Real(0)
        ? std::sqrt(sq / s.activeLen) : Real(0);
      return s;
    };

    Real ePrev = surfaceState(u.getData()).energy;
    using Clock = std::chrono::steady_clock;
    auto secondsSince = [](Clock::time_point t0) -> Real
    {
      return std::chrono::duration<Real>(Clock::now() - t0).count();
    };

    struct AndersonState
    {
      std::size_t memory = 0;
      std::vector<Math::Vector<Real>> u;
      std::vector<Math::Vector<Real>> g;
      std::vector<Math::Vector<Real>> f;

      void push(
          const Math::Vector<Real>& uk,
          const Math::Vector<Real>& gk)
      {
        u.push_back(uk);
        g.push_back(gk);
        f.push_back(gk - uk);
        while (u.size() > memory + 1)
        {
          u.erase(u.begin());
          g.erase(g.begin());
          f.erase(f.begin());
        }
      }

      bool canAccelerate() const
      {
        return memory > 0 && f.size() >= 2;
      }

      Math::Vector<Real> candidate() const
      {
        const std::size_t p = std::min<std::size_t>(memory, f.size() - 1);
        const std::size_t start = f.size() - p - 1;
        const Eigen::Index n = f.back().size();
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> dF(n, p);
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> dG(n, p);
        for (std::size_t j = 0; j < p; ++j)
        {
          dF.col(static_cast<Eigen::Index>(j)) =
            f[start + j + 1] - f[start + j];
          dG.col(static_cast<Eigen::Index>(j)) =
            g[start + j + 1] - g[start + j];
        }
        const Eigen::Matrix<Real, Eigen::Dynamic, 1> gamma =
          dF.colPivHouseholderQr().solve(f.back());
        return g.back() - dG * gamma;
      }
    };

    AndersonState anderson;
    anderson.memory = p.andersonMemory;

    // =================================================================
    // Constant bulk metric, pre-assembled once per frame.
    // =================================================================
    // M_bulk = γ_M ∫ v·z
    //        + ℓ² ∫ γ_dev dev eps(v):dev eps(z)
    //              + γ_div div(v) div(z)
    // is independent of u, so it is assembled and Dirichlet-eliminated once
    // and reused across all iterations. This first-derivative metric is
    // FES-independent and distributes interface motion by penalizing shear
    // and compression modes rather than raw gradient components.
    if (!bulkFormAssembled)
    {
      const auto epsTrial = Real(0.5)
        * (Variational::Jacobian(duStep)
           + Variational::Jacobian(duStep).T());
      const auto epsTest = Real(0.5)
        * (Variational::Jacobian(vStep)
           + Variational::Jacobian(vStep).T());
      const auto divTrial = Variational::Trace(epsTrial);
      const auto divTest = Variational::Trace(epsTest);
      const auto devTrial =
        epsTrial - (Real(1) / d) * divTrial * Variational::IdentityMatrix(meshDim);
      const auto devTest =
        epsTest - (Real(1) / d) * divTest * Variational::IdentityMatrix(meshDim);
      bulkForm =
          Variational::Integral(gammaM * duStep, vStep)
        + Variational::Integral(
            gammaH * ellM * ellM * devTrial, devTest)
        + Variational::Integral(
            gammaDiv * ellM * ellM * divTrial, divTest);
      bulkForm.assemble();
      bulkFormAssembled = true;
    }

    // =================================================================
    // Nonlinear iteration.
    // =================================================================
    for (; rep.iterations < p.maxIterations; ++rep.iterations)
    {
      // ---- One variational assembly: bulk + surface + admissibility ----
      auto tic = Clock::now();
      auto zeroBoundary = Variational::VectorFunction(
          meshDim,
          [&](const Geometry::Point&)
          {
            return zeroVec(meshDim);
          });

      Detail::WNGIRAdmissibilityMetric admMetric(duStep, vStep, u, p);
      Detail::WNGIRAdmissibilityGradient admGradient(vStep, u, p);
      Detail::WNGIRSurfaceObservationMetric obsMetric(
          phi, grad, duStep, vStep, u, p, sigma2);
      obsMetric.over(p.interfaceAttribute);
      Detail::WNGIRSurfaceForce surfaceForce(
          phi, grad, vStep, u, p, sigma2);
      surfaceForce.over(p.interfaceAttribute);

      // Only the u-dependent terms are reassembled here. The main fit RHS is
      // Welsch + E_qual by default. Soft quality terms enter the metric only
      // when includeQualityMetric is enabled. If includeQualityGradient is
      // disabled, the main RHS is Welsch-only. In split mode a second solve
      // uses the quality first variation as a
      // quality-recovery RHS, and the final direction is v_fit + eta v_qual.
      const bool useGradientForce =
           p.includeAdmissibilityGradient
        || p.includeQualityGradient;
      const bool useSplitQuality =
        p.splitQualityDirection
        && (p.gammaQual > Real(0)
            || p.gammaJ > Real(0)
            || p.gammaSize > Real(0));

      Math::Vector<Real> vK;
      std::size_t linearIterations = 0;
      Real linearError = std::numeric_limits<Real>::infinity();

      if (useSplitQuality)
      {
        typename ProblemT::ProblemBodyType body(bulkForm);
        body = body + obsMetric + admMetric - surfaceForce
             + Variational::DirichletBC(duStep, zeroBoundary);
        stepProblem = body;
        stepProblem.assemble();
        rep.tAssembly += secondsSince(tic);

        tic = Clock::now();
        Math::Vector<Real> vFit;
        std::size_t fitIterations = 0;
        Real fitError = std::numeric_limits<Real>::infinity();
        if (!linearSolver.solve(
              stepProblem, vFit, fitIterations, fitError, p))
        {
          rep.exitReason = "solve-fit-linear-failed";
          break;
        }
        rep.tSolve += secondsSince(tic);

        tic = Clock::now();
        typename ProblemT::ProblemBodyType qualityBody(bulkForm);
        qualityBody = qualityBody + obsMetric + admMetric - admGradient
                    + Variational::DirichletBC(duStep, zeroBoundary);
        stepProblem = qualityBody;
        stepProblem.assemble();
        rep.tAssembly += secondsSince(tic);

        tic = Clock::now();
        Math::Vector<Real> vQual;
        std::size_t qualIterations = 0;
        Real qualError = std::numeric_limits<Real>::infinity();
        if (!linearSolver.solve(
              stepProblem, vQual, qualIterations, qualError, p))
        {
          rep.exitReason = "solve-quality-linear-failed";
          break;
        }
        rep.tSolve += secondsSince(tic);

        const Real fitNorm2 =
          linearSolver.metricDot(stepProblem, vFit, vFit);
        if (fitNorm2 > Real(1e-30) && std::isfinite(fitNorm2))
        {
          const Real qualFit =
            linearSolver.metricDot(stepProblem, vQual, vFit);
          if (std::isfinite(qualFit))
            vQual -= (qualFit / fitNorm2) * vFit;
        }
        vK = vFit + p.qualityDirectionWeight * vQual;
        linearIterations = fitIterations + qualIterations;
        linearError = std::max(fitError, qualError);
      }
      else
      {
        typename ProblemT::ProblemBodyType body(bulkForm);
        if (useGradientForce)
        {
          body = body + obsMetric + admMetric - surfaceForce - admGradient
               + Variational::DirichletBC(duStep, zeroBoundary);
        }
        else
        {
          body = body + obsMetric + admMetric - surfaceForce
               + Variational::DirichletBC(duStep, zeroBoundary);
        }
        stepProblem = body;
        stepProblem.assemble();
        rep.tAssembly += secondsSince(tic);

        // ---- Solve the single backend-native assembled Problem ----
        tic = Clock::now();
        if (!linearSolver.solve(
              stepProblem, vK, linearIterations, linearError, p))
        {
          rep.exitReason = "solve-linear-failed";
          break;
        }
        rep.tSolve += secondsSince(tic);
      }
      rep.linearIterations += linearIterations;
      rep.linearError = linearError;
      if (!vK.allFinite())
      {
        rep.exitReason = "solve-nonfinite";
        break;
      }

      // ---- 6. Optimal step rescale β = ⟨d,v⟩_Γ/⟨v,v⟩_Γ ----
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
            const Vec y = fq.X + fieldAtFacetQP(u.getData(), ft, fq);
            const auto src = sourcePoint(ft, fq);
            const Vec displacement = y - fq.X;
            const Real r =
              Detail::evaluateTranslatedPoint(
                  phi, src, displacement, p.pointLocationTolerance);
            const Vec g =
              Detail::evaluateTranslatedPoint(
                  grad, src, displacement, p.pointLocationTolerance);
            const Real obsWeight =
                g.dot(g) + epsG
              + (p.residualStabilizedObservationMetric
                  ? (r * r) / sigma2 : Real(0));
            const Real omega = std::exp(-r * r / sigma2);
            const Vec dVec = (-omega * r / obsWeight) * g;
            const Vec v = fieldAtFacetQP(vK, ft, fq);
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

      const Real maxStep = vK.cwiseAbs().maxCoeff();
      if (maxStep <= stepTol)
      {
        rep.exitReason = "best-effort-step-stagnation";
        break;
      }

      // ---- 7. Nonlinear line search on TRUE geometry ----
      tic = Clock::now();
      const Math::Vector<Real> previousU = u.getData();
      Real alpha = Real(1);
      bool accepted = false;
      std::size_t backtracks = 0;
      FastAdm adm{};
      Real eTrial = std::numeric_limits<Real>::infinity();
      while (alpha >= p.alphaMin)
      {
        const Math::Vector<Real> uTrial = previousU + alpha * vK;
        adm = fastAdmissibility(uTrial);
        const bool jOK =
          adm.inadmissibleCount == 0
          && adm.minJ > p.jLineSearchRatio;
        const bool qOK = adm.maxQ < p.qMax;
        bool eOK = true;
        if (jOK && qOK && p.energyLineSearch)
        {
          eTrial = surfaceState(uTrial).energy;
          eOK = std::isfinite(eTrial) && eTrial <= ePrev;
        }
        if (jOK && qOK && eOK)
        {
          u.getData() = uTrial;
          accepted = true;
          break;
        }
        alpha *= Real(0.5);
        ++backtracks;
      }
      rep.tLineSearch += secondsSince(tic);
      if (!accepted)
      {
        u.getData() = previousU;
        rep.exitReason = "line-search-failure";
        break;
      }

      Math::Vector<Real> acceptedU = u.getData();
      FastAdm acceptedAdm = adm;
      SurfaceState acceptedSurf = surfaceState(acceptedU);
      Real acceptedEnergy =
        p.energyLineSearch && std::isfinite(eTrial)
          ? eTrial : acceptedSurf.energy;

      // ---- 8. Safeguarded Anderson acceleration of the accepted map ----
      bool aaTried = false;
      bool aaAccepted = false;
      Real aaTheta = Real(0);
      anderson.push(previousU, acceptedU);
      if (anderson.canAccelerate()
          && rep.iterations + 1 >= p.andersonStart)
      {
        aaTried = true;
        ++rep.andersonTried;
        const Math::Vector<Real> aaFull = anderson.candidate();
        const Math::Vector<Real> aaDelta = aaFull - acceptedU;
        Real theta = std::clamp(
            p.andersonDamping, Real(0), Real(1));
        const Real thetaMin = std::clamp(
            p.andersonMinDamping, Real(0), theta);
        while (theta >= thetaMin && theta > Real(0))
        {
          const Math::Vector<Real> aaTrial =
            acceptedU + theta * aaDelta;
          const FastAdm aaAdm = fastAdmissibility(aaTrial);
          const bool jOK =
            aaAdm.inadmissibleCount == 0
            && aaAdm.minJ > p.jLineSearchRatio;
          const bool qOK = aaAdm.maxQ < p.qMax;
          if (jOK && qOK)
          {
            const SurfaceState aaSurf = surfaceState(aaTrial);
            const bool eOK =
              !p.energyLineSearch
              || (std::isfinite(aaSurf.energy)
                  && aaSurf.energy <= acceptedEnergy);
            if (eOK)
            {
              acceptedU = aaTrial;
              acceptedAdm = aaAdm;
              acceptedSurf = aaSurf;
              acceptedEnergy = aaSurf.energy;
              aaAccepted = true;
              aaTheta = theta;
              u.getData() = acceptedU;
              ++rep.andersonAccepted;
              break;
            }
          }
          theta *= Real(0.5);
        }
      }

      rep.lastAlpha = alpha;
      rep.lastAndersonTheta = aaTheta;
      rep.acceptedStep = (acceptedU - previousU).cwiseAbs().maxCoeff();
      rep.minJ = acceptedAdm.minJ;
      rep.maxQRel = acceptedAdm.maxQ;

      // ---- Diagnostics + stopping ----
      const auto surf = acceptedSurf;
      const Real eNow = acceptedEnergy;
      rep.activeRMS = surf.activeRMS;
      rep.activeSup = surf.activeSup;
      rep.activeFraction = surf.totalLen > Real(0)
        ? surf.activeLen / surf.totalLen : Real(0);
      rep.energy = eNow;
      if (p.trace)
        std::cout << "      wngir it=" << std::setw(3) << rep.iterations
                  << "  E=" << std::scientific << std::setprecision(3)
                  << eNow
                  << "  actRMS=" << surf.activeRMS
                  << "  actRMS/h=" << (h > Real(0) ? surf.activeRMS / h : Real(0))
                  << "  actSup=" << surf.activeSup
                  << "  actFrac=" << rep.activeFraction
                  << "  dΓ/h=" << (h > Real(0) ? rawDemonsRMS / h : Real(0))
                  << "  vΓ/h=" << (h > Real(0) ? liftedTraceRMS / h : Real(0))
                  << "  β=" << beta
                  << "  βvΓ/h=" << (h > Real(0) ? scaledTraceRMS / h : Real(0))
                  << "  step/h=" << (h > Real(0) ? rep.acceptedStep / h : Real(0))
                  << "  aa=" << (aaAccepted ? 1 : 0)
                  << "  aaTry=" << (aaTried ? 1 : 0)
                  << "  aaθ=" << aaTheta
                  << "  linIt=" << rep.linearIterations
                  << "  linErr=" << rep.linearError
                  << "  qE=" << (p.includeQualityGradient ? 1 : 0)
                  << "  qM=" << (p.includeQualityMetric ? 1 : 0)
                  << "  admM=" << (p.includeAdmissibilityMetric ? 1 : 0)
                  << "  splitQ=" << (useSplitQuality ? 1 : 0)
                  << "  α=" << alpha
                  << "  bt=" << backtracks
                  << "  min_j=" << rep.minJ
                  << "  max_Q=" << rep.maxQRel
                  << '\n';

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
      if (activeRMSOverHTol > Real(0) && h > Real(0)
          && surf.activeRMS / h <= activeRMSOverHTol)
      {
        rep.exitReason = "geometric-rms-converged";
        ++rep.iterations;
        break;
      }
      if (activeSupOverHTol > Real(0) && h > Real(0)
          && surf.activeSup / h <= activeSupOverHTol)
      {
        rep.exitReason = "geometric-sup-converged";
        ++rep.iterations;
        break;
      }
      const Real eRel =
        std::abs(ePrev - eNow) / std::max(ePrev, Real(1e-30));
      if (eRel < p.energyStagTol)
      {
        rep.exitReason = "best-effort-energy-stagnation";
        ++rep.iterations;
        break;
      }
      ePrev = eNow;
    }
    return rep;
  }
}

#endif
