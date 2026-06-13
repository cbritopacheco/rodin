/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXAMPLES_GEOMETRY_WNGIR_H
#define RODIN_EXAMPLES_GEOMETRY_WNGIR_H

//
// Welsch natural-gradient interface registration (WNGIR).
//
// Classify-then-displace: the classifier provides a topological-
// reference skeleton Γ_ψ,h; the target geometry is the level set φ
// with ∇φ; the displacement u ∈ V_{h,0} (vector Pk, zero trace on ∂Ω)
// moves the skeleton onto φ = 0.
//
// Per frame (cold start u = 0, no warm start, no σ continuation):
//   σ = max(3h, quantile_90(|φ| on Γ_ψ,h at u = 0)),  fixed.
//
// Per nonlinear iteration k:
//   1. Welsch data at skeleton quadrature points X_q:
//        y_q = X_q + u^k(X_q),  r_q = φ(y_q),  g_q = ∇φ(y_q),
//        ω_q = exp(−r_q²/σ²).
//   2. Assemble the negative first variation
//        FΓ(v) = −∫_Γ ω_q r_q g_q·trace(v) dS.
//   3. H¹ bulk metric plus residual-stabilized observation metric
//        M(v,z) = γ_M ∫ v·z + γ_H ℓ_M² ∫ ∇v:∇z.
//        Mobs(v,z) = γ_obs ∫_Γ (|g|² + ε_g + r²/σ²) v·z dS.
//   4. Linearised admissibility at quadrature/validation points (K,q):
//        F = I + ∇u^k,  j = det F,  Q = ‖F‖²_F/(d·det F^{2/d}),
//        a_j(v) = j F^{-T}:∇v,
//        dQdF = (2/d) j^{−2/d}(F − (F:F)/d · F^{-T}),
//        a_Q(v) = dQdF:∇v,
//        s_j(v) = j − j_safe + a_j(v),  s_Q(v) = Q_max − Q − a_Q(v).
//   5. One-SPD active log-barrier metric:
//        B(s) = −log(s/s0) + s/s0 − 1 on (0, s0), else 0,
//        B'(s) = −1/s + 1/s0,  B''(s) = 1/s²,
//      assembled as Gauss–Newton pullback:
//        (M + Mobs + K_adm) v = FΓ.
//   6. Optimal 1-D step rescale β = ⟨d,v⟩_Γ/⟨v,v⟩_Γ ∈ [1, βmax].
//   7. Nonlinear line search on the TRUE geometry:
//      accept the largest α ∈ (0,1] with j > j_ls, Q < Q_max at all
//      validation points; optionally also E_Welsch decrease.
//   8. u^{k+1} = u^k + α β v.
//   9. Stop on: small step, energy stagnation, active RMS/sup below
//      tolerance, iteration budget, or line-search failure.
//
// Order-generic: the displacement space may be vector P1 or vector
// H1 of order k ≥ 2. All basis values/gradients are queried through
// the FES finite-element API and tabulated once per frame at cell and
// facet quadrature points; per-iteration work is pure coefficient
// contractions. The mesh geometry itself is straight (affine cells);
// curved isoparametric geometry would require per-point Jacobians in
// the tabulation pass (single localised change).
//
// Welsch is for partial matching / topology mismatch: the active
// subset (ω ≥ ω_min) has geometric control through the level-set
// residual; the inactive subset is topological-outlier skeleton. The
// residual is φ(X + u(X)) on Γ_ψ,h — ψ never enters.
//

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <Rodin/Math.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Types.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace Rodin::Examples
{
  struct WNGIRParameters
  {
    Real h = 0;                 ///< reference mesh size (required).
    Real gammaM = 0;            ///< L² weight; ≤0 ⇒ 1/h.
    Real gammaH = 0;            ///< H¹ weight; ≤0 ⇒ 1/h.
    Real ellM = 0;              ///< Sobolev length; ≤0 ⇒ 3h.
    Real gammaObs = 1;          ///< surface observation metric weight.
    bool residualStabilizedObservationMetric = true;
    Real gammaJ = 1;            ///< j-barrier weight.
    Real gammaQ = 1;            ///< Q-barrier weight.
    Real jSafe = 1e-2;          ///< barrier floor on j.
    Real qMax = 10;             ///< barrier + line-search ceiling on Q.
    Real s0J = 0.25;            ///< j-barrier activation width.
    Real s0Q = 2;               ///< Q-barrier activation width.
    Real omegaMin = 0.1;        ///< active-set threshold on ω.
    Real alphaMin = 1e-4;       ///< line-search floor.
    bool energyLineSearch = true;
    Real jMinRatio = 1e-8;      ///< hard inadmissibility floor.
    Real jLineSearchRatio = 1e-2;
    Real activeRMSTol = 0;      ///< ≤0 ⇒ 4h².
    Real activeSupTol = 0;      ///< ≤0 ⇒ 10h².
    Real energyStagTol = 1e-4;
    Real stepTol = 0;           ///< ≤0 ⇒ 1e-4·h.
    std::size_t maxIterations = 200;
    std::size_t quadratureOrder = 0; ///< 0 ⇒ 2·(FE order).
    bool trace = false;
    /// If true, also add the nonlinear barrier first variation to the
    /// RHS. Default false: admissibility enters as metric stiffness and
    /// true-geometry line-search validation, not as a competing energy.
    bool includeAdmissibilityGradient = false;
    /// Optimal 1-D rescale of the lifted step along itself:
    ///   β = ⟨d, v⟩_Γ / ⟨v, v⟩_Γ  (surface inner products),
    /// clamped to [1, betaMax]; line search starts at β·v instead of
    /// v. The H¹ lift systematically under-scales the skeleton trace
    /// (gain ≈ surface-weight / M-diagonal ≈ 1/20 at default γ), so
    /// without β the iteration is linearly convergent with ρ ≈ 0.95.
    /// β recovers Newton-matched magnitude while preserving the lift's
    /// smooth admissibility-aware shape. Since β only scales the same
    /// descent direction, the nonlinear line search remains the final
    /// admissibility and energy-decrease guard.
    Real betaMax = 50;
  };

  struct WNGIRReport
  {
    std::size_t iterations = 0;
    Real sigma = 0;
    Real lastAlpha = 0;
    Real acceptedStep = 0;
    Real minJ = 1;
    Real maxQRel = 1;
    Real activeRMS = 0;
    Real activeSup = 0;
    Real activeFraction = 0;
    Real energy = 0;
    const char* exitReason = "iter-budget";
    // Wall-clock breakdown (seconds, accumulated over iterations).
    Real tForce = 0;      ///< skeleton force assembly (Step 1+2).
    Real tBarrier = 0;    ///< linearised barrier assembly (Step 4+5).
    Real tFactor = 0;     ///< sparse LDLT factorization.
    Real tSolve = 0;      ///< triangular solves.
    Real tLineSearch = 0; ///< true-geometry admissibility + energy LS.
  };

  /// PhiFn  : (const Math::SpatialVector<Real>&) -> Real
  /// GradFn : (const Math::SpatialVector<Real>&) -> Math::SpatialVector<Real>
  template <class Mesh, class Displacement, class PhiFn, class GradFn>
  WNGIRReport solveWNGIR(
      const Mesh& mesh,
      Displacement& u,
      const std::vector<Rodin::Index>& interfaceFacets,
      PhiFn&& phiAt,
      GradFn&& gradAt,
      const WNGIRParameters& p)
  {
    using Vec2 = Math::SpatialVector<Real>;
    using Mat2 = Math::SpatialMatrix<Real>;
    using Rodin::Index;
    auto vec2 = [](Real x, Real y) {
      Vec2 out(2); out(0) = x; out(1) = y; return out;
    };

    WNGIRReport rep;
    const Real h = p.h;
    const Real gammaM = p.gammaM > Real(0) ? p.gammaM : Real(1) / h;
    const Real gammaH = p.gammaH > Real(0) ? p.gammaH : Real(1) / h;
    const Real ellM = p.ellM > Real(0) ? p.ellM : Real(3) * h;
    const Real activeRMSTol =
      p.activeRMSTol > Real(0) ? p.activeRMSTol : Real(4) * h * h;
    const Real activeSupTol =
      p.activeSupTol > Real(0) ? p.activeSupTol : Real(10) * h * h;
    const Real stepTol =
      p.stepTol > Real(0) ? p.stepTol : Real(1e-4) * h;
    constexpr Real epsG = Real(1e-12);

    const auto& fes = u.getFiniteElementSpace();
    const std::size_t nDofs = static_cast<std::size_t>(u.getData().size());
    const std::size_t meshDim = mesh.getDimension();

    // =================================================================
    // Per-frame tabulation (geometry is fixed; only u changes).
    // =================================================================

    // ---- Facet tables: trace basis at facet quadrature points ----
    struct FacetQP
    {
      Real w = 0;                 ///< quadrature weight × |J_facet|.
      Vec2 X;                     ///< physical coordinates.
      std::vector<Vec2> val;      ///< trace basis values (per local).
    };
    struct FacetTable
    {
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
        fq.X = vec2(pc(0), pc(1));
        const auto& rc = pt.getReferenceCoordinates();
        fq.val.resize(nLocal);
        for (std::size_t l = 0; l < nLocal; ++l)
        {
          const auto bv = fe.getBasis(l)(rc);
          fq.val[l] = vec2(bv(0), bv(1));
        }
      }
      facetTables.push_back(std::move(ft));
    }

    // ---- Fixed Welsch scale: σ = max(3h, q90(|φ| on Γ at u=0)) ----
    std::vector<Real> r0abs;
    for (const auto& ft : facetTables)
      for (const auto& fq : ft.qps)
        r0abs.push_back(std::abs(phiAt(fq.X)));
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

    // ---- Boundary dofs (strong elimination ⇒ u = 0 on ∂Ω) ----
    std::vector<bool> isBoundaryDof(nDofs, false);
    for (auto bIt = mesh.getBoundary(); bIt; ++bIt)
    {
      const Index bIdx = bIt->getIndex();
      const auto& bfe = fes.getFiniteElement(meshDim - 1, bIdx);
      for (std::size_t l = 0; l < bfe.getCount(); ++l)
        isBoundaryDof[fes.getGlobalIndex({meshDim - 1, bIdx}, l)] = true;
    }

    // ---- Cell tables: basis physical Jacobians at quadrature points,
    //      plus the constant M-metric triplets (assembled in the same
    //      pass; basis values are not stored).
    struct CellQP
    {
      Real w = 0;                  ///< quadrature weight × |J_cell|.
      std::vector<Mat2> jac;       ///< physical basis Jacobians.
    };
    struct CellTable
    {
      std::vector<Index> dofs;
      std::vector<CellQP> qps;
    };
    std::vector<CellTable> cellTables;
    cellTables.reserve(mesh.getCellCount());
    std::vector<Eigen::Triplet<Real>> mTriplets;
    const Real ellM2 = ellM * ellM;
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
      std::vector<Vec2> vals(nLocal);
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
          const auto bv = fe.getBasis(l)(rc);
          vals[l] = vec2(bv(0), bv(1));
          const auto jref = fe.getBasis(l).getJacobian()(rc);
          Mat2 jp(2, 2);
          for (int r = 0; r < 2; ++r)
            for (int c = 0; c < 2; ++c)
              jp(r, c) = jref(r, 0) * Jinv(0, c) + jref(r, 1) * Jinv(1, c);
          cq.jac[l] = jp;
        }
        // M-metric contributions at this quadrature point.
        for (std::size_t i = 0; i < nLocal; ++i)
          for (std::size_t j = 0; j < nLocal; ++j)
          {
            const Real mass = gammaM * cq.w * vals[i].dot(vals[j]);
            Real stiff = 0;
            for (int r = 0; r < 2; ++r)
              for (int c = 0; c < 2; ++c)
                stiff += cq.jac[i](r, c) * cq.jac[j](r, c);
            stiff *= gammaH * ellM2 * cq.w;
            mTriplets.emplace_back(ct.dofs[i], ct.dofs[j], mass + stiff);
          }
      }
      cellTables.push_back(std::move(ct));
    }

    // ---- Pre-factor the boundary-eliminated M once per frame ----
    Eigen::SparseMatrix<Real> mSys(nDofs, nDofs);
    {
      std::vector<Eigen::Triplet<Real>> all;
      all.reserve(mTriplets.size() + nDofs);
      for (const auto& t : mTriplets)
        if (!isBoundaryDof[t.row()] && !isBoundaryDof[t.col()])
          all.push_back(t);
      for (std::size_t i = 0; i < nDofs; ++i)
        if (isBoundaryDof[i])
          all.emplace_back(i, i, Real(1));
      mSys.setFromTriplets(all.begin(), all.end());
    }
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>> mLdlt(mSys);
    const bool mLdltOK = (mLdlt.info() == Eigen::Success);

    // =================================================================
    // Per-iteration helpers (pure contractions of the tables).
    // =================================================================

    // u (or any coefficient vector) evaluated at a facet QP.
    auto fieldAtFacetQP =
      [&](const Math::Vector<Real>& coef,
          const FacetTable& ft, const FacetQP& fq) -> Vec2
    {
      Vec2 out = vec2(0, 0);
      for (std::size_t l = 0; l < ft.dofs.size(); ++l)
      {
        out(0) += coef(ft.dofs[l]) * fq.val[l](0);
        out(1) += coef(ft.dofs[l]) * fq.val[l](1);
      }
      return out;
    };

    // F = I + ∇u at a cell QP.
    auto deformationAtCellQP =
      [&](const Math::Vector<Real>& coef,
          const CellTable& ct, const CellQP& cq) -> Mat2
    {
      Mat2 F = Mat2::Identity(2, 2);
      for (std::size_t l = 0; l < ct.dofs.size(); ++l)
      {
        const Real c = coef(ct.dofs[l]);
        for (int r = 0; r < 2; ++r)
          for (int cc = 0; cc < 2; ++cc)
            F(r, cc) += c * cq.jac[l](r, cc);
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
          const Mat2 F = deformationAtCellQP(uData, ct, cq);
          const Real j = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
          if (j < a.minJ) a.minJ = j;
          if (j <= p.jMinRatio) ++a.inadmissibleCount;
          if (j > Real(0))
          {
            const Real q =
              (F(0, 0) * F(0, 0) + F(0, 1) * F(0, 1)
               + F(1, 0) * F(1, 0) + F(1, 1) * F(1, 1))
              / (Real(2) * j);
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
          const Vec2 y = fq.X + fieldAtFacetQP(uData, ft, fq);
          const Real r = phiAt(y);
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

    // =================================================================
    // Nonlinear iteration.
    // =================================================================
    for (; rep.iterations < p.maxIterations; ++rep.iterations)
    {
      // ---- 1+2. Welsch force + observation metric on Γ ----
      auto tic = Clock::now();
      Math::Vector<Real> rhs = Math::Vector<Real>::Zero(nDofs);
      std::vector<Eigen::Triplet<Real>> obsTriplets;
      for (const auto& ft : facetTables)
        for (const auto& fq : ft.qps)
        {
          const Vec2 y = fq.X + fieldAtFacetQP(u.getData(), ft, fq);
          const Real r = phiAt(y);
          const Vec2 g = gradAt(y);
          const Real omega = std::exp(-r * r / sigma2);
          const Real coef = -omega * r;
          const Real obsWeight =
              g.dot(g) + epsG
            + (p.residualStabilizedObservationMetric
                ? (r * r) / sigma2 : Real(0));
          for (std::size_t l = 0; l < ft.dofs.size(); ++l)
            rhs(ft.dofs[l]) +=
              fq.w * coef * (g(0) * fq.val[l](0) + g(1) * fq.val[l](1));
          if (p.gammaObs > Real(0))
            for (std::size_t i = 0; i < ft.dofs.size(); ++i)
              for (std::size_t j = 0; j < ft.dofs.size(); ++j)
                obsTriplets.emplace_back(
                    ft.dofs[i], ft.dofs[j],
                    p.gammaObs * fq.w * obsWeight
                      * fq.val[i].dot(fq.val[j]));
        }
      rep.tForce += secondsSince(tic);

      // ---- 4+5. Linearised admissibility metric K_adm ----
      tic = Clock::now();
      std::vector<Eigen::Triplet<Real>> kTriplets;
      for (const auto& ct : cellTables)
        for (const auto& cq : ct.qps)
        {
          const Mat2 F = deformationAtCellQP(u.getData(), ct, cq);
          const Real jK = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
          if (jK <= Real(0))
            continue; // inverted validation point: leave to line search.
          const Real frob2 =
            F(0, 0) * F(0, 0) + F(0, 1) * F(0, 1)
            + F(1, 0) * F(1, 0) + F(1, 1) * F(1, 1);
          const Real qK = frob2 / (Real(2) * jK);

          const Real sJ0 = jK - p.jSafe;
          const Real sQ0 = p.qMax - qK;
          const bool jActive = (sJ0 > Real(0) && sJ0 < p.s0J);
          const bool qActive = (sQ0 > Real(0) && sQ0 < p.s0Q);
          if (!jActive && !qActive)
            continue;

          Mat2 FinvT(2, 2);
          FinvT(0, 0) =  F(1, 1) / jK;
          FinvT(0, 1) = -F(1, 0) / jK;
          FinvT(1, 0) = -F(0, 1) / jK;
          FinvT(1, 1) =  F(0, 0) / jK;

          Mat2 dQdF(2, 2);
          for (int r = 0; r < 2; ++r)
            for (int c = 0; c < 2; ++c)
              dQdF(r, c) =
                (F(r, c) - (frob2 / Real(2)) * FinvT(r, c)) / jK;

          const std::size_t nLocal = ct.dofs.size();
          std::vector<Real> aJ(nLocal), aQ(nLocal);
          for (std::size_t l = 0; l < nLocal; ++l)
          {
            Real accJ = 0, accQ = 0;
            for (int r = 0; r < 2; ++r)
              for (int c = 0; c < 2; ++c)
              {
                accJ += jK * FinvT(r, c) * cq.jac[l](r, c);
                accQ += dQdF(r, c) * cq.jac[l](r, c);
              }
            aJ[l] = accJ;
            aQ[l] = accQ;
          }

          auto addRank1 =
            [&](const std::vector<Real>& a, Real bpp, Real bp,
                Real gamma, Real rhsSign)
          {
            if (bpp != Real(0))
              for (std::size_t i = 0; i < nLocal; ++i)
                for (std::size_t j = 0; j < nLocal; ++j)
                  kTriplets.emplace_back(
                      ct.dofs[i], ct.dofs[j], gamma * bpp * a[i] * a[j]);
            if (p.includeAdmissibilityGradient && bp != Real(0))
              for (std::size_t i = 0; i < nLocal; ++i)
                rhs(ct.dofs[i]) -= rhsSign * gamma * bp * a[i];
          };

          if (jActive)
          {
            const Real bp = -Real(1) / sJ0 + Real(1) / p.s0J;
            const Real bpp = Real(1) / (sJ0 * sJ0);
            addRank1(aJ, bpp, bp, p.gammaJ, Real(1));
          }
          if (qActive)
          {
            const Real bp = -Real(1) / sQ0 + Real(1) / p.s0Q;
            const Real bpp = Real(1) / (sQ0 * sQ0);
            addRank1(aQ, bpp, bp, p.gammaQ, Real(-1));
          }
        }
      rep.tBarrier += secondsSince(tic);

      // ---- Solve: cached M factorization when the dynamic metric is empty ----
      tic = Clock::now();
      for (std::size_t i = 0; i < nDofs; ++i)
        if (isBoundaryDof[i])
          rhs(i) = Real(0);
      Math::Vector<Real> vK;
      if (obsTriplets.empty() && kTriplets.empty() && mLdltOK)
      {
        rep.tFactor += secondsSince(tic);
        tic = Clock::now();
        vK = mLdlt.solve(rhs);
      }
      else
      {
        std::vector<Eigen::Triplet<Real>> all;
        all.reserve(
            mTriplets.size() + obsTriplets.size() + kTriplets.size() + nDofs);
        for (const auto& t : mTriplets)
          if (!isBoundaryDof[t.row()] && !isBoundaryDof[t.col()])
            all.push_back(t);
        for (const auto& t : obsTriplets)
          if (!isBoundaryDof[t.row()] && !isBoundaryDof[t.col()])
            all.push_back(t);
        for (const auto& t : kTriplets)
          if (!isBoundaryDof[t.row()] && !isBoundaryDof[t.col()])
            all.push_back(t);
        for (std::size_t i = 0; i < nDofs; ++i)
          if (isBoundaryDof[i])
            all.emplace_back(i, i, Real(1));
        Eigen::SparseMatrix<Real> sys(nDofs, nDofs);
        sys.setFromTriplets(all.begin(), all.end());
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>> ldlt(sys);
        if (ldlt.info() != Eigen::Success)
        {
          rep.exitReason = "solve-factorization-failed";
          break;
        }
        rep.tFactor += secondsSince(tic);
        tic = Clock::now();
        vK = ldlt.solve(rhs);
      }
      if (!vK.allFinite())
      {
        rep.exitReason = "solve-nonfinite";
        break;
      }
      rep.tSolve += secondsSince(tic);

      // ---- 6. Optimal step rescale β = ⟨d,v⟩_Γ/⟨v,v⟩_Γ ----
      if (p.betaMax > Real(1))
      {
        Real bNum = 0, bDen = 0;
        for (const auto& ft : facetTables)
          for (const auto& fq : ft.qps)
          {
            const Vec2 y = fq.X + fieldAtFacetQP(u.getData(), ft, fq);
            const Real r = phiAt(y);
            const Vec2 g = gradAt(y);
            const Real obsWeight =
                g.dot(g) + epsG
              + (p.residualStabilizedObservationMetric
                  ? (r * r) / sigma2 : Real(0));
            const Real omega = std::exp(-r * r / sigma2);
            const Vec2 d = vec2(
                -omega * r * g(0) / obsWeight,
                -omega * r * g(1) / obsWeight);
            const Vec2 v = fieldAtFacetQP(vK, ft, fq);
            bNum += fq.w * d.dot(v);
            bDen += fq.w * v.dot(v);
          }
        if (bDen > Real(0) && std::isfinite(bNum))
        {
          const Real beta =
            std::clamp(bNum / bDen, Real(1), p.betaMax);
          vK *= beta;
        }
      }

      const Real maxStep = vK.cwiseAbs().maxCoeff();
      if (maxStep <= stepTol)
      {
        rep.exitReason = "step-below-stepTol";
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

      rep.lastAlpha = alpha;
      rep.acceptedStep = alpha * maxStep;
      rep.minJ = adm.minJ;
      rep.maxQRel = adm.maxQ;

      // ---- Diagnostics + stopping ----
      const auto surf = surfaceState(u.getData());
      const Real eNow = p.energyLineSearch && std::isfinite(eTrial)
        ? eTrial : surf.energy;
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
                  << "  actSup=" << surf.activeSup
                  << "  actFrac=" << rep.activeFraction
                  << "  α=" << alpha
                  << "  bt=" << backtracks
                  << "  min_j=" << rep.minJ
                  << "  max_Q=" << rep.maxQRel
                  << '\n';

      if (surf.activeRMS <= activeRMSTol)
      {
        rep.exitReason = "active-rms-converged";
        ++rep.iterations;
        break;
      }
      if (surf.activeSup <= activeSupTol)
      {
        rep.exitReason = "active-sup-converged";
        ++rep.iterations;
        break;
      }
      const Real eRel =
        std::abs(ePrev - eNow) / std::max(ePrev, Real(1e-30));
      if (eRel < p.energyStagTol)
      {
        rep.exitReason = "energy-stagnation";
        ++rep.iterations;
        break;
      }
      ePrev = eNow;
    }
    return rep;
  }
}

#endif
