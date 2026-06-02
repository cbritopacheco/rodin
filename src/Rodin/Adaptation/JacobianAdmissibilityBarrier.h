/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_JACOBIANADMISSIBILITYBARRIER_H
#define RODIN_ADAPTATION_JACOBIANADMISSIBILITYBARRIER_H

/**
 * @file
 * @brief Intrinsic shape quality energy `E_shape` — fully analytical
 *        residual and tangent for the 2D triangular affine prototype.
 *
 * Scope.
 *   PROTOTYPE-SPECIALISED: 2D, planar mesh, affine triangle cells,
 *   P1 vector dofs (three vertices, two components each — six local
 *   dofs per cell). NOT FES-independent.
 *
 * Mathematics.
 *   Branch-aware admissibility scalar:
 *     j_K^u(xhat) = sigma_K * det(A_K^u(xhat)) / J_K_scale,
 *     j_K^u > j_min  (only condition; sigma_K carries the initial
 *                     orientation, j_K^u = 1 has no special meaning).
 *
 *   Intrinsic shape quality (independent of orientation):
 *     Q_shape(A_K^u) = (1/d) ||A_K^u||_F^2 / (sigma_K det(A_K^u))^(2/d),
 *     well-defined whenever sigma_K det(A_K^u) > 0.
 *
 *   Per-cell energy (area-weighted at assembly):
 *     E_cell = w_s * (Q_shape - 1),
 *   with w_s = gamma * |K|/|Omega|. The global energy is `E_shape(u)`.
 *
 *   The line search of `LSRAdmissibility` is the sole admissibility
 *   mechanism: with the safety margin xi > 1, every accepted iterate
 *   stays strictly above j_safe and the older floor-barrier B(j) that
 *   used to live here is identically zero. It has been removed; the
 *   singular-cell fallback below remains as a defensive measure for
 *   the (rare) case where a cell becomes formally singular during
 *   assembly outside of a line-search-protected step (e.g. inside the
 *   Hilbert lift linear solve).
 *
 * Form-language API.
 *   gamma (shape weight) is passed as a
 *   `Variational::RealFunctionBase<...>` object, mirroring `lambda`
 *   /`mu` in `Rodin::Variational::LinearElasticityIntegral`. The
 *   integrator evaluates it at the cell centroid (in physical
 *   coordinates), pulling out the per-cell scalar value before calling
 *   the closed-form `evaluateBarrierLocal`.
 *
 *   `jMin` (singular floor ratio) and `domainMeasure` are algorithmic
 *   scalars carried on `BarrierParameters`. `jMin` is dimensionless on
 *   j_K^u = sigma_K det(A_K^u) / J_K_scale (recommended ~ 1e-8).
 *
 *   The user-facing helper `JacobianAdmissibilityBarrier` mirrors
 *   `LSRRegistration`: variational skeleton + cell cache captured at
 *   construction; gamma and params supplied at `operator()` /
 *   `Tangent` / `Residual`. The class name is retained for backward
 *   compatibility; conceptually the energy it produces is `E_shape`.
 */

#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/RealFunction.h"

#include "CellGeomCache.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Scalar parameters for the shape-quality energy.
   *
   *   jMin          = dimensionless singularity floor ratio for the
   *                   defensive singular-cell fallback (strictly
   *                   positive, e.g. 1e-8).
   *   domainMeasure = |Omega| (mesh-global normalisation)
   *
   * The spatial weight gamma (shape) is NOT carried here — it is a
   * `Variational::RealFunctionBase` object passed alongside.
   */
  struct BarrierParameters
  {
    Real jMin = 1e-8;  ///< dimensionless ratio, not raw det(A_K^u)
    Real domainMeasure = 0;

    /// Smooth Q-shape barrier (interior-point penalty on cell
    /// anisotropy). Inactive by default: when `qBarrierWeight = 0` or
    /// `qBarrierMax = +inf`, the barrier contributes nothing.
    ///
    /// Active form on (qBarrierAct, qBarrierMax):
    ///   B_Q(Q) = -log((Qmax - Q) / (Qmax - Qact)) - (Q - Qact)/(Qmax - Qact)
    /// with B_Q(Qact) = B_Q'(Qact) = 0 and B_Q(Q) -> +oo as Q -> Qmax.
    /// Identically zero for Q <= Qact.
    ///
    /// Per-cell density: w_Q = qBarrierWeight * |K| / domainMeasure,
    /// energy contribution w_Q * B_Q(Q_K^u).
    Real qBarrierWeight    = 0;
    Real qBarrierAct       = std::numeric_limits<Real>::infinity();
    Real qBarrierMax       = std::numeric_limits<Real>::infinity();
  };

  struct BarrierLocal
  {
    Real energy = 0;
    Real j = 0;
    Real qShape = 0;            ///< intrinsic shape quality (>= 1)
    bool valid = true;          ///< j > jMin (cell is on its initial branch)
    bool qBarrierActive = false; ///< Q > Qact (B_Q contribution nonzero)
    Math::Vector<Real> grad;
    Math::Matrix<Real> hess;
  };

  /**
   * @brief Process-wide counter for cells that hit j_K^u <= jMin during
   *        assembly (singularity hit). The Newton driver reads it after
   *        each iteration.
   */
  inline std::atomic<unsigned long>& barrierInadmissibleCount()
  {
    static std::atomic<unsigned long> count{0};
    return count;
  }

  /**
   * @brief Diagonal magnitude used for the local block when a cell hits
   *        j_K^u <= jMin during assembly.
   *
   * The fallback freezes the cell's six DOFs for this Newton step by
   * contributing a large multiple of the identity to the local Hessian
   * (and zero to the residual). The penalty must dominate any other
   * contribution to that cell — most notably the LSR rank-one block,
   * which scales like `rho_s * params.normalizer * area` and can easily
   * reach `O(10)` per cell on realistic mesh / LSR weights.
   *
   * `1e10` is large enough that the resulting Newton update in the
   * frozen DOFs is below floating-point noise:
   *     du ~ R / 1e10 ~ 1e-10 * R,
   * which is at most `1e-8` for any realistic `R`. The remainder of the
   * mesh continues to take its proper Newton step. The cell is allowed
   * back into the iteration as soon as its `j_K^u` re-enters the
   * admissible region, which happens automatically because the Newton
   * step at neighbouring cells continues unobstructed.
   *
   * The previous fallback used a unit diagonal (`1`), which is the same
   * magnitude as the LSR contribution on typical meshes and therefore
   * failed to freeze the cell — producing massive Newton updates in the
   * singular DOFs (we observed `step ~ 14` on a unit-square mesh in the
   * wavy-circle sweep) that immediately drove tens of neighbouring
   * cells through the singular floor and triggered a cascade.
   */
  inline constexpr Real kBarrierSingularPenalty = Real(1e10);

  /**
   * @brief Process-wide minimum of j_K^u over all assembled cells
   *        (running min across threads; reset externally if needed).
   *        Stored as the bit-cast of a `double` in an atomic, updated
   *        via compare-exchange.
   */
  inline std::atomic<Real>& barrierMinJ()
  {
    static std::atomic<Real> value{std::numeric_limits<Real>::infinity()};
    return value;
  }

  inline void barrierUpdateMinJ(Real candidate)
  {
    Real cur = barrierMinJ().load(std::memory_order_relaxed);
    while (candidate < cur
        && !barrierMinJ().compare_exchange_weak(
              cur, candidate,
              std::memory_order_relaxed,
              std::memory_order_relaxed))
    {
      // cur updated by compare_exchange_weak; loop.
    }
  }

  /**
   * @brief Closed-form local shape-quality energy + first and second
   *        derivatives w.r.t. the six nodal coefficients. `gamma` is
   *        the evaluated scalar weight for this cell.
   */
  inline BarrierLocal evaluateBarrierLocal(
      const CellGeomCache& cell,
      const std::array<Math::SpatialVector<Real>, 3>& uv,
      Real gamma,
      const BarrierParameters& params)
  {
    BarrierLocal out;
    out.grad = Math::Vector<Real>::Zero(6);
    out.hess = Math::Matrix<Real>::Zero(6, 6);

    Math::SpatialMatrix<Real> U(2, 3);
    for (std::size_t k = 0; k < 3; ++k)
    {
      U(0, k) = uv[k](0);
      U(1, k) = uv[k](1);
    }
    Math::SpatialMatrix<Real> F =
      Math::SpatialMatrix<Real>::Identity(2, 2) + U * cell.gradN;

    const Real detF = F.determinant();
    const Real sigDetA = static_cast<Real>(cell.sigmaK) * cell.detAK;
    const Real j = sigDetA * detF / cell.Jscale;
    out.j = j;
    if (j <= params.jMin)
    {
      out.valid = false;
      return out;
    }

    const Math::SpatialMatrix<Real> M = F * cell.A;
    const Real n2 = M.squaredNorm();
    const Real sigDetM = static_cast<Real>(cell.sigmaK) * M.determinant();
    constexpr Real d = 2;
    const Real D = std::pow(sigDetM, Real(2) / d);
    const Real qShape = n2 / (d * D);

    const Real areaWeight = cell.area / params.domainMeasure;
    const Real w_s = gamma * areaWeight;
    const Real w_Q = params.qBarrierWeight * areaWeight;

    out.qShape = qShape;

    // Q-shape barrier B_Q (interior-point penalty on cell anisotropy).
    // Active on (Qact, Qmax); identically zero with vanishing first and
    // second derivatives at Q <= Qact. The derivation enforces
    //   B_Q(Qact) = 0,  B_Q'(Qact) = 0,  B_Q(Q) -> +oo as Q -> Qmax-.
    //
    // For Qact < Q < Qmax:
    //   Bp = 1/(Qmax-Q) - 1/(Qmax-Qact)  (>= 0)
    //   Bpp = 1/(Qmax-Q)^2               (> 0)
    Real B_Q   = 0;
    Real B_Qp  = 0;
    Real B_Qpp = 0;
    const bool qBarrierEligible =
         params.qBarrierWeight > Real(0)
      && std::isfinite(params.qBarrierMax)
      && std::isfinite(params.qBarrierAct)
      && params.qBarrierAct < params.qBarrierMax;
    if (qBarrierEligible && qShape > params.qBarrierAct)
    {
      if (qShape >= params.qBarrierMax)
      {
        // Outside the barrier domain; report infeasibility.
        out.valid = false;
        out.qBarrierActive = true;
        return out;
      }
      out.qBarrierActive = true;
      const Real span = params.qBarrierMax - params.qBarrierAct;
      const Real toMax = params.qBarrierMax - qShape;
      const Real invToMax = Real(1) / toMax;
      const Real invSpan  = Real(1) / span;
      B_Q   = -std::log(toMax * invSpan) - (qShape - params.qBarrierAct) * invSpan;
      B_Qp  = invToMax - invSpan;
      B_Qpp = invToMax * invToMax;
    }

    out.energy = w_s * (qShape - Real(1)) + w_Q * B_Q;

    // Combined effective shape weight for gradient + scaled Hessian.
    // Density: w_s*(Q-1) + w_Q*B_Q(Q); gradient: w_eff * dQ/dF; Hessian:
    //   w_eff * d2Q/dF^2  +  w_Q * B_Qpp * dQ/dF tensor dQ/dF
    // where w_eff = w_s + w_Q * B_Qp.
    const Real w_eff = w_s + w_Q * B_Qp;

    const Math::SpatialMatrix<Real> FinvT = F.inverse().transpose();
    const Math::SpatialMatrix<Real> FC = F * cell.C;
    const Real weffOverD = w_eff / D;
    Math::SpatialMatrix<Real> S =
        weffOverD * (FC - Real(0.5) * n2 * FinvT);

    for (std::size_t k = 0; k < 3; ++k)
      for (std::size_t l = 0; l < 2; ++l)
      {
        Real sum = 0;
        for (std::size_t jj = 0; jj < 2; ++jj)
          sum += S(l, jj) * cell.gradN(k, jj);
        out.grad(k * 2 + l) = sum;
      }

    std::array<Math::SpatialVector<Real>, 3> gN;
    std::array<Math::SpatialVector<Real>, 3> FCgk;
    std::array<Math::SpatialVector<Real>, 3> FinvTgk;
    for (std::size_t k = 0; k < 3; ++k)
    {
      gN[k] = Math::SpatialVector<Real>(2);
      gN[k](0) = cell.gradN(k, 0);
      gN[k](1) = cell.gradN(k, 1);
      FCgk[k] = FC * gN[k];
      FinvTgk[k] = FinvT * gN[k];
    }

    const Real halfN = Real(0.5) * n2;

    // Precompute dQ/dF . grad N_k per local dof (k*2+l): dQ/dF = (1/D) * (FC - (N/2)*FinvT)
    // dQdF_dot_g[k][l] = sum_j (FC - (N/2)FinvT)(l, j) * gN[k](j) / D
    //                 = (FCgk[k](l) - (N/2) * FinvTgk[k](l)) / D.
    std::array<std::array<Real, 2>, 3> dQdF_dot_g;
    if (out.qBarrierActive)
    {
      const Real invD = Real(1) / D;
      for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t l = 0; l < 2; ++l)
          dQdF_dot_g[k][l] =
              invD * (FCgk[k](l) - halfN * FinvTgk[k](l));
    }

    for (std::size_t k = 0; k < 3; ++k)
      for (std::size_t m = 0; m < 3; ++m)
      {
        const Real gCg = gN[k].dot(cell.C * gN[m]);
        for (std::size_t l = 0; l < 2; ++l)
          for (std::size_t p = 0; p < 2; ++p)
          {
            const Real shape =
                (l == p ? gCg : Real(0))
              - FCgk[k](l) * FinvTgk[m](p)
              - FCgk[m](p) * FinvTgk[k](l)
              + halfN * (FinvTgk[k](l) * FinvTgk[m](p)
                       + FinvTgk[k](p) * FinvTgk[m](l));

            Real H = weffOverD * shape;
            if (out.qBarrierActive)
              H += w_Q * B_Qpp * dQdF_dot_g[k][l] * dQdF_dot_g[m][p];
            out.hess(k * 2 + l, m * 2 + p) = H;
          }
      }

    return out;
  }

  /// Read the 3-vertex displacement vectors of a cell from a vector P1
  /// GridFunction using the FES dof-mapping API.
  template <class StateType>
  std::array<Math::SpatialVector<Real>, 3> extractCellU(
      const CellGeomCache& cell, const StateType& u)
  {
    const auto& fes = u.getFiniteElementSpace();
    const auto& uData = u.getData();
    std::array<Math::SpatialVector<Real>, 3> uv;
    for (std::size_t i = 0; i < 3; ++i)
    {
      const Index v = cell.vertices[i];
      const auto& dofs = fes.getDOFs(0, v);
      uv[i] = Math::SpatialVector<Real>(2);
      uv[i](0) = uData(dofs[0]);
      uv[i](1) = uData(dofs[1]);
    }
    return uv;
  }

  namespace detail
  {
    /// Reference centroid of the affine triangle; for this prototype
    /// the barrier energy is per-cell, evaluated at the centroid.
    inline Math::SpatialPoint triangleReferenceCentroid()
    {
      Math::SpatialPoint rc(2);
      rc(0) = Real(1) / Real(3);
      rc(1) = Real(1) / Real(3);
      return rc;
    }
  }

  /**
   * @brief Linear-form integrator: shape-quality residual.
   *
   * `gamma` is stored as `unique_ptr<RealFunctionBase<...>>` after
   * `.copy()` — same pattern as `LinearElasticityIntegrator`'s
   * `m_lambda` / `m_mu`. It is queried at the cell centroid per cell.
   */
  template <class GammaDerived,
            class TestType, class StateType>
  class BarrierResidualIntegrator final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using GammaType = Variational::RealFunctionBase<GammaDerived>;

      BarrierResidualIntegrator(
          const GammaType& gamma,
          const std::vector<CellGeomCache>& cellCache,
          const std::unordered_map<Index, std::size_t>& cellToLocal,
          const TestType& v,
          const StateType& u,
          BarrierParameters params)
        : Variational::LinearFormIntegratorBase<Real>(v),
          m_gamma(gamma.copy()),
          m_cellCache(cellCache), m_cellToLocal(cellToLocal),
          m_test(v), m_state(u), m_params(params)
      {}

      BarrierResidualIntegrator(const BarrierResidualIntegrator& other)
        : Variational::LinearFormIntegratorBase<Real>(other),
          m_gamma(other.m_gamma ? other.m_gamma->copy() : nullptr),
          m_cellCache(other.m_cellCache), m_cellToLocal(other.m_cellToLocal),
          m_test(other.m_test), m_state(other.m_state),
          m_params(other.m_params),
          m_polytope(other.m_polytope),
          m_elemVec()
      {}

      BarrierResidualIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const Index cellIdx = polytope.getIndex();
        m_elemVec.resize(6);
        m_elemVec.setZero();
        const auto it = m_cellToLocal.get().find(cellIdx);
        if (it == m_cellToLocal.get().end())
          throw std::runtime_error(
              "BarrierResidualIntegrator: cell "
              + std::to_string(cellIdx) + " not in cache.");
        const auto& cell = m_cellCache.get()[it->second];
        const auto uv = extractCellU(cell, m_state.get());

        // Evaluate gamma at the cell centroid.
        const Geometry::Point centroid(
            polytope, detail::triangleReferenceCentroid());
        const Real gammaVal = m_gamma->getValue(centroid);

        const auto local = evaluateBarrierLocal(
            cell, uv, gammaVal, m_params);
        if (!local.valid)
        {
          barrierInadmissibleCount().fetch_add(1, std::memory_order_relaxed);
          return *this;
        }
        barrierUpdateMinJ(local.j);
        m_elemVec = local.grad;
        return *this;
      }

      ScalarType integrate(std::size_t te) final override
      { return m_elemVec(te); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      BarrierResidualIntegrator* copy() const noexcept final override
      { return new BarrierResidualIntegrator(*this); }

    private:
      std::unique_ptr<GammaType> m_gamma;
      std::reference_wrapper<const std::vector<CellGeomCache>> m_cellCache;
      std::reference_wrapper<const std::unordered_map<Index, std::size_t>>
        m_cellToLocal;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_state;
      BarrierParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  /**
   * @brief Bilinear-form integrator: shape-quality tangent.
   */
  template <class GammaDerived,
            class TrialType, class TestType, class StateType>
  class BarrierTangentIntegrator final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using GammaType = Variational::RealFunctionBase<GammaDerived>;

      BarrierTangentIntegrator(
          const GammaType& gamma,
          const std::vector<CellGeomCache>& cellCache,
          const std::unordered_map<Index, std::size_t>& cellToLocal,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          BarrierParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_gamma(gamma.copy()),
          m_cellCache(cellCache), m_cellToLocal(cellToLocal),
          m_trial(du), m_test(v), m_state(u), m_params(params)
      {}

      BarrierTangentIntegrator(const BarrierTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_gamma(other.m_gamma ? other.m_gamma->copy() : nullptr),
          m_cellCache(other.m_cellCache), m_cellToLocal(other.m_cellToLocal),
          m_trial(other.m_trial), m_test(other.m_test), m_state(other.m_state),
          m_params(other.m_params),
          m_polytope(other.m_polytope),
          m_matrix()
      {}

      BarrierTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const Index cellIdx = polytope.getIndex();
        m_matrix.resize(6, 6);
        m_matrix.setZero();
        const auto it = m_cellToLocal.get().find(cellIdx);
        if (it == m_cellToLocal.get().end())
          throw std::runtime_error(
              "BarrierTangentIntegrator: cell "
              + std::to_string(cellIdx) + " not in cache.");
        const auto& cell = m_cellCache.get()[it->second];
        const auto uv = extractCellU(cell, m_state.get());

        const Geometry::Point centroid(
            polytope, detail::triangleReferenceCentroid());
        const Real gammaVal = m_gamma->getValue(centroid);

        const auto local = evaluateBarrierLocal(
            cell, uv, gammaVal, m_params);
        if (!local.valid)
        {
          barrierInadmissibleCount().fetch_add(1, std::memory_order_relaxed);
          // Freeze this cell's six DOFs for one Newton step by a large
          // identity penalty. See `kBarrierSingularPenalty` for the
          // rationale and the value chosen.
          for (std::size_t i = 0; i < 6; ++i)
            m_matrix(i, i) = kBarrierSingularPenalty;
          return *this;
        }
        barrierUpdateMinJ(local.j);
        m_matrix = local.hess;
        return *this;
      }

      ScalarType integrate(std::size_t tr, std::size_t te) final override
      { return m_matrix(te, tr); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      BarrierTangentIntegrator* copy() const noexcept final override
      { return new BarrierTangentIntegrator(*this); }

    private:
      std::unique_ptr<GammaType> m_gamma;
      std::reference_wrapper<const std::vector<CellGeomCache>> m_cellCache;
      std::reference_wrapper<const std::unordered_map<Index, std::size_t>>
        m_cellToLocal;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_state;
      BarrierParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  /**
   * @brief High-level helper for the shape-quality energy, modelled on
   *        `LSRRegistration` and `LinearElasticityIntegral`.
   *
   * Holds the variational skeleton `(du, v, u)` together with the
   * per-cell geometry cache and the mesh-cell-index -> local-index map
   * — both of which are tied to the mesh of `(du, v, u)` and don't
   * change between calls. `gamma` (energy weight, possibly spatially
   * varying) and `BarrierParameters` (the scalar jMin / domainMeasure)
   * arrive at `operator()` / `Tangent` / `Residual`.
   *
   * Usage:
   *
   *     JacobianAdmissibilityBarrier barrier(du, v, u, cellCache, cellToLocal);
   *
   *     newton = barrier(gamma, params);
   *
   *     // Decomposed access:
   *     newton = barrier.Tangent(gamma, params)
   *            + barrier.Residual(gamma, params);
   */
  template <class Trial, class Test, class State>
  class JacobianAdmissibilityBarrier
  {
    public:
      JacobianAdmissibilityBarrier(
          const Trial& du, const Test& v, const State& u,
          const std::vector<CellGeomCache>& cellCache,
          const std::unordered_map<Index, std::size_t>& cellToLocal)
        : m_du(du), m_v(v), m_u(u),
          m_cellCache(cellCache), m_cellToLocal(cellToLocal)
      {}

      // ---- Decomposed: Tangent ----------------------------------------------
      template <class GammaDerived>
      auto Tangent(
          const Variational::RealFunctionBase<GammaDerived>& gamma,
          BarrierParameters params) const
      {
        return BarrierTangentIntegrator<
                  GammaDerived,
                  Trial, Test, State>(
              gamma,
              m_cellCache.get(), m_cellToLocal.get(),
              m_du.get(), m_v.get(), m_u.get(),
              params);
      }

      // ---- Decomposed: Residual ---------------------------------------------
      template <class GammaDerived>
      auto Residual(
          const Variational::RealFunctionBase<GammaDerived>& gamma,
          BarrierParameters params) const
      {
        return BarrierResidualIntegrator<
                  GammaDerived,
                  Test, State>(
              gamma,
              m_cellCache.get(), m_cellToLocal.get(),
              m_v.get(), m_u.get(),
              params);
      }

      // ---- Composite --------------------------------------------------------
      template <class GammaDerived>
      auto operator()(
          const Variational::RealFunctionBase<GammaDerived>& gamma,
          BarrierParameters params) const
      {
        return Tangent(gamma, params) + Residual(gamma, params);
      }

    private:
      std::reference_wrapper<const Trial> m_du;
      std::reference_wrapper<const Test>  m_v;
      std::reference_wrapper<const State> m_u;
      std::reference_wrapper<const std::vector<CellGeomCache>> m_cellCache;
      std::reference_wrapper<const std::unordered_map<Index, std::size_t>>
        m_cellToLocal;
  };

  template <class Trial, class Test, class State>
  JacobianAdmissibilityBarrier(
      const Trial&, const Test&, const State&,
      const std::vector<CellGeomCache>&,
      const std::unordered_map<Index, std::size_t>&)
    -> JacobianAdmissibilityBarrier<Trial, Test, State>;
}

#endif
