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
 * @brief Intrinsic shape quality energy + Jacobian admissibility /
 *        nondegeneracy / singularity barrier — fully analytical residual
 *        and tangent for the 2D triangular affine prototype.
 *
 * Scope.
 *   PROTOTYPE-SPECIALISED: 2D, planar mesh, affine triangle cells,
 *   P1 vector dofs (three vertices, two components each — six local dofs
 *   per cell). NOT FES-independent.
 *
 * Mathematics.
 *   Intrinsic admissibility:
 *     j_K^u(xhat) = sigma_K * det(A_K^u(xhat)) / J_K_scale,
 *     admissibility condition: j_K^u > j_min.
 *
 *   Intrinsic shape quality:
 *     Q_shape(A_K^u) = (1/d) ||A_K^u||_F^2 / (sigma_K det(A_K^u))^(2/d),
 *     valid only when sigma_K det(A_K^u) > 0.
 *
 *   Per-cell energy (area-weighted at assembly):
 *     E_cell = w_s * (Q_shape - 1) - w_d * log(j_K^u - j_min)
 *   with w_s = gamma * |K|/|Omega|, w_d = beta * |K|/|Omega|.
 *
 * Form-language API.
 *   gamma and beta — the shape and admissibility energy weights — are
 *   passed as `Variational::RealFunctionBase<...>` objects, mirroring
 *   the (lambda, mu) Lamé parameters in
 *   `Rodin::Variational::LinearElasticityIntegral`. The integrator
 *   evaluates them at the cell centroid (in physical coordinates),
 *   pulling out the per-cell scalar values before calling the
 *   closed-form `evaluateBarrierLocal`.
 *
 *   `jMin` and `domainMeasure` are algorithmic / mesh-global scalars
 *   that don't naturally vary in space and stay as `Real` fields on
 *   `BarrierParameters`.
 *
 *   The user-facing helper `JacobianAdmissibilityBarrier` mirrors
 *   `SignedDistanceRegistration`: variational skeleton + cell cache
 *   captured at construction; gamma, beta and params supplied at
 *   `operator()` / `Tangent` / `Residual`.
 *
 * Naming.
 *   The log(j_K^u - j_min) term is the "Jacobian admissibility /
 *   nondegeneracy / singularity barrier". Never an "orientation barrier".
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
   * @brief Scalar parameters for the shape + admissibility energy.
   *
   *   jMin          = admissibility floor (strictly positive)
   *   domainMeasure = |Omega| (mesh-global normalisation)
   *
   * The spatial weights gamma and beta are NOT carried here — they are
   * `Variational::RealFunctionBase` objects passed alongside.
   */
  struct BarrierParameters
  {
    Real jMin = 1e-3;
    Real domainMeasure = 0;
  };

  struct BarrierLocal
  {
    Real energy = 0;
    Real j = 0;
    bool valid = true;
    Math::Vector<Real> grad;
    Math::Matrix<Real> hess;
  };

  /**
   * @brief Process-wide counter for cells that hit j_K^u <= j_min during
   *        assembly. The Newton driver reads it after each iteration.
   */
  inline std::atomic<unsigned long>& barrierInadmissibleCount()
  {
    static std::atomic<unsigned long> count{0};
    return count;
  }

  /**
   * @brief Closed-form local barrier energy + first and second
   *        derivatives w.r.t. the six nodal coefficients. `gamma` and
   *        `beta` are evaluated scalar values for this cell.
   */
  inline BarrierLocal evaluateBarrierLocal(
      const CellGeomCache& cell,
      const std::array<Math::SpatialVector<Real>, 3>& uv,
      Real gamma, Real beta,
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
    const Real w_d = beta  * areaWeight;

    out.energy =
        w_s * (qShape - Real(1))
      - w_d * std::log(j - params.jMin);

    const Math::SpatialMatrix<Real> FinvT = F.inverse().transpose();
    const Math::SpatialMatrix<Real> FC = F * cell.C;
    const Real h = j / (j - params.jMin);
    const Real Kd = params.jMin * j
                  / ((j - params.jMin) * (j - params.jMin));
    const Real wsOverD = w_s / D;
    Math::SpatialMatrix<Real> S =
        wsOverD * (FC - Real(0.5) * n2 * FinvT)
      - w_d * h * FinvT;

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

            const Real det =
                Kd * FinvTgk[k](l) * FinvTgk[m](p)
              + h  * FinvTgk[k](p) * FinvTgk[m](l);

            out.hess(k * 2 + l, m * 2 + p) =
                wsOverD * shape + w_d * det;
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
   * @brief Linear-form integrator: shape + Jacobian admissibility residual.
   *
   * `gamma` and `beta` are stored as `unique_ptr<RealFunctionBase<...>>`
   * after `.copy()` — same pattern as `LinearElasticityIntegrator`'s
   * `m_lambda` / `m_mu`. They're queried at the cell centroid per cell.
   */
  template <class GammaDerived, class BetaDerived,
            class TestType, class StateType>
  class BarrierResidualIntegrator final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using GammaType = Variational::RealFunctionBase<GammaDerived>;
      using BetaType  = Variational::RealFunctionBase<BetaDerived>;

      BarrierResidualIntegrator(
          const GammaType& gamma, const BetaType& beta,
          const std::vector<CellGeomCache>& cellCache,
          const std::unordered_map<Index, std::size_t>& cellToLocal,
          const TestType& v,
          const StateType& u,
          BarrierParameters params)
        : Variational::LinearFormIntegratorBase<Real>(v),
          m_gamma(gamma.copy()), m_beta(beta.copy()),
          m_cellCache(cellCache), m_cellToLocal(cellToLocal),
          m_test(v), m_state(u), m_params(params)
      {}

      BarrierResidualIntegrator(const BarrierResidualIntegrator& other)
        : Variational::LinearFormIntegratorBase<Real>(other),
          m_gamma(other.m_gamma ? other.m_gamma->copy() : nullptr),
          m_beta(other.m_beta ? other.m_beta->copy() : nullptr),
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

        // Evaluate gamma and beta at the cell centroid.
        const Geometry::Point centroid(
            polytope, detail::triangleReferenceCentroid());
        const Real gammaVal = m_gamma->getValue(centroid);
        const Real betaVal  = m_beta->getValue(centroid);

        const auto local = evaluateBarrierLocal(
            cell, uv, gammaVal, betaVal, m_params);
        if (!local.valid)
        {
          barrierInadmissibleCount().fetch_add(1, std::memory_order_relaxed);
          return *this;
        }
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
      std::unique_ptr<BetaType>  m_beta;
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
   * @brief Bilinear-form integrator: shape + Jacobian admissibility tangent.
   */
  template <class GammaDerived, class BetaDerived,
            class TrialType, class TestType, class StateType>
  class BarrierTangentIntegrator final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using GammaType = Variational::RealFunctionBase<GammaDerived>;
      using BetaType  = Variational::RealFunctionBase<BetaDerived>;

      BarrierTangentIntegrator(
          const GammaType& gamma, const BetaType& beta,
          const std::vector<CellGeomCache>& cellCache,
          const std::unordered_map<Index, std::size_t>& cellToLocal,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          BarrierParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_gamma(gamma.copy()), m_beta(beta.copy()),
          m_cellCache(cellCache), m_cellToLocal(cellToLocal),
          m_trial(du), m_test(v), m_state(u), m_params(params)
      {}

      BarrierTangentIntegrator(const BarrierTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_gamma(other.m_gamma ? other.m_gamma->copy() : nullptr),
          m_beta(other.m_beta ? other.m_beta->copy() : nullptr),
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
        const Real betaVal  = m_beta->getValue(centroid);

        const auto local = evaluateBarrierLocal(
            cell, uv, gammaVal, betaVal, m_params);
        if (!local.valid)
        {
          barrierInadmissibleCount().fetch_add(1, std::memory_order_relaxed);
          for (std::size_t i = 0; i < 6; ++i)
            m_matrix(i, i) = Real(1);
          return *this;
        }
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
      std::unique_ptr<BetaType>  m_beta;
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
   * @brief High-level helper for the shape + Jacobian admissibility
   *        barrier, modelled on `SignedDistanceRegistration` and
   *        `LinearElasticityIntegral`.
   *
   * Holds the variational skeleton `(du, v, u)` together with the
   * per-cell geometry cache and the mesh-cell-index -> local-index map
   * — both of which are tied to the mesh of `(du, v, u)` and don't
   * change between calls. `gamma` and `beta` (energy weights, possibly
   * spatially varying) and `BarrierParameters` (the scalar jMin /
   * domainMeasure) arrive at `operator()` / `Tangent` / `Residual`.
   *
   * Usage:
   *
   *     JacobianAdmissibilityBarrier barrier(du, v, u, cellCache, cellToLocal);
   *
   *     newton = barrier(gamma, beta, params);
   *
   *     // Mixing form-language gammaFn with a constant beta:
   *     newton = barrier(gammaFn, 0.01, params);
   *
   *     // Decomposed access:
   *     newton = barrier.Tangent(gamma, beta, params)
   *            + barrier.Residual(gamma, beta, params);
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
      template <class GammaDerived, class BetaDerived>
      auto Tangent(
          const Variational::RealFunctionBase<GammaDerived>& gamma,
          const Variational::RealFunctionBase<BetaDerived>& beta,
          BarrierParameters params) const
      {
        return BarrierTangentIntegrator<
                  GammaDerived, BetaDerived,
                  Trial, Test, State>(
              gamma, beta,
              m_cellCache.get(), m_cellToLocal.get(),
              m_du.get(), m_v.get(), m_u.get(),
              params);
      }

      // ---- Decomposed: Residual ---------------------------------------------
      template <class GammaDerived, class BetaDerived>
      auto Residual(
          const Variational::RealFunctionBase<GammaDerived>& gamma,
          const Variational::RealFunctionBase<BetaDerived>& beta,
          BarrierParameters params) const
      {
        return BarrierResidualIntegrator<
                  GammaDerived, BetaDerived,
                  Test, State>(
              gamma, beta,
              m_cellCache.get(), m_cellToLocal.get(),
              m_v.get(), m_u.get(),
              params);
      }

      // ---- Composite --------------------------------------------------------
      //
      // For constant weights, wrap a `Real` value via
      // `Variational::RealFunction<Real>(value)` at the call site —
      // exactly the same idiom Rodin's other form-language helpers
      // expect for scalar constants.
      template <class GammaDerived, class BetaDerived>
      auto operator()(
          const Variational::RealFunctionBase<GammaDerived>& gamma,
          const Variational::RealFunctionBase<BetaDerived>& beta,
          BarrierParameters params) const
      {
        return Tangent(gamma, beta, params)
             + Residual(gamma, beta, params);
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
