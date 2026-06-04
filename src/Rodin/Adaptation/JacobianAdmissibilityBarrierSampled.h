/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_JACOBIANADMISSIBILITYBARRIERSAMPLED_H
#define RODIN_ADAPTATION_JACOBIANADMISSIBILITYBARRIERSAMPLED_H

/**
 * @file
 * @brief Pk/curved-cell sampled relative-distortion energy `E_shape` —
 *        residual
 *        and tangent assembled by quadrature, FES-agnostic.
 *
 * Scope.
 *   Replaces the affine P1 closed-form prototype in
 *   `JacobianAdmissibilityBarrier.h`. Evaluates the relative deformation
 *   gradient F = I + grad u_h(x_q), Q_rel(F_q), and the optional smooth
 *   relative-Q barrier B_Q(Q_rel) at every quadrature point of every cell.
 *   Residual and tangent contributions are accumulated via the
 *   physical-space gradients of the trial and test basis at the quadrature
 *   points. Absolute physical element quality Q_abs(F_q A_q) is deliberately
 *   not part of this energy; it is reported by the admissibility diagnostic.
 *
 *   On affine triangles (constant gradN, constant F) the
 *   per-quadrature-point integrand is itself constant, so any positive
 *   quadrature order integrates the expression exactly.
 *
 * Mathematics.
 *   Per-quadrature-point density:
 *     e(F; gamma, w_Q, w_j) =
 *          (gamma/|Omega|) (Q_rel(F) - 1)
 *        + (w_Q/|Omega|)   B_Q(Q_rel(F))
 *        + (w_j/|Omega|)   B_j(det F)
 *        + (w_v/|Omega|)   0.5 (log det F)^2,
 *     N = ||F||_F^2,  D = det(F)^(2/d),  Q_rel = N / (d D).
 *
 *   Cell energy:
 *     E_K(u_h) = int_K e(F,A; gamma, w_Q) dx
 *              = sum_q w_q |det J_q| * e(F_q, A_q; gamma_q, w_Q).
 *
 *   First derivative w.r.t. F (per qpt):
 *     dQ/dF = (2/(d D)) [F - (N/d) F^{-T}].
 *   Residual contribution to test dof i = (basis b, component alpha):
 *     R_i = sum_q w_q |det J_q| * w_eff(F_q) * <dQ/dF, jac_i_q>,
 *     w_eff = (gamma/|Omega|) + (w_Q/|Omega|) B_Q'(Q),
 *     jac_i_q = e_alpha (grad phi_b(x_q))^T.
 *
 *   Tangent contribution (next turn):
 *     K_ij = sum_q w_q |det J_q| * [
 *              w_eff * d2Q/dF^2 . jac_i_q . jac_j_q
 *            + (w_Q/|Omega|) B_Q''(Q) * <dQ/dF, jac_i_q> <dQ/dF, jac_j_q>
 *           ],
 *     d2Q/dF^2 . (H1, H2) =
 *         (2/(d D)) [<H1, H2> + (N/d) <H1, F^{-T} H2^T F^{-T}>]
 *       - (4/(d^2 D)) [<F, H1><F^{-T}, H2> + <F^{-T}, H1><F, H2>]
 *       + (4 N/(d^3 D)) <F^{-T}, H1><F^{-T}, H2>,
 *
 * Singular cells.
 *   When `det(F_q) <= jMin` at some qpt of cell K, the qpt is
 *   reported via `barrierInadmissibleCount()` and dropped from the
 *   cell-local energy/residual. (The tangent will substitute a frozen
 *   diagonal block exactly as the closed-form does, see
 *   `kBarrierSingularPenalty`.)
 */

#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include <Eigen/Eigenvalues>

#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/RealFunction.h"

#include "JacobianAdmissibilityBarrier.h"  // BarrierParameters, kBarrierSingularPenalty, barrierInadmissibleCount, barrierUpdateMinJ
#include "LSRIntegrators.h"  // lsrQuadOrderFor

namespace Rodin::Adaptation
{
  inline Math::Matrix<Real> projectSymmetricMatrixPSD(
      const Math::Matrix<Real>& matrix)
  {
    Math::Matrix<Real> symmetric =
      Real(0.5) * (matrix + matrix.transpose());
    Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eig(symmetric);
    if (eig.info() != Eigen::Success)
      return symmetric;
    auto values = eig.eigenvalues();
    for (Eigen::Index i = 0; i < values.size(); ++i)
      if (values(i) < Real(0))
        values(i) = Real(0);
    return eig.eigenvectors()
      * values.asDiagonal()
      * eig.eigenvectors().transpose();
  }

  /**
   * @brief Per-quadrature-point sampled-barrier evaluation.
   *
   * Returns the scalar fields needed to assemble both the energy and the
   * cell-local residual contraction. `valid = false` indicates the qpt
   * is dropped (j_K^u <= jMin); the caller decides what to do (skip,
   * freeze tangent block, etc.). All derivative quantities are computed
   * w.r.t. the deformation gradient F.
   */
  struct BarrierSampledQpt
  {
    Math::SpatialMatrix<Real> F;          ///< I + grad u_h
    Math::SpatialMatrix<Real> A;          ///< D x = cell jacobian at qpt
    Math::SpatialMatrix<Real> M;          ///< relative map used by Q_rel (= F)
    Math::SpatialMatrix<Real> FinvT;      ///< F^{-T}
    Math::SpatialMatrix<Real> MAt;        ///< F, retained for shared formulas
    Math::SpatialMatrix<Real> C;          ///< Identity in the relative formula
    Real detA = 0;
    Real sigma = 1;
    Real detM = 0;
    Real sigDetM = 0;                     ///< det F in the relative convention
    Real N = 0;                           ///< ||M||_F^2
    Real D = 0;                           ///< (sigma det M)^(2/d)
    Real Q = 0;                           ///< Q_rel = N / (d D)
    Real j = 0;                           ///< det F, the dimensionless
                                          ///< line-search ratio
    bool valid = true;
  };

  /// Build the per-qpt sampled-barrier evaluation. The state-gradient
  /// matrix `gradU` is the value of `Variational::Jacobian(u_state)` at
  /// the qpt, the cell jacobian `A` is `cell.getTransformation().jacobian(rc)`.
  inline BarrierSampledQpt evaluateBarrierSampledQpt(
      const Math::SpatialMatrix<Real>& gradU,
      const Math::SpatialMatrix<Real>& A,
      Real jMin)
  {
    BarrierSampledQpt out;
    const std::size_t d = static_cast<std::size_t>(A.rows());
    out.A = A;
    out.F = Math::SpatialMatrix<Real>::Identity(
        static_cast<Eigen::Index>(d),
        static_cast<Eigen::Index>(d))
      + gradU;
    out.detA = A.determinant();
    out.sigma = Real(1);
    const Real detF = out.F.determinant();
    out.j = detF;
    if (out.j <= jMin)
    {
      out.valid = false;
      return out;
    }
    out.M = out.F;
    out.detM = out.M.determinant();
    out.sigDetM = out.detM;
    if (out.sigDetM <= Real(0))
    {
      out.valid = false;
      return out;
    }
    out.N = out.M.squaredNorm();
    out.D = std::pow(out.sigDetM, Real(2) / static_cast<Real>(d));
    out.Q = out.N / (static_cast<Real>(d) * out.D);
    out.FinvT = out.F.inverse().transpose();
    out.MAt = out.M;
    out.C = Math::SpatialMatrix<Real>::Identity(
        static_cast<Eigen::Index>(d),
        static_cast<Eigen::Index>(d));
    return out;
  }

  /// Smooth Q-barrier B_Q(Q), B_Q'(Q), B_Q''(Q). Identically zero for
  /// Q <= Qact. Outside [Qact, Qmax) returns infeasible (active=true,
  /// but the caller should treat the qpt as invalid).
  struct BarrierSampledQBarrier
  {
    Real B = 0;
    Real Bp = 0;
    Real Bpp = 0;
    bool active = false;
    bool feasible = true;
  };

  inline BarrierSampledQBarrier evaluateQBarrier(
      Real Q, const BarrierParameters& params)
  {
    BarrierSampledQBarrier out;
    const bool eligible =
         params.qBarrierWeight > Real(0)
      && std::isfinite(params.qBarrierMax)
      && std::isfinite(params.qBarrierAct)
      && params.qBarrierAct < params.qBarrierMax;
    if (!eligible || Q <= params.qBarrierAct)
      return out;
    if (Q >= params.qBarrierMax)
    {
      out.active = true;
      out.feasible = false;
      return out;
    }
    out.active = true;
    const Real span = params.qBarrierMax - params.qBarrierAct;
    const Real toMax = params.qBarrierMax - Q;
    const Real invToMax = Real(1) / toMax;
    const Real invSpan  = Real(1) / span;
    out.B   = -std::log(toMax * invSpan) - (Q - params.qBarrierAct) * invSpan;
    out.Bp  = invToMax - invSpan;
    out.Bpp = invToMax * invToMax;
    return out;
  }

  /// Smooth admissibility barrier B_j(j), B_j'(j), B_j''(j).
  /// Active on (jMin, jSafe); identically zero for j >= jSafe.
  /// With t = (j - jMin)/(jSafe - jMin),
  ///   B_j(t) = -log(t) + t - 1,
  /// so B_j(1) = B_j'(1) = 0 and B_j -> +oo as j -> jMin+.
  struct BarrierSampledJBarrier
  {
    Real B = 0;
    Real Bp = 0;
    Real Bpp = 0;
    bool active = false;
    bool feasible = true;
  };

  inline BarrierSampledJBarrier evaluateJBarrier(
      Real j, const BarrierParameters& params)
  {
    BarrierSampledJBarrier out;
    const bool eligible =
         params.jBarrierWeight > Real(0)
      && params.jBarrierSafeRatio > params.jMin;
    if (!eligible || j >= params.jBarrierSafeRatio)
      return out;
    if (j <= params.jMin)
    {
      out.active = true;
      out.feasible = false;
      return out;
    }
    out.active = true;
    const Real span = params.jBarrierSafeRatio - params.jMin;
    const Real t = (j - params.jMin) / span;
    out.B = -std::log(t) + t - Real(1);
    out.Bp = (Real(1) - Real(1) / t) / span;
    out.Bpp = Real(1) / (t * t * span * span);
    return out;
  }

  /// Centered volume tether B_v(j) = 0.5 (log j)^2 and derivatives
  /// with respect to j. Unlike the lower j-barrier, this has
  /// B_v'(1) = 0 and B_v''(1) = 1, so the identity has zero residual
  /// but an active volumetric tangent.
  struct BarrierSampledVolumeTether
  {
    Real B = 0;
    Real Bp = 0;
    Real Bpp = 0;
    bool active = false;
    bool feasible = true;
  };

  inline BarrierSampledVolumeTether evaluateVolumeTether(
      Real j, const BarrierParameters& params)
  {
    BarrierSampledVolumeTether out;
    if (params.jVolumeTetherWeight <= Real(0))
      return out;
    if (j <= Real(0))
    {
      out.active = true;
      out.feasible = false;
      return out;
    }
    out.active = true;
    const Real logj = std::log(j);
    out.B = Real(0.5) * logj * logj;
    out.Bp = logj / j;
    out.Bpp = (Real(1) - logj) / (j * j);
    return out;
  }

  /// dQ/dF as a vdim x dim matrix. `q` already carries `valid = true`.
  inline Math::SpatialMatrix<Real> dQdF(const BarrierSampledQpt& q)
  {
    const std::size_t d = static_cast<std::size_t>(q.F.rows());
    const Real coeff = Real(2) / (static_cast<Real>(d) * q.D);
    const Real ndt = q.N / static_cast<Real>(d);
    return coeff * (q.MAt - ndt * q.FinvT);
  }

  /// Per-quadrature-point energy density `e(F, A; gamma, w_Q)`.
  inline Real energyDensity(
      const BarrierSampledQpt& q,
      Real gamma, Real domainMeasure,
      const BarrierParameters& params,
      const BarrierSampledQBarrier& bq,
      const BarrierSampledJBarrier& bj,
      const BarrierSampledVolumeTether& bv)
  {
    if (!q.valid || !bq.feasible || !bj.feasible || !bv.feasible)
      return Real(0);
    const Real w_s = gamma / domainMeasure;
    const Real w_Q = params.qBarrierWeight / domainMeasure;
    const Real w_j = params.jBarrierWeight / domainMeasure;
    const Real w_v = params.jVolumeTetherWeight / domainMeasure;
    return w_s * (q.Q - Real(1)) + w_Q * bq.B + w_j * bj.B + w_v * bv.B;
  }

  /// Per-qpt "first Piola" tensor S = dE_density/dF (vdim x dim).
  inline Math::SpatialMatrix<Real> firstPiola(
      const BarrierSampledQpt& q,
      Real gamma, Real domainMeasure,
      const BarrierParameters& params,
      const BarrierSampledQBarrier& bq,
      const BarrierSampledJBarrier& bj,
      const BarrierSampledVolumeTether& bv)
  {
    const Real w_s = gamma / domainMeasure;
    const Real w_Q = params.qBarrierWeight / domainMeasure;
    const Real w_j = params.jBarrierWeight / domainMeasure;
    const Real w_v = params.jVolumeTetherWeight / domainMeasure;
    const Real w_eff = w_s + (bq.active ? w_Q * bq.Bp : Real(0));
    Math::SpatialMatrix<Real> S = w_eff * dQdF(q);
    if (bj.active)
      S += w_j * bj.Bp * q.j * q.FinvT;
    if (bv.active)
      S += w_v * bv.Bp * q.j * q.FinvT;
    return S;
  }

  /**
   * @brief Cell-local sampled energy and residual: standalone primitives.
   *
   * These free functions iterate the quadrature points of one cell,
   * evaluate `evaluateBarrierSampledQpt` and the Q-barrier at each, and
   * accumulate energy or residual contributions. They are used both by
   * the form-language integrator and directly by the FD unit test.
   *
   * Template parameter `StateType` must be a vector GridFunction whose
   * `Variational::Jacobian` supports `.getValue(IntegrationPoint)`.
   */

  /// Per-cell energy E_K(u) = int_K e(F, A) dx (quadrature approximation).
  template <class StateType>
  Real computeBarrierSampledCellEnergy(
      const Geometry::Polytope& cell,
      const StateType& u,
      Real gamma,
      const BarrierParameters& params,
      std::size_t quadratureOrder = 0)
  {
    using Variational::Jacobian;
    using Variational::IntegrationPoint;

    const std::size_t d = cell.getDimension();
    const auto& fes     = u.getFiniteElementSpace();
    const auto& fe      = fes.getFiniteElement(d, cell.getIndex());

    const std::size_t order =
      quadratureOrder > 0
        ? quadratureOrder
        : lsrQuadOrderFor(fe.getOrder());

    const auto& qf =
      QF::PolytopeQuadratureFormula::get(order, cell.getGeometry());
    const auto& quadrature = cell.getQuadrature(qf);

    auto gradU = Jacobian(u);
    Math::SpatialMatrix<Real> A;
    Real energy = Real(0);
    for (std::size_t q = 0; q < quadrature.getSize(); ++q)
    {
      const auto& p  = quadrature.getPoint(q);
      const IntegrationPoint ip(p, &qf, q);
      const auto& rc = qf.getPoint(q);
      cell.getTransformation().jacobian(A, rc);
      const auto qpt = evaluateBarrierSampledQpt(
          gradU.getValue(ip), A, params.jMin);
      if (!qpt.valid) continue;
      const auto bq = evaluateQBarrier(qpt.Q, params);
      const auto bj = evaluateJBarrier(qpt.j, params);
      const auto bv = evaluateVolumeTether(qpt.j, params);
      if (!bq.feasible || !bj.feasible || !bv.feasible) continue;
      const Real wdet = qf.getWeight(q) * p.getDistortion();
      energy += wdet * energyDensity(
          qpt, gamma, params.domainMeasure, params, bq, bj, bv);
    }
    return energy;
  }

  /// Per-cell residual R^e[i] = int_K (S : jac_i) dx (quadrature approximation).
  /// The returned vector has `fe.getCount()` entries — one per FE local dof.
  template <class StateType>
  Math::Vector<Real> computeBarrierSampledCellResidual(
      const Geometry::Polytope& cell,
      const StateType& u,
      Real gamma,
      const BarrierParameters& params,
      std::size_t quadratureOrder = 0)
  {
    using Variational::Jacobian;
    using Variational::IntegrationPoint;

    const std::size_t d = cell.getDimension();
    const auto& fes     = u.getFiniteElementSpace();
    const auto& fe      = fes.getFiniteElement(d, cell.getIndex());
    const std::size_t nte = fe.getCount();

    const std::size_t order =
      quadratureOrder > 0
        ? quadratureOrder
        : lsrQuadOrderFor(fe.getOrder());

    const auto& qf =
      QF::PolytopeQuadratureFormula::get(order, cell.getGeometry());
    const auto& quadrature = cell.getQuadrature(qf);

    auto gradU = Jacobian(u);
    Math::SpatialMatrix<Real> A;
    Math::SpatialMatrix<Real> jac_i(
        static_cast<Eigen::Index>(d), static_cast<Eigen::Index>(d));

    Math::Vector<Real> residual =
      Math::Vector<Real>::Zero(static_cast<Eigen::Index>(nte));

    for (std::size_t q = 0; q < quadrature.getSize(); ++q)
    {
      const auto& p  = quadrature.getPoint(q);
      const IntegrationPoint ip(p, &qf, q);
      const auto& rc = qf.getPoint(q);
      cell.getTransformation().jacobian(A, rc);
      const auto qpt = evaluateBarrierSampledQpt(
          gradU.getValue(ip), A, params.jMin);
      if (!qpt.valid) continue;
      const auto bq = evaluateQBarrier(qpt.Q, params);
      const auto bj = evaluateJBarrier(qpt.j, params);
      const auto bv = evaluateVolumeTether(qpt.j, params);
      if (!bq.feasible || !bj.feasible || !bv.feasible) continue;

      const Real wdet = qf.getWeight(q) * p.getDistortion();
      const Math::SpatialMatrix<Real> S =
        firstPiola(qpt, gamma, params.domainMeasure, params, bq, bj, bv);

      for (std::size_t i = 0; i < nte; ++i)
      {
        const auto& basis_i = fe.getBasis(i);
        for (std::size_t r = 0; r < d; ++r)
          for (std::size_t c = 0; c < d; ++c)
            jac_i(r, c) = basis_i.template getDerivative<1>(r, c)(rc);
        jac_i *= p.getJacobianInverse();
        residual(static_cast<Eigen::Index>(i)) +=
          wdet * S.dot(jac_i);
      }
    }
    return residual;
  }

  /**
   * @brief Linear-form integrator: sampled shape-quality residual.
   *
   * FES-agnostic. Works for any Pk vector displacement space where
   * `Variational::Jacobian(state)` is available and the test FE exposes
   * `getBasis(i).getDerivative<1>(r,c)(rc)`.
   */
  template <class GammaDerived, class TestType, class StateType>
  class BarrierResidualIntegratorSampled final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using GammaType = Variational::RealFunctionBase<GammaDerived>;

      BarrierResidualIntegratorSampled(
          const GammaType& gamma,
          const TestType& v,
          const StateType& u,
          BarrierParameters params,
          std::size_t quadratureOrder = 0)
        : Variational::LinearFormIntegratorBase<Real>(v),
          m_gamma(gamma.copy()),
          m_test(v), m_state(u),
          m_params(params),
          m_quadratureOrder(quadratureOrder)
      {}

      BarrierResidualIntegratorSampled(const BarrierResidualIntegratorSampled& other)
        : Variational::LinearFormIntegratorBase<Real>(other),
          m_gamma(other.m_gamma ? other.m_gamma->copy() : nullptr),
          m_test(other.m_test), m_state(other.m_state),
          m_params(other.m_params),
          m_quadratureOrder(other.m_quadratureOrder),
          m_polytope(other.m_polytope),
          m_elemVec()
      {}

      BarrierResidualIntegratorSampled& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const std::size_t d = polytope.getDimension();
        const Index idx     = polytope.getIndex();

        const auto& testfes = m_test.get().getFiniteElementSpace();
        const auto& testfe  = testfes.getFiniteElement(d, idx);
        const std::size_t nte = testfe.getCount();

        m_elemVec.resize(static_cast<Eigen::Index>(nte));
        m_elemVec.setZero();

        // Evaluate gamma at the cell centroid (consistent with the
        // closed-form prototype; the sampled energy itself is fully
        // quadrature-integrated, so gamma spatial variation would only
        // need a qpt evaluation in computeBarrierSampledCellResidual —
        // for now a centroid sample matches the existing convention).
        const Geometry::Point centroid(
            polytope, detail::triangleReferenceCentroid());
        const Real gammaVal = m_gamma->getValue(centroid);

        const auto R = computeBarrierSampledCellResidual(
            polytope, m_state.get(), gammaVal, m_params, m_quadratureOrder);

        if (R.size() != static_cast<Eigen::Index>(nte))
        {
          // Singular cell or dof-count mismatch; count and skip.
          barrierInadmissibleCount().fetch_add(1, std::memory_order_relaxed);
          return *this;
        }
        m_elemVec = R;
        return *this;
      }

      ScalarType integrate(std::size_t te) final override
      { return m_elemVec(te); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      BarrierResidualIntegratorSampled* copy() const noexcept final override
      { return new BarrierResidualIntegratorSampled(*this); }

    private:
      std::unique_ptr<GammaType> m_gamma;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_state;
      BarrierParameters m_params;
      std::size_t m_quadratureOrder;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  /**
   * @brief Per-cell tangent matrix assembled by quadrature.
   *
   * Returns an (nte x ntr) matrix of contributions
   *   K^e[i,j] = int_K [ w_eff * d2Q/dF^2.(jac_i, jac_j)
   *                     + (w_Q/|Omega|) * B_Q''(Q) * <dQ/dF, jac_i> <dQ/dF, jac_j>] dx.
   *
   * Indices of the returned matrix: row = test dof, col = trial dof.
   * When a qpt is inadmissible, a large diagonal penalty of
   * `kBarrierSingularPenalty` is contributed to the diagonal of ALL
   * local dofs for that qpt — same policy as the closed-form.
   */
  template <class StateType>
  Math::Matrix<Real> computeBarrierSampledCellTangent(
      const Geometry::Polytope& cell,
      const StateType& u,
      Real gamma,
      const BarrierParameters& params,
      std::size_t quadratureOrder = 0,
      bool projectPSD = false)
  {
    using Variational::Jacobian;
    using Variational::IntegrationPoint;

    const std::size_t d   = cell.getDimension();
    const auto& fes       = u.getFiniteElementSpace();
    const auto& fe        = fes.getFiniteElement(d, cell.getIndex());
    const std::size_t n   = fe.getCount();
    const Real rd         = static_cast<Real>(d);

    const std::size_t order =
      quadratureOrder > 0
        ? quadratureOrder
        : lsrQuadOrderFor(fe.getOrder());

    const auto& qf =
      QF::PolytopeQuadratureFormula::get(order, cell.getGeometry());
    const auto& quadrature = cell.getQuadrature(qf);

    auto gradU = Jacobian(u);
    Math::SpatialMatrix<Real> A;
    Math::Matrix<Real> K =
      Math::Matrix<Real>::Zero(
          static_cast<Eigen::Index>(n),
          static_cast<Eigen::Index>(n));

    bool anyInadmissible = false;

    for (std::size_t q = 0; q < quadrature.getSize(); ++q)
    {
      const auto& p  = quadrature.getPoint(q);
      const IntegrationPoint ip(p, &qf, q);
      const auto& rc = qf.getPoint(q);
      cell.getTransformation().jacobian(A, rc);

      const auto qpt = evaluateBarrierSampledQpt(
          gradU.getValue(ip), A, params.jMin);
      if (!qpt.valid)
      {
        anyInadmissible = true;
        continue;
      }

      const auto bq = evaluateQBarrier(qpt.Q, params);
      const auto bj = evaluateJBarrier(qpt.j, params);
      const auto bv = evaluateVolumeTether(qpt.j, params);
      if (!bq.feasible || !bj.feasible || !bv.feasible)
      {
        anyInadmissible = true;
        continue;
      }

      const Real wdet = qf.getWeight(q) * p.getDistortion();
      const Real w_s  = gamma / params.domainMeasure;
      const Real w_Q  = params.qBarrierWeight / params.domainMeasure;
      const Real w_j  = params.jBarrierWeight / params.domainMeasure;
      const Real w_v  = params.jVolumeTetherWeight / params.domainMeasure;
      const Real w_eff = w_s + (bq.active ? w_Q * bq.Bp : Real(0));

      // Basis-function physical-space Jacobian per local dof.
      // jac[i] is a (d x d) SpatialMatrix — row = component, col = direction.
      std::vector<Math::SpatialMatrix<Real>> jac(n,
          Math::SpatialMatrix<Real>(
              static_cast<Eigen::Index>(d),
              static_cast<Eigen::Index>(d)));
      for (std::size_t i = 0; i < n; ++i)
      {
        const auto& basis_i = fe.getBasis(i);
        for (std::size_t r = 0; r < d; ++r)
          for (std::size_t c = 0; c < d; ++c)
            jac[i](r, c) = basis_i.template getDerivative<1>(r, c)(rc);
        jac[i] *= p.getJacobianInverse();
      }

      // Precompute per-basis contractions with dQ/dF-related matrices.
      // MAts[i] = <F, jac_i>  (scalar)
      // FinvTs[i] = <F⁻ᵀ, jac_i>  (scalar)
      // MAts_vec[i] = F · jac_i  (would be meaningless — we need dots)
      // For the general-d Hessian formula we also need MAtjac and FinvTjac
      // contracted into jac_j via a matrix product, so we precompute those.
      //
      //   Term 1: (2/(d D)) <H1, H2>
      //   Term 2: (2/(d D))(N/d) <H1, F⁻ᵀ H2ᵀ F⁻ᵀ>
      //   Terms 3/4: −(4/(d^2 D)) [<F,H1><F⁻ᵀ,H2> + <F⁻ᵀ,H1><F,H2>]
      //   Term 5: +(4N/(d^3 D)) <F⁻ᵀ,H1><F⁻ᵀ,H2>
      //
      // Plus the Q-barrier rank-1 term: w_Q B_Q'' <dQ/dF,H1><dQ/dF,H2>.

      const Real twoDinvD   = Real(2) / (rd * qpt.D);
      const Real fourd2invD = Real(4) / (rd * rd * qpt.D);
      const Real nod        = qpt.N / rd;

      const Math::SpatialMatrix<Real> GdQ = dQdF(qpt); // (vdim x dim)

      // Precompute scalar contractions per basis.
      const Math::SpatialMatrix<Real> Finv = qpt.FinvT.transpose();
      const Math::SpatialMatrix<Real> Gdet = qpt.j * qpt.FinvT;

      std::vector<Real> sc_MAt(n), sc_FinvT(n), sc_dQ(n), sc_det(n);
      std::vector<Math::SpatialMatrix<Real>> FinvJac(n,
          Math::SpatialMatrix<Real>(
              static_cast<Eigen::Index>(d),
              static_cast<Eigen::Index>(d)));
      for (std::size_t i = 0; i < n; ++i)
      {
        sc_MAt[i]   = qpt.MAt.dot(jac[i]);
        sc_FinvT[i] = qpt.FinvT.dot(jac[i]);
        sc_dQ[i]    = GdQ.dot(jac[i]);
        sc_det[i]   = Gdet.dot(jac[i]);
        FinvJac[i]  = Finv * jac[i];
      }

      // Precompute per-basis contracted vectors for Terms 1 and 2.
      // Term 1 inner product: <jac_i, jac_j C> = Tr(jac_i^T jac_j C).
      // Term 2 inner product: <jac_i, F⁻ᵀ jac_j^T F⁻ᵀ>.
      //   Let v_j = F⁻ᵀ (col j of jac_j) for each direction...
      //   Actually: F⁻ᵀ jac_j^T F⁻ᵀ = (F⁻¹ jac_j F⁻ᵀ)^T — no simpler.
      //   We use the matrix inner product directly:
      //   <jac_i, F⁻ᵀ jac_j^T F⁻ᵀ> = sum_aA jac_i_aA (F⁻ᵀ jac_j^T F⁻ᵀ)_aA
      //     = Tr(jac_i^T F⁻ᵀ jac_j^T F⁻ᵀ).
      //   For the accumulation, we precompute FinvT_jac_i = F⁻ᵀ jac_i^T,
      //   then Term2[i,j] = Tr(FinvT_jac_i^T jac_j F⁻ᵀ^T) = (FinvT_jac_i).dot(jac_j * FinvT).
      //   This avoids a full rank-4 tensor.

      // FiTji[i] = F⁻ᵀ · jac_i^T  (dim x dim matrix)
      std::vector<Math::SpatialMatrix<Real>> FiTji(n,
          Math::SpatialMatrix<Real>(
              static_cast<Eigen::Index>(d),
              static_cast<Eigen::Index>(d)));
      for (std::size_t i = 0; i < n; ++i)
        FiTji[i] = qpt.FinvT * jac[i].transpose();

      // jacC[i] = jac_i · C  (jac_i: vdim x dim, C: dim x dim → vdim x dim)
      std::vector<Math::SpatialMatrix<Real>> jacC(n,
          Math::SpatialMatrix<Real>(
              static_cast<Eigen::Index>(d),
              static_cast<Eigen::Index>(d)));
      for (std::size_t i = 0; i < n; ++i)
        jacC[i] = jac[i] * qpt.C;

      // jacFinvT[i] = jac_i · F⁻ᵀ^T = jac_i · F⁻¹  (for Term 2 contract)
      // Actually I need: Tr(FiTji[j] · jac_i^T) for term-2 with roles (i,j)
      // which is FiTji[j].dot(jac[i]).
      // Symmetrically: Tr(FiTji[i] · jac_j^T) = FiTji[i].dot(jac[j]).
      // Note: Tr(A B^T) = A.dot(B) for matrices, so:
      //   Tr(jac_i^T F⁻ᵀ jac_j^T F⁻ᵀ)
      //     = Tr((F⁻ᵀ jac_i jac_j^T)^T F⁻ᵀ)  -- no...
      //   More carefully: (F⁻ᵀ jac_j^T F⁻ᵀ)_{aA} = sum_{bB} F⁻ᵀ_{ab} jac_j^T_{bc} F⁻ᵀ_{cA}
      //     = sum_{bc} FinvT_{ab} jac_j_{cb} FinvT_{cA}.
      //   <jac_i, F⁻ᵀ jac_j^T F⁻ᵀ> = sum_{aA} jac_i_{aA} sum_{bc} FinvT_{ab} jac_j_{cb} FinvT_{cA}
      //     = sum_{b,A,c} (sum_a jac_i_{aA} FinvT_{ab}) jac_j_{cb} FinvT_{cA}
      //     = sum_{b,A,c} (jac_i^T FinvT)_{Ab} jac_j_{cb} FinvT_{cA}
      //     = sum_{b,c} (jac_i^T FinvT)_{·b} · (jac_j FinvT^T)_{c·} contracted on b... messy.
      //   Simplest: just compute FiTji[j].dot(jac[i]) by construction:
      //     FiTji[j] = F⁻ᵀ jac_j^T  (dim x dim)
      //     FiTji[j].dot(jac[i]) = sum_{rA} (FinvT jac_j^T)_{rA} jac_i_{rA}  -- wrong shape
      //   Actually FiTji[i] (dim x dim):  F⁻ᵀ jac_i^T.
      //   Want: sum_{aA} jac_i_{aA} (F⁻ᵀ jac_j^T F⁻ᵀ)_{aA}
      //       = Tr( jac_i^T · F⁻ᵀ · jac_j^T · F⁻ᵀ )
      //       = (jac_i^T F⁻ᵀ) : (jac_j^T F⁻ᵀ)  since Tr(X Y) = X^T . Y  ??? no.
      //   Tr(A B) = Tr(B A) cyclic. Let A = jac_i^T F⁻ᵀ, B = jac_j^T F⁻ᵀ.
      //   Tr(A B) = A.dot(B^T) = (jac_i^T F⁻ᵀ).dot(F⁻¹ jac_j) = FiTji[i]^T.dot(F⁻¹ jac_j).
      //   Even simpler: precompute JiFiT[i] = jac_i^T · F⁻ᵀ (dim x dim),
      //   then Term2 contribution = JiFiT[i].dot(JiFiT[j]^T)... no.
      //   Final answer: Tr(jac_i^T F⁻ᵀ jac_j^T F⁻ᵀ) =
      //     (F⁻ᵀ jac_i · F⁻¹ jac_j^T) trace? Still no.
      //   Use cyclic: Tr(jac_i^T F⁻ᵀ jac_j^T F⁻ᵀ) = Tr(F⁻ᵀ jac_i^T · F⁻ᵀ jac_j^T)^T wait.
      //   No. Tr(ABCD) = Tr(BCDA). So Tr(jac_i^T F⁻ᵀ jac_j^T F⁻ᵀ)
      //     = Tr(F⁻ᵀ jac_i^T F⁻ᵀ jac_j^T).
      //   Let P = F⁻ᵀ jac_i^T (= FiTji[i]) and Q = F⁻ᵀ jac_j^T (= FiTji[j]).
      //   = Tr(P Q) = P.dot(Q^T) by definition of Tr(AB) = sum_ij A_ij B_ji.
      //   So Term2[i,j] = FiTji[i].dot(FiTji[j].transpose()).
      //   Let me verify: P = FiTji[i] = F⁻ᵀ jac_i^T (dim x dim).
      //                  Q = FiTji[j] = F⁻ᵀ jac_j^T (dim x dim).
      //   Tr(P Q) = Tr(F⁻ᵀ jac_i^T F⁻ᵀ jac_j^T). ✓
      //   P.dot(Q^T) = sum_{rc} P_{rc} Q^T_{rc} = sum_{rc} P_{rc} Q_{cr}
      //              = sum_{rc} P_{rc} Q_{cr} = Tr(P Q) ✓.

      for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i; j < n; ++j)
        {
          // Term 1: (2/(d D)) <jac_i, jac_j C>
          const Real t1 = twoDinvD * jac[i].dot(jacC[j]);

          // Term 2: (2/(d D))(N/d) Tr(F⁻ᵀ jac_i^T F⁻ᵀ jac_j^T) = nod * FiTji[i].dot(FiTji[j]^T)
          const Real t2 = twoDinvD * nod * FiTji[i].dot(FiTji[j].transpose());

	          // Terms 3+4: −(4/(d^2 D)) [<F,jac_i><F⁻ᵀ,jac_j> + sym]
          const Real t34 = -fourd2invD
              * (sc_MAt[i] * sc_FinvT[j] + sc_FinvT[i] * sc_MAt[j]);

          // Term 5: +(4N/(d^3 D)) <F⁻ᵀ,jac_i><F⁻ᵀ,jac_j>
          const Real t5 = (fourd2invD / rd) * qpt.N * sc_FinvT[i] * sc_FinvT[j];

          // Q-barrier rank-1 term.
          const Real tQ = bq.active
              ? w_Q * bq.Bpp * sc_dQ[i] * sc_dQ[j]
              : Real(0);

          Real tJ = Real(0);
          if (bj.active)
          {
            const Real d2det = qpt.j
              * (sc_FinvT[i] * sc_FinvT[j]
                 - (FinvJac[i] * FinvJac[j]).trace());
            tJ = w_j * (bj.Bpp * sc_det[i] * sc_det[j] + bj.Bp * d2det);
          }
          if (bv.active)
          {
            const Real d2det = qpt.j
              * (sc_FinvT[i] * sc_FinvT[j]
                 - (FinvJac[i] * FinvJac[j]).trace());
            tJ += w_v * (bv.Bpp * sc_det[i] * sc_det[j] + bv.Bp * d2det);
          }

          const Real H = wdet * (w_eff * (t1 + t2 + t34 + t5) + tQ + tJ);
          K(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) += H;
          if (i != j)
            K(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) += H;
        }
    }

    if (anyInadmissible)
    {
      // Freeze all local dofs: large diagonal penalty, zero off-diagonal.
      K.setZero();
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(n); ++i)
        K(i, i) = kBarrierSingularPenalty;
    }

    return projectPSD ? projectSymmetricMatrixPSD(K) : K;
  }

  /**
   * @brief Bilinear-form integrator: sampled shape-quality tangent.
   */
  template <class GammaDerived, class TrialType, class TestType, class StateType>
  class BarrierTangentIntegratorSampled final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using GammaType = Variational::RealFunctionBase<GammaDerived>;

      BarrierTangentIntegratorSampled(
          const GammaType& gamma,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          BarrierParameters params,
          std::size_t quadratureOrder = 0,
          bool projectPSD = false)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_gamma(gamma.copy()),
          m_trial(du), m_test(v), m_state(u),
          m_params(params),
          m_quadratureOrder(quadratureOrder),
          m_projectPSD(projectPSD)
      {}

      BarrierTangentIntegratorSampled(const BarrierTangentIntegratorSampled& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_gamma(other.m_gamma ? other.m_gamma->copy() : nullptr),
          m_trial(other.m_trial), m_test(other.m_test), m_state(other.m_state),
          m_params(other.m_params),
          m_quadratureOrder(other.m_quadratureOrder),
          m_projectPSD(other.m_projectPSD),
          m_polytope(other.m_polytope),
          m_matrix()
      {}

      BarrierTangentIntegratorSampled& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const std::size_t d = polytope.getDimension();
        const Index idx     = polytope.getIndex();

        const auto& testfes = m_test.get().getFiniteElementSpace();
        const auto& testfe  = testfes.getFiniteElement(d, idx);
        const std::size_t nte = testfe.getCount();

        m_matrix.resize(
            static_cast<Eigen::Index>(nte),
            static_cast<Eigen::Index>(nte));
        m_matrix.setZero();

        const Geometry::Point centroid(
            polytope, detail::triangleReferenceCentroid());
        const Real gammaVal = m_gamma->getValue(centroid);

        m_matrix = computeBarrierSampledCellTangent(
            polytope, m_state.get(), gammaVal, m_params, m_quadratureOrder,
            m_projectPSD);

        return *this;
      }

      ScalarType integrate(std::size_t tr, std::size_t te) final override
      { return m_matrix(te, tr); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      BarrierTangentIntegratorSampled* copy() const noexcept final override
      { return new BarrierTangentIntegratorSampled(*this); }

    private:
      std::unique_ptr<GammaType> m_gamma;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType>  m_test;
      std::reference_wrapper<const StateType> m_state;
      BarrierParameters m_params;
      std::size_t m_quadratureOrder;
      bool m_projectPSD;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  /**
   * @brief High-level helper for the sampled shape-quality energy.
   *
   * Mirrors `JacobianAdmissibilityBarrier`'s API.
   */
  template <class Trial, class Test, class State>
  class JacobianAdmissibilityBarrierSampled
  {
    public:
      JacobianAdmissibilityBarrierSampled(
          const Trial& du, const Test& v, const State& u,
          std::size_t quadratureOrder = 0)
        : m_du(du), m_v(v), m_u(u), m_quadratureOrder(quadratureOrder)
      {}

      template <class GammaDerived>
      auto Tangent(
          const Variational::RealFunctionBase<GammaDerived>& gamma,
          BarrierParameters params) const
      {
        return BarrierTangentIntegratorSampled<
                   GammaDerived, Trial, Test, State>(
              gamma, m_du.get(), m_v.get(), m_u.get(),
	              params, m_quadratureOrder);
      }

      template <class GammaDerived>
      auto TangentPSDProjected(
          const Variational::RealFunctionBase<GammaDerived>& gamma,
          BarrierParameters params) const
      {
        return BarrierTangentIntegratorSampled<
                   GammaDerived, Trial, Test, State>(
              gamma, m_du.get(), m_v.get(), m_u.get(),
              params, m_quadratureOrder, true);
      }

      template <class GammaDerived>
      auto Residual(
          const Variational::RealFunctionBase<GammaDerived>& gamma,
          BarrierParameters params) const
      {
        return BarrierResidualIntegratorSampled<
                   GammaDerived, Test, State>(
              gamma, m_v.get(), m_u.get(), params, m_quadratureOrder);
      }

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
      std::size_t m_quadratureOrder;
  };

  template <class Trial, class Test, class State>
  JacobianAdmissibilityBarrierSampled(
      const Trial&, const Test&, const State&, std::size_t = 0)
    -> JacobianAdmissibilityBarrierSampled<Trial, Test, State>;
}

#endif
