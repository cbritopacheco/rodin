/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_LSRINTEGRATORS_H
#define RODIN_ADAPTATION_LSRINTEGRATORS_H

/**
 * @file
 * @brief Form-language integrators for the LSR (level-set registration)
 *        data term.
 *
 * The level set is described by three Rodin form-language objects of
 * equal status:
 *
 *   - phi:  `Variational::RealFunctionBase<...>`     (always required)
 *   - grad: `Variational::VectorFunctionBase<Real, ...>` (always required)
 *   - hess: `Variational::MatrixFunctionBase<Real, ...>` (Newton mode only)
 *
 * The integrator stores each as `std::unique_ptr<...Base>` after
 * `.copy()`, mirroring `LinearElasticityIntegrator`'s storage of
 * `lambda` and `mu`.
 *
 * Evaluation.
 *   The LSR penalty needs `phi(X + u_h(X))`. At each quadrature point of
 *   the parent mesh we construct a `Geometry::Point` whose physical
 *   coordinates are `X + u_h(X)` and whose reference coordinates are
 *   obtained by localising that physical point in the background mesh,
 *   then call `phi->getValue(movedPoint)`,
 *   `grad->getValue(movedPoint)`, `hess->getValue(movedPoint)`. The
 *   integrator does NOT invoke `Grad(phi)` or `Hess(phi)` itself — the
 *   derivatives are the user's inputs.
 *
 * The LSR tangent flavour is selected at the TYPE level via the
 * non-type template parameter `LSRIntegratorTangentMode`.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/MatrixFunction.h"
#include "Rodin/Variational/RealFunction.h"
#include "Rodin/Variational/Flow.h"
#include "Rodin/Variational/VectorFunction.h"

namespace Rodin::Adaptation
{
  enum class LSRIntegratorTangentMode
  {
    GaussNewton,
    Newton,
    PSDProjectedNewton  ///< Full Newton with the per-quadrature-point
                        ///< correction `r * hess(phi)` clamped to its
                        ///< positive-semidefinite part (eigenvalue >= 0).
                        ///< This removes the indefinite component of the
                        ///< second-order correction; it is therefore a
                        ///< safeguarded Newton model, not the exact Newton
                        ///< tangent whenever the projection is active.
  };

  inline const char* lsrIntegratorTangentModeName(LSRIntegratorTangentMode m)
  {
    switch (m)
    {
      case LSRIntegratorTangentMode::GaussNewton:        return "GaussNewton";
      case LSRIntegratorTangentMode::Newton:             return "Newton";
      case LSRIntegratorTangentMode::PSDProjectedNewton: return "PSDProjectedNewton";
    }
    return "?";
  }

  namespace detail
  {
    /**
     * @brief Project a symmetric 2x2 matrix onto its positive-semidefinite
     *        part by eigenvalue clamping.
     *
     * Given M with spectral decomposition M = lp P+ + lm P-, returns
     * `(M)+ = max(0, lp) P+ + max(0, lm) P-`. Using the spectral identity
     * `lp P+ = (lp / (lp - lm)) (M - lm I)` (and the symmetric variant for
     * lm < 0 < lp), the projection is computed without an explicit
     * eigendecomposition.
     */
    inline Math::SpatialMatrix<Real> psdProject2x2(
        const Math::SpatialMatrix<Real>& M)
    {
      const Real a = M(0, 0);
      const Real c = M(1, 1);
      const Real b = Real(0.5) * (M(0, 1) + M(1, 0));   ///< symmetric part
      const Real tr = a + c;
      const Real disc = (a - c) * (a - c) * Real(0.25) + b * b;
      const Real s = std::sqrt(std::max(disc, Real(0)));
      const Real lp = Real(0.5) * tr + s;               ///< larger eigenvalue
      const Real lm = Real(0.5) * tr - s;               ///< smaller eigenvalue

      Math::SpatialMatrix<Real> out(2, 2);
      if (lp <= Real(0))
      {
        out.setZero();
        return out;
      }
      if (lm >= Real(0))
      {
        out = M;
        // Re-symmetrise to scrub off any numerical asymmetry from caller.
        const Real off = Real(0.5) * (out(0, 1) + out(1, 0));
        out(0, 1) = off;
        out(1, 0) = off;
        return out;
      }
      // Mixed signs: lm < 0 < lp. Keep only the lp eigenspace.
      // M+ = (lp / (lp - lm)) * (M - lm I).
      const Real scale = lp / (lp - lm);
      out = scale * M;
      out(0, 0) -= scale * lm;
      out(1, 1) -= scale * lm;
      const Real off = Real(0.5) * (out(0, 1) + out(1, 0));
      out(0, 1) = off;
      out(1, 0) = off;
      return out;
    }
  }

  struct LSRIntegratorParameters
  {
    enum class FieldEvaluation
    {
      PhysicalPoint,
      FlowTrace
    };

    Real rhoS = 1;
    Real deltaW = 0;
    Real hRef = 0;
    Real normalizer = 1;
    std::size_t quadratureOrder = 0; ///< 0 selects the default FE-based rule.
    FieldEvaluation fieldEvaluation = FieldEvaluation::PhysicalPoint;
  };

  inline std::size_t lsrQuadOrderFor(std::size_t feOrder)
  {
    // Need (#qpts * vdim) >= n_dof_per_cell with a comfortable margin
    // so the basis-at-qpt evaluation map B is well-conditioned on every
    // cell DOF, including P2 edge-midpoint internal modes. The
    // polynomial-exactness lower bound is 2*feOrder; 4*feOrder gives an
    // overdetermined B on P2 vector cells without changing P1.
    return std::max<std::size_t>(2, 4 * feOrder);
  }

  inline std::size_t lsrQuadOrderFor(
      std::size_t feOrder, const LSRIntegratorParameters& params)
  {
    return params.quadratureOrder > 0
      ? params.quadratureOrder
      : lsrQuadOrderFor(feOrder);
  }

  namespace detail
  {
    inline Real referenceCellMargin(
        Geometry::Polytope::Type geometry,
        const Math::SpatialPoint& rc,
        std::size_t& mostViolatedFace)
    {
      const Geometry::Polytope::Traits traits(geometry);
      const auto& hs = traits.getHalfSpace();
      Real margin = std::numeric_limits<Real>::infinity();
      mostViolatedFace = 0;
      for (std::size_t j = 0; j < static_cast<std::size_t>(hs.vector.size()); ++j)
      {
        const Real phi = hs.vector[j] - rc.dot(hs.matrix.row(j).transpose());
        if (phi < margin)
        {
          margin = phi;
          mostViolatedFace = j;
        }
      }
      return margin;
    }

    /// Build a Geometry::Point at the deformed location y = X + u(X).
    ///
    /// GridFunction evaluation is reference-coordinate based when the point
    /// belongs to the same mesh. Consequently the moved physical point must
    /// be localised before it is passed to phi, grad, or hess. The routine
    /// first inverts the original cell and then walks through adjacent cells
    /// through the most violated reference half-space. If localisation fails
    /// at the boundary, it returns the parent-cell extrapolation.
    inline Geometry::Point makeMovedPoint(
        const Geometry::Point& original,
        const Math::SpatialVector<Real>& yPhysical)
    {
      const auto& parent = original.getPolytope();
      const auto& mesh = parent.getMesh();
      const std::size_t cd = parent.getDimension();
      Index cell = parent.getIndex();
      Math::SpatialPoint movedReference;
      const Real tol = Real(64) * std::numeric_limits<Real>::epsilon();

      for (std::size_t hop = 0; hop < 64; ++hop)
      {
        mesh.getPolytopeTransformation(cd, cell).inverse(
            movedReference, yPhysical);
        std::size_t mostViolatedFace = 0;
        const Real margin = referenceCellMargin(
            mesh.getGeometry(cd, cell), movedReference, mostViolatedFace);
        if (margin >= -tol)
        {
          const auto it = mesh.getPolytope(cd, cell);
          return Geometry::Point(*it, movedReference, yPhysical);
        }

        const auto& conn = mesh.getConnectivity();
        const auto& faces = conn.getIncidence(cd, cd - 1).at(cell);
        if (mostViolatedFace >= faces.size())
          break;

        const Index face = faces[mostViolatedFace];
        if (mesh.isBoundary(face))
          break;

        const auto& nbrs = conn.getIncidence(cd - 1, cd).at(face);
        if (nbrs.size() != 2)
          break;

        const Index next = (nbrs[0] == cell) ? nbrs[1] : nbrs[0];
        if (next == cell)
          break;
        cell = next;
      }

      parent.getTransformation().inverse(movedReference, yPhysical);
      return Geometry::Point(
          parent,
          movedReference,
          yPhysical);
    }

    template <class Phi>
    struct TracedPoint
    {
      Geometry::Point point;
      Real correction = 0;
      bool exited = false;
    };

    template <class Phi>
    TracedPoint<Phi> traceMovedPoint(
        const Geometry::Point& original,
        const Math::SpatialVector<Real>& displacement,
        const Phi& phi,
        LSRIntegratorParameters::FieldEvaluation mode)
    {
      if (mode == LSRIntegratorParameters::FieldEvaluation::PhysicalPoint)
      {
        Math::SpatialVector<Real> y(original.getPolytope().getDimension());
        for (std::size_t c = 0; c < static_cast<std::size_t>(y.size()); ++c)
          y(c) = original.getCoordinates()(c) + displacement(c);
        return { makeMovedPoint(original, y), Real(0), false };
      }

      auto velocity =
        [displacement](const Geometry::Point&)
        {
          return displacement;
        };
      Variational::Flow flow(Real(1), phi, velocity);
      const auto trace = flow.trace(original);
      return { trace.getPoint(), trace.getCorrection(), trace.exited() };
    }
  }

  /**
   * @brief Linear-form integrator: LSR residual.
   *
   *   R[te] = sum_q wq * |J_q| * rhoS * W(s_q) * normalizer
   *           * (phi(y_q) - s_q) * (grad(y_q) . v_te).
   *
   * Both `phi` and `grad` are evaluated at the deformed point y_q.
   */
  template <class PhiDerived, class GradDerived,
            class PsiType, class TestType, class StateType>
  class LSRResidualIntegrator final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using PhiType  = Variational::RealFunctionBase<PhiDerived>;
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      LSRResidualIntegrator(
          const PhiType& phi, const GradType& grad,
          const PsiType& psi,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LinearFormIntegratorBase<Real>(v),
          m_phi(phi.copy()), m_grad(grad.copy()),
          m_psi(psi), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRResidualIntegrator(const LSRResidualIntegrator& other)
        : Variational::LinearFormIntegratorBase<Real>(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_psi(other.m_psi), m_test(other.m_test), m_state(other.m_state),
          m_params(other.m_params),
          m_polytope(other.m_polytope),
          m_elemVec()
      {}

      LSRResidualIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const std::size_t d = polytope.getDimension();
        const Index cellIdx = polytope.getIndex();
        const auto& testFES = m_test.get().getFiniteElementSpace();
        const auto& testFE = testFES.getFiniteElement(d, cellIdx);
        const std::size_t testDofs = testFE.getCount();
        const std::size_t vdim = testFES.getVectorDimension();

        m_elemVec.resize(testDofs);
        m_elemVec.setZero();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(
              lsrQuadOrderFor(testFE.getOrder(), m_params),
              polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const std::size_t nqp = quadrature.getSize();

        for (std::size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const Real wq = qf.getWeight(q);
          const Real distortion = pt.getDistortion();
          const auto& rc = pt.getReferenceCoordinates();

          const Real s = m_psi.get().getValue(ip);
          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> displacement(d);
          displacement.setZero();
          for (std::size_t c = 0; c < vdim; ++c)
            displacement(c) = uqRange(c);

          const auto traced =
            detail::traceMovedPoint(
                pt, displacement, *m_phi, m_params.fieldEvaluation);
          if (traced.exited)
            continue;

          const Real phi_y =
            m_phi->getValue(traced.point) + traced.correction;
          const auto gradPhi = m_grad->getValue(traced.point);

          const Real r = phi_y - s;
          const Real weight =
            std::exp(-s * s / (2 * m_params.deltaW * m_params.deltaW));
          const Real coef =
            m_params.rhoS * weight * r * m_params.normalizer;
          const Real measure = wq * distortion;

          for (std::size_t te = 0; te < testDofs; ++te)
          {
            const auto testValue = testFE.getBasis(te)(rc);
            Real dot = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              dot += gradPhi(c) * testValue(c);
            m_elemVec(te) += measure * coef * dot;
          }
        }
        return *this;
      }

      ScalarType integrate(std::size_t te) final override
      { return m_elemVec(te); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      LSRResidualIntegrator* copy() const noexcept final override
      { return new LSRResidualIntegrator(*this); }

    private:
      std::unique_ptr<PhiType>  m_phi;
      std::unique_ptr<GradType> m_grad;
      std::reference_wrapper<const PsiType> m_psi;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_state;
      LSRIntegratorParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  template <class P, class G, class S, class V, class U>
  LSRResidualIntegrator(
      const Variational::RealFunctionBase<P>&,
      const Variational::VectorFunctionBase<Real, G>&,
      const S&, const V&, const U&, LSRIntegratorParameters)
    -> LSRResidualIntegrator<P, G, S, V, U>;

  template <LSRIntegratorTangentMode Mode,
            class PhiDerived, class GradDerived, class HessDerived,
            class PsiType, class TrialType, class TestType, class StateType>
  class LSRTangentIntegrator;

  // ---- GaussNewton specialisation ------------------------------------------
  //
  //   K_GN[te, tr] = sum_q wq * |J_q| * rhoS * W(s_q) * normalizer
  //                  * (grad(y_q) . v_te) * (grad(y_q) . u_tr).
  template <class PhiDerived, class GradDerived, class HessDerived,
            class PsiType, class TrialType, class TestType, class StateType>
  class LSRTangentIntegrator<
            LSRIntegratorTangentMode::GaussNewton,
            PhiDerived, GradDerived, HessDerived,
            PsiType, TrialType, TestType, StateType>
    final : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using PhiType  = Variational::RealFunctionBase<PhiDerived>;
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;
      static constexpr LSRIntegratorTangentMode Mode = LSRIntegratorTangentMode::GaussNewton;

      LSRTangentIntegrator(
          const PhiType& phi,
          const GradType& grad,
          const PsiType& psi,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_phi(phi.copy()), m_grad(grad.copy()),
          m_psi(psi), m_trial(du), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRTangentIntegrator(const LSRTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_psi(other.m_psi), m_trial(other.m_trial),
          m_test(other.m_test), m_state(other.m_state),
          m_params(other.m_params), m_polytope(other.m_polytope),
          m_matrix()
      {}

      LSRTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const std::size_t d = polytope.getDimension();
        const Index cellIdx = polytope.getIndex();
        const auto& trialFES = m_trial.get().getFiniteElementSpace();
        const auto& testFES  = m_test.get().getFiniteElementSpace();
        const auto& trialFE  = trialFES.getFiniteElement(d, cellIdx);
        const auto& testFE   = testFES.getFiniteElement(d, cellIdx);
        const std::size_t trialDofs = trialFE.getCount();
        const std::size_t testDofs  = testFE.getCount();
        const std::size_t vdim = testFES.getVectorDimension();

        m_matrix.resize(testDofs, trialDofs);
        m_matrix.setZero();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(
              lsrQuadOrderFor(
                  std::max(testFE.getOrder(), trialFE.getOrder()), m_params),
              polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const std::size_t nqp = quadrature.getSize();

        for (std::size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const Real wq = qf.getWeight(q);
          const Real distortion = pt.getDistortion();
          const auto& rc = pt.getReferenceCoordinates();

          const Real s = m_psi.get().getValue(ip);
          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> displacement(d);
          displacement.setZero();
          for (std::size_t c = 0; c < vdim; ++c)
            displacement(c) = uqRange(c);

          const auto traced =
            detail::traceMovedPoint(
                pt, displacement, *m_phi, m_params.fieldEvaluation);
          if (traced.exited)
            continue;
          const auto gradPhi = m_grad->getValue(traced.point);

          const Real weight =
            std::exp(-s * s / (2 * m_params.deltaW * m_params.deltaW));
          const Real coef = m_params.rhoS * weight * m_params.normalizer;
          const Real measure = wq * distortion;

          std::vector<Real> gpDotV(testDofs);
          std::vector<Real> gpDotU(trialDofs);
          for (std::size_t te = 0; te < testDofs; ++te)
          {
            const auto basis = testFE.getBasis(te)(rc);
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPhi(c) * basis(c);
            gpDotV[te] = acc;
          }
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
          {
            const auto basis = trialFE.getBasis(tr)(rc);
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPhi(c) * basis(c);
            gpDotU[tr] = acc;
          }
          for (std::size_t te = 0; te < testDofs; ++te)
            for (std::size_t tr = 0; tr < trialDofs; ++tr)
              m_matrix(te, tr) +=
                measure * coef * gpDotU[tr] * gpDotV[te];
        }
        return *this;
      }

      ScalarType integrate(std::size_t tr, std::size_t te) final override
      { return m_matrix(te, tr); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      LSRTangentIntegrator* copy() const noexcept final override
      { return new LSRTangentIntegrator(*this); }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::unique_ptr<GradType> m_grad;
      std::reference_wrapper<const PsiType>  m_psi;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType>  m_test;
      std::reference_wrapper<const StateType> m_state;
      LSRIntegratorParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  // ---- Newton specialisation ------------------------------------------------
  //
  //   K_N = K_GN + sum_q wq * |J_q| * rhoS * W(s_q) * normalizer
  //                * (phi(y_q) - s_q) * v_te^T * hess(y_q) * u_tr.
  template <class PhiDerived, class GradDerived, class HessDerived,
            class PsiType, class TrialType, class TestType, class StateType>
  class LSRTangentIntegrator<
            LSRIntegratorTangentMode::Newton,
            PhiDerived, GradDerived, HessDerived,
            PsiType, TrialType, TestType, StateType>
    final : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using PhiType  = Variational::RealFunctionBase<PhiDerived>;
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;
      using HessType = Variational::MatrixFunctionBase<Real, HessDerived>;
      static constexpr LSRIntegratorTangentMode Mode = LSRIntegratorTangentMode::Newton;

      LSRTangentIntegrator(
          const PhiType& phi, const GradType& grad, const HessType& hess,
          const PsiType& psi,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_phi(phi.copy()), m_grad(grad.copy()), m_hess(hess.copy()),
          m_psi(psi), m_trial(du), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRTangentIntegrator(const LSRTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_hess(other.m_hess ? other.m_hess->copy() : nullptr),
          m_psi(other.m_psi), m_trial(other.m_trial),
          m_test(other.m_test), m_state(other.m_state),
          m_params(other.m_params), m_polytope(other.m_polytope),
          m_matrix()
      {}

      LSRTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const std::size_t d = polytope.getDimension();
        const Index cellIdx = polytope.getIndex();
        const auto& trialFES = m_trial.get().getFiniteElementSpace();
        const auto& testFES  = m_test.get().getFiniteElementSpace();
        const auto& trialFE  = trialFES.getFiniteElement(d, cellIdx);
        const auto& testFE   = testFES.getFiniteElement(d, cellIdx);
        const std::size_t trialDofs = trialFE.getCount();
        const std::size_t testDofs  = testFE.getCount();
        const std::size_t vdim = testFES.getVectorDimension();

        m_matrix.resize(testDofs, trialDofs);
        m_matrix.setZero();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(
              lsrQuadOrderFor(
                  std::max(testFE.getOrder(), trialFE.getOrder()), m_params),
              polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const std::size_t nqp = quadrature.getSize();

        for (std::size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const Real wq = qf.getWeight(q);
          const Real distortion = pt.getDistortion();
          const auto& rc = pt.getReferenceCoordinates();

          const Real s = m_psi.get().getValue(ip);
          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> displacement(d);
          displacement.setZero();
          for (std::size_t c = 0; c < vdim; ++c)
            displacement(c) = uqRange(c);

          const auto traced =
            detail::traceMovedPoint(
                pt, displacement, *m_phi, m_params.fieldEvaluation);
          if (traced.exited)
            continue;

          const Real r =
            m_phi->getValue(traced.point) + traced.correction - s;
          const auto gradPhi = m_grad->getValue(traced.point);
          const auto hessPhi = m_hess->getValue(traced.point);

          const Real weight =
            std::exp(-s * s / (2 * m_params.deltaW * m_params.deltaW));
          const Real coef = m_params.rhoS * weight * m_params.normalizer;
          const Real measure = wq * distortion;

          std::vector<Math::SpatialVector<Real>> testValues(testDofs);
          std::vector<Math::SpatialVector<Real>> trialValues(trialDofs);
          for (std::size_t te = 0; te < testDofs; ++te)
            testValues[te] = testFE.getBasis(te)(rc);
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
            trialValues[tr] = trialFE.getBasis(tr)(rc);

          std::vector<Real> gpDotV(testDofs);
          std::vector<Real> gpDotU(trialDofs);
          for (std::size_t te = 0; te < testDofs; ++te)
          {
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPhi(c) * testValues[te](c);
            gpDotV[te] = acc;
          }
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
          {
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPhi(c) * trialValues[tr](c);
            gpDotU[tr] = acc;
          }
          for (std::size_t te = 0; te < testDofs; ++te)
            for (std::size_t tr = 0; tr < trialDofs; ++tr)
              m_matrix(te, tr) +=
                measure * coef * gpDotU[tr] * gpDotV[te];

          std::vector<Math::SpatialVector<Real>> HU(trialDofs);
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
            HU[tr] = hessPhi * trialValues[tr];
          for (std::size_t te = 0; te < testDofs; ++te)
            for (std::size_t tr = 0; tr < trialDofs; ++tr)
            {
              Real vTHu = 0;
              for (std::size_t c = 0; c < vdim; ++c)
                vTHu += testValues[te](c) * HU[tr](c);
              m_matrix(te, tr) += measure * coef * r * vTHu;
            }
        }
        return *this;
      }

      ScalarType integrate(std::size_t tr, std::size_t te) final override
      { return m_matrix(te, tr); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      LSRTangentIntegrator* copy() const noexcept final override
      { return new LSRTangentIntegrator(*this); }

    private:
      std::unique_ptr<PhiType>  m_phi;
      std::unique_ptr<GradType> m_grad;
      std::unique_ptr<HessType> m_hess;
      std::reference_wrapper<const PsiType>  m_psi;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType>  m_test;
      std::reference_wrapper<const StateType> m_state;
      LSRIntegratorParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  // ---- PSDProjectedNewton specialisation -----------------------------------
  //
  // Identical to the full-Newton specialisation up to the per-quadrature-point
  // assembly of the second-order correction `r * hess(phi)`. Here that 2x2
  // matrix is projected to its PSD part (`(M)+`) before being contracted with
  // the test/trial basis. This produces a safeguarded local model when the
  // residual-weighted Hessian is indefinite. If the projection is inactive on
  // the active quadrature set, the tangent coincides with full Newton. If the
  // projection is active, the method is a modified Newton method and no
  // quadratic residual convergence is implied by the tangent alone.
  //
  template <class PhiDerived, class GradDerived, class HessDerived,
            class PsiType, class TrialType, class TestType, class StateType>
  class LSRTangentIntegrator<
            LSRIntegratorTangentMode::PSDProjectedNewton,
            PhiDerived, GradDerived, HessDerived,
            PsiType, TrialType, TestType, StateType>
    final : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using PhiType  = Variational::RealFunctionBase<PhiDerived>;
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;
      using HessType = Variational::MatrixFunctionBase<Real, HessDerived>;
      static constexpr LSRIntegratorTangentMode Mode = LSRIntegratorTangentMode::PSDProjectedNewton;

      LSRTangentIntegrator(
          const PhiType& phi, const GradType& grad, const HessType& hess,
          const PsiType& psi,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_phi(phi.copy()), m_grad(grad.copy()), m_hess(hess.copy()),
          m_psi(psi), m_trial(du), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRTangentIntegrator(const LSRTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_hess(other.m_hess ? other.m_hess->copy() : nullptr),
          m_psi(other.m_psi), m_trial(other.m_trial),
          m_test(other.m_test), m_state(other.m_state),
          m_params(other.m_params), m_polytope(other.m_polytope),
          m_matrix()
      {}

      LSRTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const std::size_t d = polytope.getDimension();
        const Index cellIdx = polytope.getIndex();
        const auto& trialFES = m_trial.get().getFiniteElementSpace();
        const auto& testFES  = m_test.get().getFiniteElementSpace();
        const auto& trialFE  = trialFES.getFiniteElement(d, cellIdx);
        const auto& testFE   = testFES.getFiniteElement(d, cellIdx);
        const std::size_t trialDofs = trialFE.getCount();
        const std::size_t testDofs  = testFE.getCount();
        const std::size_t vdim = testFES.getVectorDimension();

        m_matrix.resize(testDofs, trialDofs);
        m_matrix.setZero();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(
              lsrQuadOrderFor(
                  std::max(testFE.getOrder(), trialFE.getOrder()), m_params),
              polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const std::size_t nqp = quadrature.getSize();

        for (std::size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const Real wq = qf.getWeight(q);
          const Real distortion = pt.getDistortion();
          const auto& rc = pt.getReferenceCoordinates();

          const Real s = m_psi.get().getValue(ip);
          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> displacement(d);
          displacement.setZero();
          for (std::size_t c = 0; c < vdim; ++c)
            displacement(c) = uqRange(c);

          const auto traced =
            detail::traceMovedPoint(
                pt, displacement, *m_phi, m_params.fieldEvaluation);
          if (traced.exited)
            continue;

          const Real r =
            m_phi->getValue(traced.point) + traced.correction - s;
          const auto gradPhi = m_grad->getValue(traced.point);
          const auto hessPhi = m_hess->getValue(traced.point);

          const Real weight =
            std::exp(-s * s / (2 * m_params.deltaW * m_params.deltaW));
          const Real coef = m_params.rhoS * weight * m_params.normalizer;
          const Real measure = wq * distortion;

          std::vector<Math::SpatialVector<Real>> testValues(testDofs);
          std::vector<Math::SpatialVector<Real>> trialValues(trialDofs);
          for (std::size_t te = 0; te < testDofs; ++te)
            testValues[te] = testFE.getBasis(te)(rc);
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
            trialValues[tr] = trialFE.getBasis(tr)(rc);

          std::vector<Real> gpDotV(testDofs);
          std::vector<Real> gpDotU(trialDofs);
          for (std::size_t te = 0; te < testDofs; ++te)
          {
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPhi(c) * testValues[te](c);
            gpDotV[te] = acc;
          }
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
          {
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPhi(c) * trialValues[tr](c);
            gpDotU[tr] = acc;
          }
          for (std::size_t te = 0; te < testDofs; ++te)
            for (std::size_t tr = 0; tr < trialDofs; ++tr)
              m_matrix(te, tr) +=
                measure * coef * gpDotU[tr] * gpDotV[te];

          // The second-order correction is `r * hess(phi)`. Project to PSD
          // before contracting with the basis — this is the entire point of
          // this specialisation, see the class docstring.
          const Math::SpatialMatrix<Real> rH = r * hessPhi;
          const Math::SpatialMatrix<Real> rHplus = detail::psdProject2x2(rH);

          std::vector<Math::SpatialVector<Real>> HU(trialDofs);
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
            HU[tr] = rHplus * trialValues[tr];
          for (std::size_t te = 0; te < testDofs; ++te)
            for (std::size_t tr = 0; tr < trialDofs; ++tr)
            {
              Real vTHu = 0;
              for (std::size_t c = 0; c < vdim; ++c)
                vTHu += testValues[te](c) * HU[tr](c);
              m_matrix(te, tr) += measure * coef * vTHu;
            }
        }
        return *this;
      }

      ScalarType integrate(std::size_t tr, std::size_t te) final override
      { return m_matrix(te, tr); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      LSRTangentIntegrator* copy() const noexcept final override
      { return new LSRTangentIntegrator(*this); }

    private:
      std::unique_ptr<PhiType>  m_phi;
      std::unique_ptr<GradType> m_grad;
      std::unique_ptr<HessType> m_hess;
      std::reference_wrapper<const PsiType>  m_psi;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType>  m_test;
      std::reference_wrapper<const StateType> m_state;
      LSRIntegratorParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };


  // ---------------------------------------------------------------------------
  // E_pull data term.
  //
  // Energy:
  //   E_pull(u) = (rhoS/2) * integral_Omega
  //                 w_delta(psi(X)) * normalizer
  //                 * ( phi(X) - psi(X - u(X)) )^2 dX.
  //
  // Notes.
  //   - phi is evaluated at X (no displacement).
  //   - psi is evaluated at the pulled-back point X - u(X). The integrator
  //     constructs a Geometry::Point whose physical coordinates are
  //     X - u(X) and whose polytope is the ORIGINAL cell. The user-supplied
  //     psi-displaced / grad-psi-displaced adapters decide how to interpret
  //     that point. The analytic case simply reads the physical
  //     coordinates; for a discrete psi the caller wraps it through
  //     Variational::Flow with a constant velocity equal to -u(X) at the
  //     starting integration point (Phi_1(X) = X - u(X) exactly).
  //   - The band weight uses psi at X. A separate `psiBand` input is taken;
  //     it can be the raw GridFunction in the discrete case.
  //
  // Residual:
  //   delta r / delta u = grad psi(X - u(X)).
  //   R[i] = sum_q w_q * |J_q| * rhoS * w_delta(psi(X)) * normalizer
  //          * ( phi(X) - psi(X - u(X)) )
  //          * ( grad psi(X - u(X)) . b_i(X) ).
  //
  // Gauss-Newton tangent (per qpt rank-1; PSD by construction):
  //   K[i, j] = sum_q w_q * |J_q| * rhoS * w_delta(psi(X)) * normalizer
  //             * ( grad psi(X - u(X)) . b_i ) ( grad psi(X - u(X)) . b_j ).
  // ---------------------------------------------------------------------------
  namespace detail
  {
    inline Geometry::Point makeDisplacedPoint(
        const Geometry::Point& original,
        const Math::SpatialVector<Real>& yPhysical)
    {
      return Geometry::Point(
          original.getPolytope(),
          original.getReferenceCoordinates(),
          yPhysical);
    }
  }

  template <class PhiDerived, class PsiBandType,
            class PsiDispDerived, class GradPsiDispDerived,
            class TestType, class StateType>
  class LSRPullResidualIntegrator final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType       = Real;
      using PhiType          = Variational::RealFunctionBase<PhiDerived>;
      using PsiDispType      = Variational::RealFunctionBase<PsiDispDerived>;
      using GradPsiDispType  = Variational::VectorFunctionBase<Real, GradPsiDispDerived>;

      LSRPullResidualIntegrator(
          const PhiType& phi,
          const PsiBandType& psiBand,
          const PsiDispType& psiDisp,
          const GradPsiDispType& gradPsiDisp,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LinearFormIntegratorBase<Real>(v),
          m_phi(phi.copy()),
          m_psiBand(psiBand),
          m_psiDisp(psiDisp.copy()),
          m_gradPsiDisp(gradPsiDisp.copy()),
          m_test(v), m_state(u),
          m_params(params)
      {}

      LSRPullResidualIntegrator(const LSRPullResidualIntegrator& other)
        : Variational::LinearFormIntegratorBase<Real>(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_psiBand(other.m_psiBand),
          m_psiDisp(other.m_psiDisp ? other.m_psiDisp->copy() : nullptr),
          m_gradPsiDisp(other.m_gradPsiDisp ? other.m_gradPsiDisp->copy() : nullptr),
          m_test(other.m_test), m_state(other.m_state),
          m_params(other.m_params),
          m_polytope(other.m_polytope),
          m_elemVec()
      {}

      LSRPullResidualIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const std::size_t d = polytope.getDimension();
        const Index cellIdx = polytope.getIndex();
        const auto& testFES = m_test.get().getFiniteElementSpace();
        const auto& testFE = testFES.getFiniteElement(d, cellIdx);
        const std::size_t testDofs = testFE.getCount();
        const std::size_t vdim = testFES.getVectorDimension();

        m_elemVec.resize(testDofs);
        m_elemVec.setZero();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(
              lsrQuadOrderFor(testFE.getOrder(), m_params),
              polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const std::size_t nqp = quadrature.getSize();

        for (std::size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const Real wq = qf.getWeight(q);
          const Real distortion = pt.getDistortion();
          const auto& rc = pt.getReferenceCoordinates();

          // Band weight uses psi at X.
          const Real s = m_psiBand.get().getValue(ip);

          // Pulled-back evaluation point: X - u(X).
          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> yMinus(d);
          for (std::size_t c = 0; c < vdim; ++c)
            yMinus(c) = pt.getCoordinates()(c) - uqRange(c);
          const Geometry::Point pulledPoint =
            detail::makeDisplacedPoint(pt, yMinus);

          const Real phi_X = m_phi->getValue(ip);
          const Real psiDisp = m_psiDisp->getValue(pulledPoint);
          const auto gradPsiDisp = m_gradPsiDisp->getValue(pulledPoint);

          const Real r = phi_X - psiDisp;
          const Real weight =
            std::exp(-s * s / (2 * m_params.deltaW * m_params.deltaW));
          const Real coef =
            m_params.rhoS * weight * r * m_params.normalizer;
          const Real measure = wq * distortion;

          for (std::size_t te = 0; te < testDofs; ++te)
          {
            const auto testValue = testFE.getBasis(te)(rc);
            Real dot = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              dot += gradPsiDisp(c) * testValue(c);
            m_elemVec(te) += measure * coef * dot;
          }
        }
        return *this;
      }

      ScalarType integrate(std::size_t te) final override
      { return m_elemVec(te); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      LSRPullResidualIntegrator* copy() const noexcept final override
      { return new LSRPullResidualIntegrator(*this); }

    private:
      std::unique_ptr<PhiType>         m_phi;
      std::reference_wrapper<const PsiBandType> m_psiBand;
      std::unique_ptr<PsiDispType>     m_psiDisp;
      std::unique_ptr<GradPsiDispType> m_gradPsiDisp;
      std::reference_wrapper<const TestType>  m_test;
      std::reference_wrapper<const StateType> m_state;
      LSRIntegratorParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  template <class P, class B, class PD, class GP, class V, class U>
  LSRPullResidualIntegrator(
      const Variational::RealFunctionBase<P>&,
      const B&,
      const Variational::RealFunctionBase<PD>&,
      const Variational::VectorFunctionBase<Real, GP>&,
      const V&, const U&, LSRIntegratorParameters)
    -> LSRPullResidualIntegrator<P, B, PD, GP, V, U>;

  template <LSRIntegratorTangentMode Mode,
            class PhiDerived, class PsiBandType,
            class PsiDispDerived, class GradPsiDispDerived,
            class HessPsiDispDerived,
            class TrialType, class TestType, class StateType>
  class LSRPullTangentIntegrator;

  // ---- E_pull Gauss-Newton tangent ------------------------------------------
  template <class PhiDerived, class PsiBandType,
            class PsiDispDerived, class GradPsiDispDerived,
            class HessPsiDispDerived,
            class TrialType, class TestType, class StateType>
  class LSRPullTangentIntegrator<
            LSRIntegratorTangentMode::GaussNewton,
            PhiDerived, PsiBandType,
            PsiDispDerived, GradPsiDispDerived, HessPsiDispDerived,
            TrialType, TestType, StateType>
    final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using GradPsiDispType =
        Variational::VectorFunctionBase<Real, GradPsiDispDerived>;
      static constexpr LSRIntegratorTangentMode Mode =
        LSRIntegratorTangentMode::GaussNewton;

      LSRPullTangentIntegrator(
          const PsiBandType& psiBand,
          const GradPsiDispType& gradPsiDisp,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_psiBand(psiBand),
          m_gradPsiDisp(gradPsiDisp.copy()),
          m_trial(du), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRPullTangentIntegrator(const LSRPullTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_psiBand(other.m_psiBand),
          m_gradPsiDisp(other.m_gradPsiDisp ? other.m_gradPsiDisp->copy() : nullptr),
          m_trial(other.m_trial), m_test(other.m_test),
          m_state(other.m_state),
          m_params(other.m_params),
          m_polytope(other.m_polytope),
          m_matrix()
      {}

      LSRPullTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const std::size_t d = polytope.getDimension();
        const Index cellIdx = polytope.getIndex();
        const auto& trialFES = m_trial.get().getFiniteElementSpace();
        const auto& testFES  = m_test.get().getFiniteElementSpace();
        const auto& trialFE  = trialFES.getFiniteElement(d, cellIdx);
        const auto& testFE   = testFES.getFiniteElement(d, cellIdx);
        const std::size_t trialDofs = trialFE.getCount();
        const std::size_t testDofs  = testFE.getCount();
        const std::size_t vdim = testFES.getVectorDimension();

        m_matrix.resize(testDofs, trialDofs);
        m_matrix.setZero();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(
              lsrQuadOrderFor(
                  std::max(testFE.getOrder(), trialFE.getOrder()), m_params),
              polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const std::size_t nqp = quadrature.getSize();

        for (std::size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const Real wq = qf.getWeight(q);
          const Real distortion = pt.getDistortion();
          const auto& rc = pt.getReferenceCoordinates();

          const Real s = m_psiBand.get().getValue(ip);

          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> yMinus(d);
          for (std::size_t c = 0; c < vdim; ++c)
            yMinus(c) = pt.getCoordinates()(c) - uqRange(c);
          const Geometry::Point pulledPoint =
            detail::makeDisplacedPoint(pt, yMinus);

          const auto gradPsiDisp = m_gradPsiDisp->getValue(pulledPoint);

          const Real weight =
            std::exp(-s * s / (2 * m_params.deltaW * m_params.deltaW));
          const Real coef = m_params.rhoS * weight * m_params.normalizer;
          const Real measure = wq * distortion;

          std::vector<Real> gpDotV(testDofs);
          std::vector<Real> gpDotU(trialDofs);
          for (std::size_t te = 0; te < testDofs; ++te)
          {
            const auto basis = testFE.getBasis(te)(rc);
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPsiDisp(c) * basis(c);
            gpDotV[te] = acc;
          }
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
          {
            const auto basis = trialFE.getBasis(tr)(rc);
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPsiDisp(c) * basis(c);
            gpDotU[tr] = acc;
          }
          for (std::size_t te = 0; te < testDofs; ++te)
            for (std::size_t tr = 0; tr < trialDofs; ++tr)
              m_matrix(te, tr) +=
                measure * coef * gpDotU[tr] * gpDotV[te];
        }
        return *this;
      }

      ScalarType integrate(std::size_t tr, std::size_t te) final override
      { return m_matrix(te, tr); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      LSRPullTangentIntegrator* copy() const noexcept final override
      { return new LSRPullTangentIntegrator(*this); }

    private:
      std::reference_wrapper<const PsiBandType> m_psiBand;
      std::unique_ptr<GradPsiDispType> m_gradPsiDisp;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType>  m_test;
      std::reference_wrapper<const StateType> m_state;
      LSRIntegratorParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  template <class B, class GP, class DU, class V, class U>
  LSRPullTangentIntegrator(
      const B&,
      const Variational::VectorFunctionBase<Real, GP>&,
      const DU&, const V&, const U&, LSRIntegratorParameters)
    -> LSRPullTangentIntegrator<
         LSRIntegratorTangentMode::GaussNewton,
         void, B, void, GP, void, DU, V, U>;

  // ---- E_pull full Newton tangent -------------------------------------------
  //
  // With r(u) = phi(X) - psi(X - u(X)), one has
  //   D r[u][w] = grad psi(X-u) . w
  // and
  //   D^2 r[u][z,w] = - z^T hess psi(X-u) w.
  // Hence the exact data Hessian is
  //   grad psi \otimes grad psi - r hess psi.
  template <class PhiDerived, class PsiBandType,
            class PsiDispDerived, class GradPsiDispDerived,
            class HessPsiDispDerived,
            class TrialType, class TestType, class StateType>
  class LSRPullTangentIntegrator<
            LSRIntegratorTangentMode::Newton,
            PhiDerived, PsiBandType,
            PsiDispDerived, GradPsiDispDerived, HessPsiDispDerived,
            TrialType, TestType, StateType>
    final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
      using PsiDispType = Variational::RealFunctionBase<PsiDispDerived>;
      using GradPsiDispType =
        Variational::VectorFunctionBase<Real, GradPsiDispDerived>;
      using HessPsiDispType =
        Variational::MatrixFunctionBase<Real, HessPsiDispDerived>;
      static constexpr LSRIntegratorTangentMode Mode =
        LSRIntegratorTangentMode::Newton;

      LSRPullTangentIntegrator(
          const PhiType& phi,
          const PsiBandType& psiBand,
          const PsiDispType& psiDisp,
          const GradPsiDispType& gradPsiDisp,
          const HessPsiDispType& hessPsiDisp,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_phi(phi.copy()),
          m_psiBand(psiBand),
          m_psiDisp(psiDisp.copy()),
          m_gradPsiDisp(gradPsiDisp.copy()),
          m_hessPsiDisp(hessPsiDisp.copy()),
          m_trial(du), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRPullTangentIntegrator(const LSRPullTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_psiBand(other.m_psiBand),
          m_psiDisp(other.m_psiDisp ? other.m_psiDisp->copy() : nullptr),
          m_gradPsiDisp(other.m_gradPsiDisp ? other.m_gradPsiDisp->copy() : nullptr),
          m_hessPsiDisp(other.m_hessPsiDisp ? other.m_hessPsiDisp->copy() : nullptr),
          m_trial(other.m_trial), m_test(other.m_test),
          m_state(other.m_state),
          m_params(other.m_params),
          m_polytope(other.m_polytope),
          m_matrix()
      {}

      LSRPullTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const std::size_t d = polytope.getDimension();
        const Index cellIdx = polytope.getIndex();
        const auto& trialFES = m_trial.get().getFiniteElementSpace();
        const auto& testFES  = m_test.get().getFiniteElementSpace();
        const auto& trialFE  = trialFES.getFiniteElement(d, cellIdx);
        const auto& testFE   = testFES.getFiniteElement(d, cellIdx);
        const std::size_t trialDofs = trialFE.getCount();
        const std::size_t testDofs  = testFE.getCount();
        const std::size_t vdim = testFES.getVectorDimension();

        m_matrix.resize(testDofs, trialDofs);
        m_matrix.setZero();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(
              lsrQuadOrderFor(
                  std::max(testFE.getOrder(), trialFE.getOrder()), m_params),
              polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);

        for (std::size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const Real wq = qf.getWeight(q);
          const Real distortion = pt.getDistortion();
          const auto& rc = pt.getReferenceCoordinates();

          const Real s = m_psiBand.get().getValue(ip);

          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> yMinus(d);
          for (std::size_t c = 0; c < vdim; ++c)
            yMinus(c) = pt.getCoordinates()(c) - uqRange(c);
          const Geometry::Point pulledPoint =
            detail::makeDisplacedPoint(pt, yMinus);

          const Real phiX = m_phi->getValue(ip);
          const Real psiDisp = m_psiDisp->getValue(pulledPoint);
          const Real r = phiX - psiDisp;
          const auto gradPsiDisp = m_gradPsiDisp->getValue(pulledPoint);
          const auto hessPsiDisp = m_hessPsiDisp->getValue(pulledPoint);

          const Real weight =
            std::exp(-s * s / (2 * m_params.deltaW * m_params.deltaW));
          const Real coef = m_params.rhoS * weight * m_params.normalizer;
          const Real measure = wq * distortion;

          std::vector<Math::SpatialVector<Real>> testValues(testDofs);
          std::vector<Math::SpatialVector<Real>> trialValues(trialDofs);
          for (std::size_t te = 0; te < testDofs; ++te)
            testValues[te] = testFE.getBasis(te)(rc);
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
            trialValues[tr] = trialFE.getBasis(tr)(rc);

          std::vector<Real> gpDotV(testDofs);
          std::vector<Real> gpDotU(trialDofs);
          for (std::size_t te = 0; te < testDofs; ++te)
          {
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPsiDisp(c) * testValues[te](c);
            gpDotV[te] = acc;
          }
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
          {
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPsiDisp(c) * trialValues[tr](c);
            gpDotU[tr] = acc;
          }
          for (std::size_t te = 0; te < testDofs; ++te)
            for (std::size_t tr = 0; tr < trialDofs; ++tr)
              m_matrix(te, tr) +=
                measure * coef * gpDotU[tr] * gpDotV[te];

          std::vector<Math::SpatialVector<Real>> HU(trialDofs);
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
            HU[tr] = hessPsiDisp * trialValues[tr];
          for (std::size_t te = 0; te < testDofs; ++te)
            for (std::size_t tr = 0; tr < trialDofs; ++tr)
            {
              Real vTHu = 0;
              for (std::size_t c = 0; c < vdim; ++c)
                vTHu += testValues[te](c) * HU[tr](c);
              m_matrix(te, tr) -= measure * coef * r * vTHu;
            }
        }
        return *this;
      }

      ScalarType integrate(std::size_t tr, std::size_t te) final override
      { return m_matrix(te, tr); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      LSRPullTangentIntegrator* copy() const noexcept final override
      { return new LSRPullTangentIntegrator(*this); }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::reference_wrapper<const PsiBandType> m_psiBand;
      std::unique_ptr<PsiDispType> m_psiDisp;
      std::unique_ptr<GradPsiDispType> m_gradPsiDisp;
      std::unique_ptr<HessPsiDispType> m_hessPsiDisp;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType>  m_test;
      std::reference_wrapper<const StateType> m_state;
      LSRIntegratorParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  // ---- E_pull PSD-projected Newton tangent ----------------------------------
  template <class PhiDerived, class PsiBandType,
            class PsiDispDerived, class GradPsiDispDerived,
            class HessPsiDispDerived,
            class TrialType, class TestType, class StateType>
  class LSRPullTangentIntegrator<
            LSRIntegratorTangentMode::PSDProjectedNewton,
            PhiDerived, PsiBandType,
            PsiDispDerived, GradPsiDispDerived, HessPsiDispDerived,
            TrialType, TestType, StateType>
    final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
      using PsiDispType = Variational::RealFunctionBase<PsiDispDerived>;
      using GradPsiDispType =
        Variational::VectorFunctionBase<Real, GradPsiDispDerived>;
      using HessPsiDispType =
        Variational::MatrixFunctionBase<Real, HessPsiDispDerived>;
      static constexpr LSRIntegratorTangentMode Mode =
        LSRIntegratorTangentMode::PSDProjectedNewton;

      LSRPullTangentIntegrator(
          const PhiType& phi,
          const PsiBandType& psiBand,
          const PsiDispType& psiDisp,
          const GradPsiDispType& gradPsiDisp,
          const HessPsiDispType& hessPsiDisp,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_phi(phi.copy()),
          m_psiBand(psiBand),
          m_psiDisp(psiDisp.copy()),
          m_gradPsiDisp(gradPsiDisp.copy()),
          m_hessPsiDisp(hessPsiDisp.copy()),
          m_trial(du), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRPullTangentIntegrator(const LSRPullTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_psiBand(other.m_psiBand),
          m_psiDisp(other.m_psiDisp ? other.m_psiDisp->copy() : nullptr),
          m_gradPsiDisp(other.m_gradPsiDisp ? other.m_gradPsiDisp->copy() : nullptr),
          m_hessPsiDisp(other.m_hessPsiDisp ? other.m_hessPsiDisp->copy() : nullptr),
          m_trial(other.m_trial), m_test(other.m_test),
          m_state(other.m_state),
          m_params(other.m_params),
          m_polytope(other.m_polytope),
          m_matrix()
      {}

      LSRPullTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const std::size_t d = polytope.getDimension();
        const Index cellIdx = polytope.getIndex();
        const auto& trialFES = m_trial.get().getFiniteElementSpace();
        const auto& testFES  = m_test.get().getFiniteElementSpace();
        const auto& trialFE  = trialFES.getFiniteElement(d, cellIdx);
        const auto& testFE   = testFES.getFiniteElement(d, cellIdx);
        const std::size_t trialDofs = trialFE.getCount();
        const std::size_t testDofs  = testFE.getCount();
        const std::size_t vdim = testFES.getVectorDimension();

        m_matrix.resize(testDofs, trialDofs);
        m_matrix.setZero();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(
              lsrQuadOrderFor(
                  std::max(testFE.getOrder(), trialFE.getOrder()), m_params),
              polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);

        for (std::size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const Real wq = qf.getWeight(q);
          const Real distortion = pt.getDistortion();
          const auto& rc = pt.getReferenceCoordinates();

          const Real s = m_psiBand.get().getValue(ip);

          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> yMinus(d);
          for (std::size_t c = 0; c < vdim; ++c)
            yMinus(c) = pt.getCoordinates()(c) - uqRange(c);
          const Geometry::Point pulledPoint =
            detail::makeDisplacedPoint(pt, yMinus);

          const Real phiX = m_phi->getValue(ip);
          const Real psiDisp = m_psiDisp->getValue(pulledPoint);
          const Real r = phiX - psiDisp;
          const auto gradPsiDisp = m_gradPsiDisp->getValue(pulledPoint);
          const auto hessPsiDisp = m_hessPsiDisp->getValue(pulledPoint);

          const Real weight =
            std::exp(-s * s / (2 * m_params.deltaW * m_params.deltaW));
          const Real coef = m_params.rhoS * weight * m_params.normalizer;
          const Real measure = wq * distortion;

          std::vector<Math::SpatialVector<Real>> testValues(testDofs);
          std::vector<Math::SpatialVector<Real>> trialValues(trialDofs);
          for (std::size_t te = 0; te < testDofs; ++te)
            testValues[te] = testFE.getBasis(te)(rc);
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
            trialValues[tr] = trialFE.getBasis(tr)(rc);

          std::vector<Real> gpDotV(testDofs);
          std::vector<Real> gpDotU(trialDofs);
          for (std::size_t te = 0; te < testDofs; ++te)
          {
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPsiDisp(c) * testValues[te](c);
            gpDotV[te] = acc;
          }
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
          {
            Real acc = 0;
            for (std::size_t c = 0; c < vdim; ++c)
              acc += gradPsiDisp(c) * trialValues[tr](c);
            gpDotU[tr] = acc;
          }
          for (std::size_t te = 0; te < testDofs; ++te)
            for (std::size_t tr = 0; tr < trialDofs; ++tr)
              m_matrix(te, tr) +=
                measure * coef * gpDotU[tr] * gpDotV[te];

          const Math::SpatialMatrix<Real> second = -r * hessPsiDisp;
          const Math::SpatialMatrix<Real> secondPlus =
            detail::psdProject2x2(second);

          std::vector<Math::SpatialVector<Real>> HU(trialDofs);
          for (std::size_t tr = 0; tr < trialDofs; ++tr)
            HU[tr] = secondPlus * trialValues[tr];
          for (std::size_t te = 0; te < testDofs; ++te)
            for (std::size_t tr = 0; tr < trialDofs; ++tr)
            {
              Real vTHu = 0;
              for (std::size_t c = 0; c < vdim; ++c)
                vTHu += testValues[te](c) * HU[tr](c);
              m_matrix(te, tr) += measure * coef * vTHu;
            }
        }
        return *this;
      }

      ScalarType integrate(std::size_t tr, std::size_t te) final override
      { return m_matrix(te, tr); }

      const Geometry::Polytope& getPolytope() const final override
      { return m_polytope->get(); }

      Geometry::Region getRegion() const final override
      { return Geometry::Region::Cells; }

      LSRPullTangentIntegrator* copy() const noexcept final override
      { return new LSRPullTangentIntegrator(*this); }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::reference_wrapper<const PsiBandType> m_psiBand;
      std::unique_ptr<PsiDispType> m_psiDisp;
      std::unique_ptr<GradPsiDispType> m_gradPsiDisp;
      std::unique_ptr<HessPsiDispType> m_hessPsiDisp;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType>  m_test;
      std::reference_wrapper<const StateType> m_state;
      LSRIntegratorParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

}

#endif
