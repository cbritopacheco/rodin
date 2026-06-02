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
 *   coordinates are overridden with `X + u_h(X)` and whose polytope and
 *   reference coordinates are copied from the original quadrature
 *   point, then call `phi->getValue(movedPoint)`,
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
                        ///< The global tangent stays PSD, so Newton
                        ///< contracts super-linearly past the GN noise
                        ///< floor without the tangential indefiniteness
                        ///< that destabilises raw full-Newton wherever
                        ///< the LSR residual r changes sign.
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
    Real rhoS = 1;
    Real deltaW = 0;
    Real hRef = 0;
    Real normalizer = 1;
  };

  inline std::size_t lsrQuadOrderFor(std::size_t feOrder)
  {
    return std::max<std::size_t>(2, 2 * feOrder);
  }

  namespace detail
  {
    /// Build a Geometry::Point at the deformed location y = X + u(X),
    /// reusing the parent polytope and reference coordinates of the
    /// original quadrature point and overriding the physical
    /// coordinates with `y`.
    inline Geometry::Point makeMovedPoint(
        const Geometry::Point& original,
        const Math::SpatialVector<Real>& yPhysical)
    {
      return Geometry::Point(
          original.getPolytope(),
          original.getReferenceCoordinates(),
          yPhysical);
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
            class SLFType, class TestType, class StateType>
  class LSRResidualIntegrator final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using PhiType  = Variational::RealFunctionBase<PhiDerived>;
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      LSRResidualIntegrator(
          const PhiType& phi, const GradType& grad,
          const SLFType& sLF,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LinearFormIntegratorBase<Real>(v),
          m_phi(phi.copy()), m_grad(grad.copy()),
          m_sLF(sLF), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRResidualIntegrator(const LSRResidualIntegrator& other)
        : Variational::LinearFormIntegratorBase<Real>(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_sLF(other.m_sLF), m_test(other.m_test), m_state(other.m_state),
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
              lsrQuadOrderFor(testFE.getOrder()),
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

          const Real s = m_sLF.get().getValue(ip);
          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> y(d);
          for (std::size_t c = 0; c < vdim; ++c)
            y(c) = pt.getCoordinates()(c) + uqRange(c);

          const Geometry::Point movedPoint = detail::makeMovedPoint(pt, y);

          const Real phi_y = m_phi->getValue(movedPoint);
          const auto gradPhi = m_grad->getValue(movedPoint);

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
      std::reference_wrapper<const SLFType> m_sLF;
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
            class SLFType, class TrialType, class TestType, class StateType>
  class LSRTangentIntegrator;

  // ---- GaussNewton specialisation ------------------------------------------
  //
  //   K_GN[te, tr] = sum_q wq * |J_q| * rhoS * W(s_q) * normalizer
  //                  * (grad(y_q) . v_te) * (grad(y_q) . u_tr).
  template <class PhiDerived, class GradDerived, class HessDerived,
            class SLFType, class TrialType, class TestType, class StateType>
  class LSRTangentIntegrator<
            LSRIntegratorTangentMode::GaussNewton,
            PhiDerived, GradDerived, HessDerived,
            SLFType, TrialType, TestType, StateType>
    final : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;
      static constexpr LSRIntegratorTangentMode Mode = LSRIntegratorTangentMode::GaussNewton;

      LSRTangentIntegrator(
          const GradType& grad,
          const SLFType& sLF,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_grad(grad.copy()),
          m_sLF(sLF), m_trial(du), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRTangentIntegrator(const LSRTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_sLF(other.m_sLF), m_trial(other.m_trial),
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
              lsrQuadOrderFor(std::max(testFE.getOrder(), trialFE.getOrder())),
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

          const Real s = m_sLF.get().getValue(ip);
          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> y(d);
          for (std::size_t c = 0; c < vdim; ++c)
            y(c) = pt.getCoordinates()(c) + uqRange(c);

          const Geometry::Point movedPoint = detail::makeMovedPoint(pt, y);
          const auto gradPhi = m_grad->getValue(movedPoint);

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
      std::unique_ptr<GradType> m_grad;
      std::reference_wrapper<const SLFType>  m_sLF;
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
            class SLFType, class TrialType, class TestType, class StateType>
  class LSRTangentIntegrator<
            LSRIntegratorTangentMode::Newton,
            PhiDerived, GradDerived, HessDerived,
            SLFType, TrialType, TestType, StateType>
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
          const SLFType& sLF,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_phi(phi.copy()), m_grad(grad.copy()), m_hess(hess.copy()),
          m_sLF(sLF), m_trial(du), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRTangentIntegrator(const LSRTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_hess(other.m_hess ? other.m_hess->copy() : nullptr),
          m_sLF(other.m_sLF), m_trial(other.m_trial),
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
              lsrQuadOrderFor(std::max(testFE.getOrder(), trialFE.getOrder())),
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

          const Real s = m_sLF.get().getValue(ip);
          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> y(d);
          for (std::size_t c = 0; c < vdim; ++c)
            y(c) = pt.getCoordinates()(c) + uqRange(c);

          const Geometry::Point movedPoint = detail::makeMovedPoint(pt, y);

          const Real r = m_phi->getValue(movedPoint) - s;
          const auto gradPhi = m_grad->getValue(movedPoint);
          const auto hessPhi = m_hess->getValue(movedPoint);

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
      std::reference_wrapper<const SLFType>  m_sLF;
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
  // matrix is projected to its PSD part (`(M)+`) BEFORE being contracted with
  // the test/trial basis. The result: the global tangent stays SPD, so Newton
  // contracts even when r changes sign across the band — which is exactly the
  // regime where raw full-Newton breaks (the tangential eigenvalue r/||p-c||
  // becomes negative and the indefinite step flips cells through the singular
  // floor).
  //
  // Where r * hess(phi) is already PSD this mode is identical to full Newton
  // and recovers quadratic convergence. Where it is mixed-sign, this mode
  // keeps only the PSD-contributing eigenspace and degenerates gracefully
  // toward GN on the indefinite directions — but does NOT cycle around the
  // GN minimum, because the PSD piece is supplied at every iteration and
  // contracts the dropped-term spectral radius.
  //
  template <class PhiDerived, class GradDerived, class HessDerived,
            class SLFType, class TrialType, class TestType, class StateType>
  class LSRTangentIntegrator<
            LSRIntegratorTangentMode::PSDProjectedNewton,
            PhiDerived, GradDerived, HessDerived,
            SLFType, TrialType, TestType, StateType>
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
          const SLFType& sLF,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          LSRIntegratorParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_phi(phi.copy()), m_grad(grad.copy()), m_hess(hess.copy()),
          m_sLF(sLF), m_trial(du), m_test(v), m_state(u),
          m_params(params)
      {}

      LSRTangentIntegrator(const LSRTangentIntegrator& other)
        : Variational::LocalBilinearFormIntegratorBase<Real>(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_hess(other.m_hess ? other.m_hess->copy() : nullptr),
          m_sLF(other.m_sLF), m_trial(other.m_trial),
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
              lsrQuadOrderFor(std::max(testFE.getOrder(), trialFE.getOrder())),
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

          const Real s = m_sLF.get().getValue(ip);
          const auto uqRange = m_state.get().getValue(ip);
          Math::SpatialVector<Real> y(d);
          for (std::size_t c = 0; c < vdim; ++c)
            y(c) = pt.getCoordinates()(c) + uqRange(c);

          const Geometry::Point movedPoint = detail::makeMovedPoint(pt, y);

          const Real r = m_phi->getValue(movedPoint) - s;
          const auto gradPhi = m_grad->getValue(movedPoint);
          const auto hessPhi = m_hess->getValue(movedPoint);

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
      std::reference_wrapper<const SLFType>  m_sLF;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType>  m_test;
      std::reference_wrapper<const StateType> m_state;
      LSRIntegratorParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

}

#endif
