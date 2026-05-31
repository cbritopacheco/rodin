/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_SDFRREGISTRATION_H
#define RODIN_ADAPTATION_SDFRREGISTRATION_H

/**
 * @file
 * @brief High-level helper for the SDFR penalty, modelled on
 *        `Rodin::Variational::LinearElasticityIntegral`.
 *
 * The composite captures the VARIATIONAL SKELETON `(du, v, u)` at
 * construction. The DATA — three Rodin form-language objects describing
 * the level set, plus the signed-distance grid function, plus parameters
 * — is passed at `operator()` / `Tangent` / `Residual` time:
 *
 *     phi    Variational::RealFunctionBase<...>           (always)
 *     grad   Variational::VectorFunctionBase<Real, ...>   (always)
 *     hess   Variational::MatrixFunctionBase<Real, ...>   (Newton only)
 *     sLF    GridFunction (signed distance / data target)
 *     params SDFRIntegratorParameters (defaulted)
 *
 * Usage:
 *
 *     SDFRRegistration sdfr(du, v, u);
 *
 *     newton = sdfr(phi, grad,        sLF, params);  // GaussNewton
 *     newton = sdfr(phi, grad, hess,  sLF, params);  // Newton
 *
 * The user supplies each derivative as a Rodin FunctionBase. The
 * integrators evaluate them at the deformed point y = X + u_h(X) using
 * `FunctionBase::getValue(Geometry::Point)`. They do NOT call
 * `Variational::Grad` or `Variational::Hess` internally — the
 * derivatives are the user's inputs.
 *
 * Overload resolution on the call signature picks the right
 * `SDFRTangentIntegrator` specialisation at compile time. There is no
 * runtime tangent-mode flag.
 */

#include <functional>

#include "Rodin/Types.h"

#include "SDFRIntegrators.h"

namespace Rodin::Adaptation
{
  template <class Trial, class Test, class State>
  class SDFRRegistration
  {
    public:
      SDFRRegistration(
          const Trial& du, const Test& v, const State& u)
        : m_du(du), m_v(v), m_u(u)
      {}

      // ---- Decomposed: Tangent ----------------------------------------------
      template <class PhiDerived, class GradDerived, class SLF>
      auto Tangent(
          const Variational::RealFunctionBase<PhiDerived>& /*phi*/,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const SLF& sLF,
          SDFRIntegratorParameters params = {}) const
      {
        return SDFRTangentIntegrator<
                  SDFRIntegratorTangentMode::GaussNewton,
                  PhiDerived, GradDerived, void,
                  SLF, Trial, Test, State>(
              grad, sLF,
              m_du.get(), m_v.get(), m_u.get(),
              params);
      }

      template <class PhiDerived, class GradDerived, class HessDerived,
                class SLF>
      auto Tangent(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Variational::MatrixFunctionBase<Real, HessDerived>& hess,
          const SLF& sLF,
          SDFRIntegratorParameters params = {}) const
      {
        return SDFRTangentIntegrator<
                  SDFRIntegratorTangentMode::Newton,
                  PhiDerived, GradDerived, HessDerived,
                  SLF, Trial, Test, State>(
              phi, grad, hess, sLF,
              m_du.get(), m_v.get(), m_u.get(),
              params);
      }

      /**
       * @brief PSD-projected full-Newton tangent.
       *
       * Same call signature as `Tangent(phi, grad, hess, sLF, params)`
       * (CTAD on the level-set Hessian), but the per-qpt second-order
       * correction `r * hess(phi)` is clamped to its PSD part before
       * being added to the local block. The global tangent stays SPD
       * even when `r` changes sign across the band, so Newton contracts
       * past the GN noise floor without the indefiniteness that
       * destabilises raw full-Newton on this objective.
       */
      template <class PhiDerived, class GradDerived, class HessDerived,
                class SLF>
      auto TangentPSDProjected(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Variational::MatrixFunctionBase<Real, HessDerived>& hess,
          const SLF& sLF,
          SDFRIntegratorParameters params = {}) const
      {
        return SDFRTangentIntegrator<
                  SDFRIntegratorTangentMode::PSDProjectedNewton,
                  PhiDerived, GradDerived, HessDerived,
                  SLF, Trial, Test, State>(
              phi, grad, hess, sLF,
              m_du.get(), m_v.get(), m_u.get(),
              params);
      }

      // ---- Decomposed: Residual ---------------------------------------------
      template <class PhiDerived, class GradDerived, class SLF>
      auto Residual(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const SLF& sLF,
          SDFRIntegratorParameters params = {}) const
      {
        return SDFRResidualIntegrator<
                  PhiDerived, GradDerived, SLF, Test, State>(
              phi, grad, sLF, m_v.get(), m_u.get(), params);
      }

      // Hess accepted for call-site symmetry with Tangent; ignored here.
      template <class PhiDerived, class GradDerived, class HessDerived,
                class SLF>
      auto Residual(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Variational::MatrixFunctionBase<Real, HessDerived>& /*hess*/,
          const SLF& sLF,
          SDFRIntegratorParameters params = {}) const
      {
        return Residual(phi, grad, sLF, params);
      }

      // ---- Composite --------------------------------------------------------
      template <class PhiDerived, class GradDerived, class SLF>
      auto operator()(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const SLF& sLF,
          SDFRIntegratorParameters params = {}) const
      {
        return Tangent(phi, grad, sLF, params)
             + Residual(phi, grad, sLF, params);
      }

      template <class PhiDerived, class GradDerived, class HessDerived,
                class SLF>
      auto operator()(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Variational::MatrixFunctionBase<Real, HessDerived>& hess,
          const SLF& sLF,
          SDFRIntegratorParameters params = {}) const
      {
        return Tangent(phi, grad, hess, sLF, params)
             + Residual(phi, grad, sLF, params);
      }

    private:
      std::reference_wrapper<const Trial> m_du;
      std::reference_wrapper<const Test>  m_v;
      std::reference_wrapper<const State> m_u;
  };

  template <class Trial, class Test, class State>
  SDFRRegistration(const Trial&, const Test&, const State&)
    -> SDFRRegistration<Trial, Test, State>;

}

#endif
