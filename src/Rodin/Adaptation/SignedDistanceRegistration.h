/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_SIGNEDDISTANCEREGISTRATION_H
#define RODIN_ADAPTATION_SIGNEDDISTANCEREGISTRATION_H

/**
 * @file
 * @brief High-level helper for the SDR penalty, modelled on
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
 *     params SDRParameters (defaulted)
 *
 * Usage:
 *
 *     SignedDistanceRegistration sdr(du, v, u);
 *
 *     newton = sdr(phi, grad,        sLF, params);   // GaussNewton
 *     newton = sdr(phi, grad, hess,  sLF, params);   // Newton
 *
 * The user supplies each derivative as a Rodin FunctionBase. The
 * integrators evaluate them at the deformed point y = X + u_h(X) using
 * `FunctionBase::getValue(Geometry::Point)`. They do NOT call
 * `Variational::Grad` or `Variational::Hess` internally — the
 * derivatives are the user's inputs.
 *
 * Overload resolution on the call signature picks the right
 * `SDRTangentIntegrator` specialisation at compile time. There is no
 * runtime tangent-mode flag.
 */

#include <functional>

#include "Rodin/Types.h"

#include "SDRIntegrators.h"

namespace Rodin::Adaptation
{
  template <class Trial, class Test, class State>
  class SignedDistanceRegistration
  {
    public:
      SignedDistanceRegistration(
          const Trial& du, const Test& v, const State& u)
        : m_du(du), m_v(v), m_u(u)
      {}

      // ---- Decomposed: Tangent ----------------------------------------------
      template <class PhiDerived, class GradDerived, class SLF>
      auto Tangent(
          const Variational::RealFunctionBase<PhiDerived>& /*phi*/,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const SLF& sLF,
          SDRParameters params = {}) const
      {
        return SDRTangentIntegrator<
                  SDRTangentMode::GaussNewton,
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
          SDRParameters params = {}) const
      {
        return SDRTangentIntegrator<
                  SDRTangentMode::Newton,
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
          SDRParameters params = {}) const
      {
        return SDRResidualIntegrator<
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
          SDRParameters params = {}) const
      {
        return Residual(phi, grad, sLF, params);
      }

      // ---- Composite --------------------------------------------------------
      template <class PhiDerived, class GradDerived, class SLF>
      auto operator()(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const SLF& sLF,
          SDRParameters params = {}) const
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
          SDRParameters params = {}) const
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
  SignedDistanceRegistration(const Trial&, const Test&, const State&)
    -> SignedDistanceRegistration<Trial, Test, State>;
}

#endif
