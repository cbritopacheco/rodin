/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_LSRREGISTRATION_H
#define RODIN_ADAPTATION_LSRREGISTRATION_H

/**
 * @file
 * @brief High-level helper for the LSR penalty, modelled on
 *        `Rodin::Variational::LinearElasticityIntegral`.
 *
 * The composite captures the VARIATIONAL SKELETON `(du, v, u)` at
 * construction. The DATA — three Rodin form-language objects describing
 * the level-set target grid function, plus parameters
 * — is passed at `operator()` / `Tangent` / `Residual` time:
 *
 *     phi    Variational::RealFunctionBase<...>           (always)
 *     grad   Variational::VectorFunctionBase<Real, ...>   (always)
 *     hess   Variational::MatrixFunctionBase<Real, ...>   (Newton only)
 *     psi    GridFunction (level-set data target)
 *     params LSRIntegratorParameters (defaulted)
 *
 * Usage:
 *
 *     LSRRegistration lsr(du, v, u);
 *
 *     newton = lsr(phi, grad,        psi, params);  // GaussNewton
 *     newton = lsr(phi, grad, hess,  psi, params);  // Newton
 *
 * The user supplies each derivative as a Rodin FunctionBase. The
 * integrators evaluate them at the deformed point y = X + u_h(X) using
 * `FunctionBase::getValue(Geometry::Point)`. They do NOT call
 * `Variational::Grad` or `Variational::Hess` internally — the
 * derivatives are the user's inputs.
 *
 * Overload resolution on the call signature picks the right
 * `LSRTangentIntegrator` specialisation at compile time. There is no
 * runtime tangent-mode flag.
 */

#include <functional>

#include "Rodin/Types.h"

#include "LSRIntegrators.h"

namespace Rodin::Adaptation
{
  template <class Trial, class Test, class State>
  class LSRRegistration
  {
    public:
      LSRRegistration(
          const Trial& du, const Test& v, const State& u)
        : m_du(du), m_v(v), m_u(u)
      {}

      // ---- Decomposed: Tangent ----------------------------------------------
      template <class PhiDerived, class GradDerived, class Psi>
      auto Tangent(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Psi& psi,
          LSRIntegratorParameters params = {}) const
      {
        return LSRTangentIntegrator<
                  LSRIntegratorTangentMode::GaussNewton,
                  PhiDerived, GradDerived, void,
                  Psi, Trial, Test, State>(
              phi, grad, psi,
              m_du.get(), m_v.get(), m_u.get(),
              params);
      }

      template <class PhiDerived, class GradDerived, class HessDerived,
                class Psi>
      auto Tangent(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Variational::MatrixFunctionBase<Real, HessDerived>& hess,
          const Psi& psi,
          LSRIntegratorParameters params = {}) const
      {
        return LSRTangentIntegrator<
                  LSRIntegratorTangentMode::Newton,
                  PhiDerived, GradDerived, HessDerived,
                  Psi, Trial, Test, State>(
              phi, grad, hess, psi,
              m_du.get(), m_v.get(), m_u.get(),
              params);
      }

      /**
       * @brief PSD-projected full-Newton tangent.
       *
       * Same call signature as `Tangent(phi, grad, hess, psi, params)`
       * (CTAD on the level-set Hessian), but the per-qpt second-order
       * correction `r * hess(phi)` is clamped to its PSD part before
       * being added to the local block. The global tangent stays SPD
       * even when `r` changes sign across the band, so Newton contracts
       * past the GN noise floor without the indefiniteness that
       * destabilises raw full-Newton on this objective.
       */
      template <class PhiDerived, class GradDerived, class HessDerived,
                class Psi>
      auto TangentPSDProjected(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Variational::MatrixFunctionBase<Real, HessDerived>& hess,
          const Psi& psi,
          LSRIntegratorParameters params = {}) const
      {
        return LSRTangentIntegrator<
                  LSRIntegratorTangentMode::PSDProjectedNewton,
                  PhiDerived, GradDerived, HessDerived,
                  Psi, Trial, Test, State>(
              phi, grad, hess, psi,
              m_du.get(), m_v.get(), m_u.get(),
              params);
      }

      // ---- Decomposed: Residual ---------------------------------------------
      template <class PhiDerived, class GradDerived, class Psi>
      auto Residual(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Psi& psi,
          LSRIntegratorParameters params = {}) const
      {
        return LSRResidualIntegrator<
                  PhiDerived, GradDerived, Psi, Test, State>(
              phi, grad, psi, m_v.get(), m_u.get(), params);
      }

      // Hess accepted for call-site symmetry with Tangent; ignored here.
      template <class PhiDerived, class GradDerived, class HessDerived,
                class Psi>
      auto Residual(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Variational::MatrixFunctionBase<Real, HessDerived>& /*hess*/,
          const Psi& psi,
          LSRIntegratorParameters params = {}) const
      {
        return Residual(phi, grad, psi, params);
      }

      template <class PhiDerived, class GradDerived>
      auto InterfaceResidual(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          LSRIntegratorParameters params = {}) const
      {
        return LSRInterfaceDistanceResidualIntegrator<
                  PhiDerived, GradDerived, Test, State>(
              phi, grad, m_v.get(), m_u.get(), params);
      }

      template <class PhiDerived, class GradDerived>
      auto InterfaceTangent(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          LSRIntegratorParameters params = {}) const
      {
        return LSRInterfaceDistanceTangentIntegrator<
                  PhiDerived, GradDerived, Trial, Test, State>(
              phi, grad, m_du.get(), m_v.get(), m_u.get(), params);
      }

      // ---- Composite --------------------------------------------------------
      template <class PhiDerived, class GradDerived, class Psi>
      auto operator()(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Psi& psi,
          LSRIntegratorParameters params = {}) const
      {
        return Tangent(phi, grad, psi, params)
             + Residual(phi, grad, psi, params);
      }

      template <class PhiDerived, class GradDerived, class HessDerived,
                class Psi>
      auto operator()(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad,
          const Variational::MatrixFunctionBase<Real, HessDerived>& hess,
          const Psi& psi,
          LSRIntegratorParameters params = {}) const
      {
        return Tangent(phi, grad, hess, psi, params)
             + Residual(phi, grad, psi, params);
      }

      // ---- E_pull data term --------------------------------------------------
      //
      // Energy form:
      //   E_pull(u) = (rhoS/2) integral W(psi(X)) ( phi(X) - psi(X - u(X)) )^2 dX.
      //
      // Inputs.
      //   phi          : RealFunctionBase, evaluated at X. (Analytic or any
      //                  function callable at original integration points.)
      //   psiBand      : the scalar field used by the band weight w(psi(X)).
      //                  Typically the raw psi GridFunction.
      //   psiDisp      : RealFunctionBase adapter returning psi at X - u(X).
      //                  For analytic psi this is the analytic closure
      //                  reading physical coords. For a discrete psi this is
      //                  a Variational::Flow wrapper with constant velocity
      //                  -u(X) at the starting integration point.
      //   gradPsiDisp  : VectorFunctionBase adapter returning grad psi at
      //                  X - u(X). Same wiring as psiDisp.
      //
      template <class PhiDerived, class PsiBand,
                class PsiDispDerived, class GradPsiDispDerived>
      auto PullResidual(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const PsiBand& psiBand,
          const Variational::RealFunctionBase<PsiDispDerived>& psiDisp,
          const Variational::VectorFunctionBase<Real, GradPsiDispDerived>& gradPsiDisp,
          LSRIntegratorParameters params = {}) const
      {
        return LSRPullResidualIntegrator<
                  PhiDerived, PsiBand,
                  PsiDispDerived, GradPsiDispDerived,
                  Test, State>(
              phi, psiBand, psiDisp, gradPsiDisp,
              m_v.get(), m_u.get(), params);
      }

      template <class PsiBand, class GradPsiDispDerived>
      auto PullTangent(
          const PsiBand& psiBand,
          const Variational::VectorFunctionBase<Real, GradPsiDispDerived>& gradPsiDisp,
          LSRIntegratorParameters params = {}) const
      {
        return LSRPullTangentIntegrator<
                  LSRIntegratorTangentMode::GaussNewton,
                  void, PsiBand,
                  void, GradPsiDispDerived, void,
                  Trial, Test, State>(
              psiBand, gradPsiDisp,
              m_du.get(), m_v.get(), m_u.get(), params);
      }

      template <class PhiDerived, class PsiBand,
                class PsiDispDerived, class GradPsiDispDerived,
                class HessPsiDispDerived>
      auto PullTangent(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const PsiBand& psiBand,
          const Variational::RealFunctionBase<PsiDispDerived>& psiDisp,
          const Variational::VectorFunctionBase<Real, GradPsiDispDerived>& gradPsiDisp,
          const Variational::MatrixFunctionBase<Real, HessPsiDispDerived>& hessPsiDisp,
          LSRIntegratorParameters params = {}) const
      {
        return LSRPullTangentIntegrator<
                  LSRIntegratorTangentMode::Newton,
                  PhiDerived, PsiBand,
                  PsiDispDerived, GradPsiDispDerived, HessPsiDispDerived,
                  Trial, Test, State>(
              phi, psiBand, psiDisp, gradPsiDisp, hessPsiDisp,
              m_du.get(), m_v.get(), m_u.get(), params);
      }

      template <class PhiDerived, class PsiBand,
                class PsiDispDerived, class GradPsiDispDerived,
                class HessPsiDispDerived>
      auto PullTangentPSDProjected(
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const PsiBand& psiBand,
          const Variational::RealFunctionBase<PsiDispDerived>& psiDisp,
          const Variational::VectorFunctionBase<Real, GradPsiDispDerived>& gradPsiDisp,
          const Variational::MatrixFunctionBase<Real, HessPsiDispDerived>& hessPsiDisp,
          LSRIntegratorParameters params = {}) const
      {
        return LSRPullTangentIntegrator<
                  LSRIntegratorTangentMode::PSDProjectedNewton,
                  PhiDerived, PsiBand,
                  PsiDispDerived, GradPsiDispDerived, HessPsiDispDerived,
                  Trial, Test, State>(
              phi, psiBand, psiDisp, gradPsiDisp, hessPsiDisp,
              m_du.get(), m_v.get(), m_u.get(), params);
      }

    private:
      std::reference_wrapper<const Trial> m_du;
      std::reference_wrapper<const Test>  m_v;
      std::reference_wrapper<const State> m_u;
  };

  template <class Trial, class Test, class State>
  LSRRegistration(const Trial&, const Test&, const State&)
    -> LSRRegistration<Trial, Test, State>;

}

#endif
