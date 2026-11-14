/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Perfusion.h
 * @brief Pore storage, effective pore pressure and shear-dependent drainage
 * of the reduced poroelastic shell.
 *
 * Exponential vascular compliance @f$ d\Phi/dp_P = C_0 e^{-p_P/p_0} @f$:
 * @f[
 *   p_P(\bar\Phi) = -p_0 \ln\Big(1 - \frac{\bar\Phi - \phi_0}{\Delta\Phi_{max}}\Big),
 *   \qquad \Delta\Phi_{max} = C_0 p_0,
 * @f]
 * effective pore pressure @f$ \tilde p = p_P(\bar\Phi) - \bar\lambda @f$, supply
 * @f$ Q_{ar} = \gamma_{ar}(p_{ar} - \tilde p) + Q_{ar}^{ext} @f$ and drainage
 * @f[
 *   Q_{ven} = \frac{G_{ven,0}}{\nu_{eff}(\dot\gamma_{eff})} (\tilde p - p_{ven}),
 *   \qquad
 *   \dot\gamma_{eff} = \alpha_{pm}
 *     \sqrt{\frac{\bar\lambda_r}{\bar\lambda_\theta^2 K \bar\Phi}}
 *     \frac{|Q_{ar}|}{A_{w0}},
 * @f]
 * with Carreau-Yasuda @f$ \nu_{eff} @f$. The shear rate is regularized as
 * @f$ \sqrt{\dot\gamma_{eff}^2 + \varepsilon_\gamma^2} @f$.
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_PHYSICS_PERFUSION_H
#define RODIN_HEART_POROELASTICSPHERE_PHYSICS_PERFUSION_H

#include <cassert>
#include <cmath>

namespace Rodin::Heart::PoroelasticSphere::Physics
{
  /**
   * @brief Carreau-Yasuda dynamic viscosity and its shear-rate derivative.
   * @tparam Params Parameter type with eta0, etaInf, lambda, n, a.
   */
  template <class Params, class Scalar>
  inline void carreauYasuda(const Params& cy, Scalar g, Scalar& nu, Scalar& dnu)
  {
    const Scalar lg = cy.lambda * g;
    const Scalar la = std::pow(lg, cy.a);
    const Scalar base = Scalar(1) + la;
    const Scalar ex = (cy.n - Scalar(1)) / cy.a;
    nu = cy.etaInf + (cy.eta0 - cy.etaInf) * std::pow(base, ex);
    dnu = (g > Scalar(0))
      ? (cy.eta0 - cy.etaInf) * (cy.n - Scalar(1)) * std::pow(base, ex - Scalar(1)) * la / g
      : Scalar(0);
  }

  /**
   * @brief Evaluates the storage law, the pore pressure and the supply and
   * drainage rates with their derivatives with respect to
   * @f$ (y, \bar\Phi, p_v, p_{ar}, e_c) @f$.
   *
   * @tparam Input Model input parameter type.
   */
  template <class Input>
  class PerfusionEvaluator
  {
    public:
      /// @brief Construct with model input parameters.
      explicit PerfusionEvaluator(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Storage pressure @f$ p_P @f$ and @f$ p_P' @f$.
       *
       * Beyond @f$ z_c @f$ (fraction of @f$ \Delta\Phi_{max} @f$) the law is
       * continued linearly so that Newton trials outside the constitutive
       * domain stay finite; converged states inside it are unaffected.
       */
      template <class Scalar>
      void storage(Scalar phi, Scalar& pP, Scalar& dpP) const
      {
        const Scalar p0 = m_input.storagePressureScale;
        const Scalar dPhiMax = m_input.storageCompliance * p0;
        assert(dPhiMax > Scalar(0));
        const Scalar zc = m_input.storageSaturationGuard;
        const Scalar z = (phi - m_input.phi0) / dPhiMax;
        if (z < zc)
        {
          pP = -p0 * std::log(Scalar(1) - z);
          dpP = p0 / (dPhiMax * (Scalar(1) - z));
        }
        else
        {
          dpP = p0 / (dPhiMax * (Scalar(1) - zc));
          pP = -p0 * std::log(Scalar(1) - zc) + dpP * dPhiMax * (z - zc);
        }
      }

      /**
       * @brief Evaluate @f$ \tilde p @f$, @f$ Q_{ar} @f$, @f$ Q_{ven} @f$.
       *
       * Requires kinematics, @f$ \bar\lambda @f$ and the external supply.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);

        Scalar pP = 0, dpP = 0;
        storage(data.phi, pP, dpP);
        data.pP = pP;

        data.pf = pP - data.lambdaBar;
        data.dPf_dy = -data.dLambdaBar_dy;
        data.dPf_dphi = dpP - data.dLambdaBar_dphi;
        data.dPf_dpv = -data.dLambdaBar_dpv;
        data.dPf_dec = -data.dLambdaBar_dec;

        // Supply.
        const Scalar gAr = m_input.gammaAr;
        data.perfusionLumped = gAr * (data.par - data.pf);
        data.qArterial = data.perfusionLumped + data.externalPerfusionInflow;
        const Scalar Qa = data.qArterial;
        const Scalar Qa_y = -gAr * data.dPf_dy;
        const Scalar Qa_phi = -gAr * data.dPf_dphi;
        const Scalar Qa_pv = -gAr * data.dPf_dpv;
        const Scalar Qa_par = gAr;
        const Scalar Qa_ec = -gAr * data.dPf_dec;

        // Effective pore shear rate.
        const Scalar K = m_input.permeability;
        assert(K > Scalar(0) && data.phi > Scalar(0));
        const Scalar S2 = data.lr / (data.lt * data.lt * K * data.phi);
        const Scalar S2_y = S2 * (data.lr_y / data.lr - Scalar(2) * data.lt_y / data.lt);
        const Scalar S2_phi = S2 * (data.lr_phi / data.lr - Scalar(1) / data.phi);
        const Scalar c = m_input.poreShearFactor * m_input.poreShearFactor /
          (data.Aw0 * data.Aw0);
        const Scalar eps = m_input.poreViscosity.shearRegularization;
        const Scalar g = std::sqrt(c * S2 * Qa * Qa + eps * eps);
        data.shearRate = g;
        const Scalar k = Scalar(0.5) * c / g;
        const Scalar g_y = k * (S2_y * Qa * Qa + Scalar(2) * S2 * Qa * Qa_y);
        const Scalar g_phi = k * (S2_phi * Qa * Qa + Scalar(2) * S2 * Qa * Qa_phi);
        const Scalar g_pv = k * Scalar(2) * S2 * Qa * Qa_pv;
        const Scalar g_par = k * Scalar(2) * S2 * Qa * Qa_par;
        const Scalar g_ec = k * Scalar(2) * S2 * Qa * Qa_ec;

        Scalar nu = 0, nu_g = 0;
        carreauYasuda(m_input.poreViscosity, g, nu, nu_g);
        data.poreViscosity = nu;

        // Drainage.
        const Scalar G = m_input.venousConductance;
        const Scalar dp = data.pf - data.pSvMid;
        data.qVenous = G / nu * dp;
        const Scalar a = -G * nu_g * dp / (nu * nu);
        const Scalar b = G / nu;

        // Net source Q = Q_ar - Q_ven and its derivatives.
        data.netPerfusion = data.qArterial - data.qVenous;
        data.dNet_dy = Qa_y - (b * data.dPf_dy + a * g_y);
        data.dNet_dphi = Qa_phi - (b * data.dPf_dphi + a * g_phi);
        data.dNet_dpv = Qa_pv - (b * data.dPf_dpv + a * g_pv);
        data.dNet_dpar = Qa_par - a * g_par;
        data.dNet_dec = Qa_ec - (b * data.dPf_dec + a * g_ec);
      }

    private:
      const Input& m_input;
  };
}

#endif
