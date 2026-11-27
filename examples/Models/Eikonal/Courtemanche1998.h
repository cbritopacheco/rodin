/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Adaptation and stepping code distributed under the Boost Software License,
 * Version 1.0 (see LICENSE). The underlying CellML model equations retain
 * the attribution and CC BY 3.0 terms stated below.
 */
/*
 * Ionic equations adapted from Courtemanche, Ramirez and Nattel (1998),
 * Physiome Model Repository CellML revision 78335420620450b2eab817577c0faf4450bbf8cb.
 * CellML curation: Penny Noble (Oxford University); COR export: Alan Garny.
 * Model material: CC BY 3.0, https://creativecommons.org/licenses/by/3.0/.
 * Source and modifications are documented in ReactionEikonal.md.
 */
#ifndef RODIN_EIKONAL_COURTEMANCHE1998_H
#define RODIN_EIKONAL_COURTEMANCHE1998_H
#include <array>
#include <cmath>
#include "Rodin/Types.h"
namespace Rodin::Examples
{
  /**
   * @brief Baseline human atrial ionic kinetics in mV, ms, mM and pA/pF.
   *
   * The external stimulus is positive for depolarization. Gates use a frozen-
   * coefficient Rush-Larsen step; voltage and concentrations use forward Euler.
   * All right-hand sides are evaluated at the supplied old state.
   */
  class Courtemanche1998
  {
    public:
      using State = std::array<Real, 21>;

      /// @brief State indices inherited from the published model equations.
      enum Component
      {
        V,
        m,
        h,
        j,
        oa,
        oi,
        ua,
        ui,
        xr,
        xs,
        d,
        f,
        f_Ca,
        u,
        v,
        w,
        Na_i,
        K_i,
        Ca_i,
        Ca_up,
        Ca_rel
      };

      /// @returns Published resting initial state, before optional pacing.
      static State initialState()
      {
        return {-81.18, 2.908e-3, 9.649e-1, 9.775e-1, 3.043e-2, 9.992e-1, 4.966e-3,
          9.986e-1, 3.296e-5, 1.869e-2, 1.367e-4, 9.996e-1, 7.755e-1, 2.35e-112, 1,
          0.9992, 1.117e+01, 1.39e+02, 1.013e-4, 1.488, 1.488};
      }

      /**
       * @brief Advances one local ionic step without diffusion.
       * @param state Old state.
       * @param dt Step in milliseconds.
       * @param stimulus Depolarizing current in pA/pF.
       * @returns Updated state.
       */
      static State step(const State& state, Real dt, Real stimulus)
      {
        const Real i_st = -100 * stimulus;
        const Real V = state[0];
        const Real m = state[1];
        const Real h = state[2];
        const Real j = state[3];
        const Real oa = state[4];
        const Real oi = state[5];
        const Real ua = state[6];
        const Real ui = state[7];
        const Real xr = state[8];
        const Real xs = state[9];
        const Real d = state[10];
        const Real f = state[11];
        const Real f_Ca = state[12];
        const Real u = state[13];
        const Real v = state[14];
        const Real w = state[15];
        const Real Na_i = state[16];
        const Real K_i = state[17];
        const Real Ca_i = state[18];
        const Real Ca_up = state[19];
        const Real Ca_rel = state[20];
        const Real R = 8.3143;
        const Real T = 310;
        const Real F = 96.4867;
        const Real Cm = 100;
        const Real g_Na = 7.8;
        const Real g_K1 = 0.09;
        const Real K_Q10 = 3;
        const Real g_to = 0.1652;
        const Real g_Kr = 0.029411765;
        const Real g_Ks = 0.12941176;
        const Real g_Ca_L = 0.12375;
        const Real Km_Na_i = 10;
        const Real Km_K_o = 1.5;
        const Real i_NaK_max = 0.59933874;
        const Real g_B_Na = 0.0006744375;
        const Real g_B_Ca = 0.001131;
        const Real g_B_K = 0;
        const Real I_NaCa_max = 1600;
        const Real K_mNa = 87.5;
        const Real K_mCa = 1.38;
        const Real K_sat = 0.1;
        const Real gamma = 0.35;
        const Real i_CaP_max = 0.275;
        const Real K_rel = 30;
        const Real tau_tr = 180;
        const Real I_up_max = 0.005;
        const Real K_up = 0.00092;
        const Real Ca_up_max = 15;
        const Real CMDN_max = 0.05;
        const Real TRPN_max = 0.07;
        const Real CSQN_max = 10;
        const Real Km_CMDN = 0.00238;
        const Real Km_TRPN = 0.0005;
        const Real Km_CSQN = 0.8;
        const Real V_cell = 20100;
        const Real Na_o = 140;
        const Real Ca_o = 1.8;
        const Real K_o = 5.4;
        const Real E_Na = (((R * T) / F) * std::log((Na_o / Na_i)));
        const Real alpha_m = ((V == (-47.13))
            ? 3.2
            : ((0.32 * (V + 47.13)) / (1.0 - std::exp(((-0.1) * (V + 47.13))))));
        const Real beta_m = (0.08 * std::exp(((-V) / 11.0)));
        const Real m_inf = (alpha_m / (alpha_m + beta_m));
        const Real tau_m = (1.0 / (alpha_m + beta_m));
        const Real alpha_h =
          ((V < (-40.0)) ? (0.135 * std::exp(((V + 80.0) / (-6.8)))) : 0.0);
        const Real beta_h = ((V < (-40.0))
            ? ((3.56 * std::exp((0.079 * V))) + (3.1e5 * std::exp((0.35 * V))))
            : (1.0 / (0.13 * (1.0 + std::exp(((V + 10.66) / (-11.1)))))));
        const Real h_inf = (alpha_h / (alpha_h + beta_h));
        const Real tau_h = (1.0 / (alpha_h + beta_h));
        const Real alpha_j =
          ((V < (-40.0)) ? (((((-1.2714e5) * std::exp((0.2444 * V))) -
                               (3.474e-5 * std::exp(((-0.04391) * V)))) *
                              (V + 37.78)) /
                             (1.0 + std::exp((0.311 * (V + 79.23)))))
                         : 0.0);
        const Real beta_j =
          ((V < (-40.0)) ? ((0.1212 * std::exp(((-0.01052) * V))) /
                             (1.0 + std::exp(((-0.1378) * (V + 40.14)))))
                         : ((0.3 * std::exp(((-2.535e-7) * V))) /
                             (1.0 + std::exp(((-0.1) * (V + 32.0))))));
        const Real j_inf = (alpha_j / (alpha_j + beta_j));
        const Real tau_j = (1.0 / (alpha_j + beta_j));
        const Real E_K = (((R * T) / F) * std::log((K_o / K_i)));
        const Real i_K1 =
          ((Cm * g_K1 * (V - E_K)) / (1.0 + std::exp((0.07 * (V + 80.0)))));
        const Real i_to = (Cm * g_to * std::pow(oa, 3.0) * oi * (V - E_K));
        const Real alpha_oa = (0.65 *
          std::pow((std::exp(((V - (-10.0)) / (-8.5))) +
                     std::exp((((V - (-10.0)) - 40.0) / (-59.0)))),
            (-1.0)));
        const Real beta_oa =
          (0.65 * std::pow((2.5 + std::exp((((V - (-10.0)) + 72.0) / 17.0))), (-1.0)));
        const Real tau_oa = (std::pow((alpha_oa + beta_oa), (-1.0)) / K_Q10);
        const Real oa_infinity =
          std::pow((1.0 + std::exp((((V - (-10.0)) + 10.47) / (-17.54)))), (-1.0));
        const Real alpha_oi =
          std::pow((18.53 + (1.0 * std::exp((((V - (-10.0)) + 103.7) / 10.95)))), (-1.0));
        const Real beta_oi = std::pow(
          (35.56 + (1.0 * std::exp((((V - (-10.0)) - 8.74) / (-7.44))))), (-1.0));
        const Real tau_oi = (std::pow((alpha_oi + beta_oi), (-1.0)) / K_Q10);
        const Real oi_infinity =
          std::pow((1.0 + std::exp((((V - (-10.0)) + 33.1) / 5.3))), (-1.0));
        const Real g_Kur = (0.005 + (0.05 / (1.0 + std::exp(((V - 15.0) / (-13.0))))));
        const Real i_Kur = (Cm * g_Kur * std::pow(ua, 3.0) * ui * (V - E_K));
        const Real alpha_ua = (0.65 *
          std::pow((std::exp(((V - (-10.0)) / (-8.5))) +
                     std::exp((((V - (-10.0)) - 40.0) / (-59.0)))),
            (-1.0)));
        const Real beta_ua =
          (0.65 * std::pow((2.5 + std::exp((((V - (-10.0)) + 72.0) / 17.0))), (-1.0)));
        const Real tau_ua = (std::pow((alpha_ua + beta_ua), (-1.0)) / K_Q10);
        const Real ua_infinity =
          std::pow((1.0 + std::exp((((V - (-10.0)) + 20.3) / (-9.6)))), (-1.0));
        const Real alpha_ui = std::pow(
          (21.0 + (1.0 * std::exp((((V - (-10.0)) - 195.0) / (-28.0))))), (-1.0));
        const Real beta_ui = (1.0 / std::exp((((V - (-10.0)) - 168.0) / (-16.0))));
        const Real tau_ui = (std::pow((alpha_ui + beta_ui), (-1.0)) / K_Q10);
        const Real ui_infinity =
          std::pow((1.0 + std::exp((((V - (-10.0)) - 109.45) / 27.48))), (-1.0));
        const Real i_Kr =
          ((Cm * g_Kr * xr * (V - E_K)) / (1.0 + std::exp(((V + 15.0) / 22.4))));
        const Real alpha_xr = ((std::abs((V + 14.1)) < 1e-10)
            ? 0.0015
            : ((0.0003 * (V + 14.1)) / (1.0 - std::exp(((V + 14.1) / (-5.0))))));
        const Real beta_xr = ((std::abs((V - 3.3328)) < 1e-10)
            ? 3.7836118e-4
            : ((7.3898e-05 * (V - 3.3328)) / (std::exp(((V - 3.3328) / 5.1237)) - 1.0)));
        const Real tau_xr = std::pow((alpha_xr + beta_xr), (-1.0));
        const Real xr_infinity =
          std::pow((1.0 + std::exp(((V + 14.1) / (-6.5)))), (-1.0));
        const Real i_Ks = (Cm * g_Ks * std::pow(xs, 2.0) * (V - E_K));
        const Real alpha_xs = ((std::abs((V - 19.9)) < 1e-10)
            ? 0.00068
            : ((4e-05 * (V - 19.9)) / (1.0 - std::exp(((V - 19.9) / (-17.0))))));
        const Real beta_xs = ((std::abs((V - 19.9)) < 1e-10)
            ? 0.000315
            : ((3.5e-05 * (V - 19.9)) / (std::exp(((V - 19.9) / 9.0)) - 1.0)));
        const Real tau_xs = (0.5 * std::pow((alpha_xs + beta_xs), (-1.0)));
        const Real xs_infinity =
          std::pow((1.0 + std::exp(((V - 19.9) / (-12.7)))), (-0.5));
        const Real i_Ca_L = (Cm * g_Ca_L * d * f * f_Ca * (V - 65.0));
        const Real d_infinity = std::pow((1.0 + std::exp(((V + 10.0) / (-8.0)))), (-1.0));
        const Real tau_d = ((std::abs((V + 10.0)) < 1e-10)
            ? (4.579 / (1.0 + std::exp(((V + 10.0) / (-6.24)))))
            : ((1.0 - std::exp(((V + 10.0) / (-6.24)))) /
                (0.035 * (V + 10.0) * (1.0 + std::exp(((V + 10.0) / (-6.24)))))));
        const Real f_infinity =
          (std::exp(((-(V + 28.0)) / 6.9)) / (1.0 + std::exp(((-(V + 28.0)) / 6.9))));
        const Real tau_f = (9.0 *
          std::pow(
            ((0.0197 * std::exp(((-std::pow(0.0337, 2.0)) * std::pow((V + 10.0), 2.0)))) +
              0.02),
            (-1.0)));
        const Real f_Ca_infinity = std::pow((1.0 + (Ca_i / 0.00035)), (-1.0));
        const Real tau_f_Ca = 2.0;
        const Real sigma = ((1.0 / 7.0) * (std::exp((Na_o / 67.3)) - 1.0));
        const Real f_NaK =
          std::pow((1.0 + (0.1245 * std::exp((((-0.1) * F * V) / (R * T)))) +
                     (0.0365 * sigma * std::exp((((-F) * V) / (R * T))))),
            (-1.0));
        const Real i_NaK =
          ((((Cm * i_NaK_max * f_NaK * 1.0) / (1.0 + std::pow((Km_Na_i / Na_i), 1.5))) *
             K_o) /
            (K_o + Km_K_o));
        const Real E_Ca = (((R * T) / (2.0 * F)) * std::log((Ca_o / Ca_i)));
        const Real i_B_Na = (Cm * g_B_Na * (V - E_Na));
        const Real i_B_Ca = (Cm * g_B_Ca * (V - E_Ca));
        const Real i_B_K = (Cm * g_B_K * (V - E_K));
        const Real i_NaCa =
          ((Cm * I_NaCa_max *
             ((std::exp(((gamma * F * V) / (R * T))) * std::pow(Na_i, 3.0) * Ca_o) -
               (std::exp((((gamma - 1.0) * F * V) / (R * T))) * std::pow(Na_o, 3.0) *
                 Ca_i))) /
            ((std::pow(K_mNa, 3.0) + std::pow(Na_o, 3.0)) * (K_mCa + Ca_o) *
              (1.0 + (K_sat * std::exp((((gamma - 1.0) * V * F) / (R * T)))))));
        const Real i_CaP = ((Cm * i_CaP_max * Ca_i) / (0.0005 + Ca_i));
        const Real i_rel = (K_rel * std::pow(u, 2.0) * v * w * (Ca_rel - Ca_i));
        const Real tau_u = 8.0;
        const Real tau_w = ((std::abs((V - 7.9)) < 1e-10)
            ? ((6.0 * 0.2) / 1.3)
            : ((6.0 * (1.0 - std::exp(((-(V - 7.9)) / 5.0)))) /
                ((1.0 + (0.3 * std::exp(((-(V - 7.9)) / 5.0)))) * 1.0 * (V - 7.9))));
        const Real w_infinity =
          (1.0 - std::pow((1.0 + std::exp(((-(V - 40.0)) / 17.0))), (-1.0)));
        const Real i_tr = ((Ca_up - Ca_rel) / tau_tr);
        const Real i_up = (I_up_max / (1.0 + (K_up / Ca_i)));
        const Real i_up_leak = ((I_up_max * Ca_up) / Ca_up_max);
        const Real V_i = (V_cell * 0.68);
        const Real V_rel = (0.0048 * V_cell);
        const Real V_up = (0.0552 * V_cell);
        const Real B1 =
          ((((2.0 * i_NaCa) - (i_CaP + i_Ca_L + i_B_Ca)) / (2.0 * V_i * F)) +
            (((V_up * (i_up_leak - i_up)) + (i_rel * V_rel)) / V_i));
        const Real B2 = (1.0 + ((TRPN_max * Km_TRPN) / std::pow((Ca_i + Km_TRPN), 2.0)) +
          ((CMDN_max * Km_CMDN) / std::pow((Ca_i + Km_CMDN), 2.0)));
        const Real i_Na = (Cm * g_Na * std::pow(m, 3.0) * h * j * (V - E_Na));
        const Real Fn = (1e3 *
          ((1e-15 * V_rel * i_rel) -
            ((1e-15 / (2.0 * F)) * ((0.5 * i_Ca_L) - (0.2 * i_NaCa)))));
        const Real u_infinity =
          std::pow((1.0 + std::exp(((-(Fn - 3.4175e-13)) / 13.67e-16))), (-1.0));
        const Real tau_v = (1.91 +
          (2.09 *
            std::pow((1.0 + std::exp(((-(Fn - 3.4175e-13)) / 13.67e-16))), (-1.0))));
        const Real v_infinity =
          (1.0 - std::pow((1.0 + std::exp(((-(Fn - 6.835e-14)) / 13.67e-16))), (-1.0)));
        State next{};
        next[0] = V +
          dt *
            ((-(i_Na + i_K1 + i_to + i_Kur + i_Kr + i_Ks + i_B_Na + i_B_Ca + i_NaK +
               i_CaP + i_NaCa + i_Ca_L + i_st)) /
              Cm);
        next[1] = m_inf + (m - m_inf) * std::exp(-dt / tau_m);
        next[2] = h_inf + (h - h_inf) * std::exp(-dt / tau_h);
        next[3] = j_inf + (j - j_inf) * std::exp(-dt / tau_j);
        next[4] = oa_infinity + (oa - oa_infinity) * std::exp(-dt / tau_oa);
        next[5] = oi_infinity + (oi - oi_infinity) * std::exp(-dt / tau_oi);
        next[6] = ua_infinity + (ua - ua_infinity) * std::exp(-dt / tau_ua);
        next[7] = ui_infinity + (ui - ui_infinity) * std::exp(-dt / tau_ui);
        next[8] = xr_infinity + (xr - xr_infinity) * std::exp(-dt / tau_xr);
        next[9] = xs_infinity + (xs - xs_infinity) * std::exp(-dt / tau_xs);
        next[10] = d_infinity + (d - d_infinity) * std::exp(-dt / tau_d);
        next[11] = f_infinity + (f - f_infinity) * std::exp(-dt / tau_f);
        next[12] = f_Ca_infinity + (f_Ca - f_Ca_infinity) * std::exp(-dt / tau_f_Ca);
        next[13] = u_infinity + (u - u_infinity) * std::exp(-dt / tau_u);
        next[14] = v_infinity + (v - v_infinity) * std::exp(-dt / tau_v);
        next[15] = w_infinity + (w - w_infinity) * std::exp(-dt / tau_w);
        next[16] =
          Na_i + dt * ((((-3.0) * i_NaK) - ((3.0 * i_NaCa) + i_B_Na + i_Na)) / (V_i * F));
        next[17] = K_i +
          dt *
            (((2.0 * i_NaK) - (i_K1 + i_to + i_Kur + i_Kr + i_Ks + i_B_K)) / (V_i * F));
        next[18] = Ca_i + dt * (B1 / B2);
        next[19] = Ca_up + dt * (i_up - (i_up_leak + ((i_tr * V_rel) / V_up)));
        next[20] = Ca_rel +
          dt *
            ((i_tr - i_rel) *
              std::pow((1.0 + ((CSQN_max * Km_CSQN) / std::pow((Ca_rel + Km_CSQN), 2.0))),
                (-1.0)));
        return next;
      }
  };
}
#endif
