/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ValveLaw.h
 * @brief Cavity flux evaluator implementing the valve branch laws.
 *
 * Models the mitral valve, aortic valve, and leakage pathways
 * in the 0D lumped-parameter representation of the left ventricle
 * (Caruel et al. 2014, §4, valve model).
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_VALVELAW_H
#define RODIN_HEART_CCMLC2014_PHYSICS_VALVELAW_H

namespace Rodin::Heart::CCMLC2014::Physics
{
  /**
   * @brief Evaluates cavity inflow/outflow fluxes through the valve branches.
   *
   * Determines the valve state (mitral open, both closed, aortic open) based
   * on the relative pressures and computes the corresponding volume flux and
   * its derivatives with respect to ventricular and arterial pressures.
   *
   * The three valve regimes are:
   * - **Mitral open**: @f$ p_v \le p_{at} @f$ — filling phase.
   * - **Both closed**: @f$ p_{at} \le p_v \le p_{ar} @f$ — isovolumetric phase.
   * - **Aortic open**: @f$ p_v > p_{ar} @f$ — ejection phase.
   *
   * @tparam Input Model input parameter type.
   */
  template <class Input>
  class CavityFluxEvaluator
  {
    public:
      /**
       * @brief Construct with model input parameters.
       * @param[in] input Model parameters (valve conductances).
       */
      explicit CavityFluxEvaluator(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Evaluate cavity fluxes at the current and previous time levels.
       *
       * Populates @p data.cavityFluxCur, @p data.cavityFluxPrev, and their
       * pressure derivatives.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in,out] data Evaluation data with pressures already set.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);

        // --- Current step valve state ---
        const bool mitralOpenCurrent = data.pv <= data.pAtCur;
        const bool bothClosedCurrent =
          data.pAtCur <= data.pv && data.pv <= data.par;

        if (mitralOpenCurrent)
        {
          data.cavityFluxCur = m_input.Kat * (data.pv - data.pAtCur);
          data.dCavityFluxCur_dPv = m_input.Kat;
          data.dCavityFluxCur_dPar = Scalar(0);
        }
        else if (bothClosedCurrent)
        {
          data.cavityFluxCur = m_input.Kp * (data.pv - data.pAtCur);
          data.dCavityFluxCur_dPv = m_input.Kp;
          data.dCavityFluxCur_dPar = Scalar(0);
        }
        else
        {
          data.cavityFluxCur =
            m_input.Kar * (data.pv - data.par)
            + m_input.Kp * (data.par - data.pAtCur);
          data.dCavityFluxCur_dPv = m_input.Kar;
          data.dCavityFluxCur_dPar = -m_input.Kar + m_input.Kp;
        }

        // --- Previous step valve state ---
        const bool mitralOpenPrevious = data.pvPrev <= data.pAtPrev;
        const bool bothClosedPrevious =
          data.pAtPrev <= data.pvPrev && data.pvPrev <= data.parPrev;

        if (mitralOpenPrevious)
        {
          data.cavityFluxPrev =
            m_input.Kat * (data.pvPrev - data.pAtPrev);
        }
        else if (bothClosedPrevious)
        {
          data.cavityFluxPrev =
            m_input.Kp * (data.pvPrev - data.pAtPrev);
        }
        else
        {
          data.cavityFluxPrev =
            m_input.Kar * (data.pvPrev - data.parPrev)
            + m_input.Kp * (data.parPrev - data.pAtPrev);
        }
      }

    private:
      const Input& m_input;
  };
}

#endif
