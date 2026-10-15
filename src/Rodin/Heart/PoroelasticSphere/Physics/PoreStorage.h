/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PoreStorage.h
 * @brief Pore storage pressure @f$ p_P(\Phi) = \Psi_P'(\Phi) @f$: linear
 * (@f$ K_\Phi @f$) or Bruinsma vascular compliance (Ruz et al., Sec. 4.4).
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_PHYSICS_PORESTORAGE_H
#define RODIN_HEART_POROELASTICSPHERE_PHYSICS_PORESTORAGE_H

#include <cmath>
#include <utility>

namespace Rodin::Heart::PoroelasticSphere::Physics
{
  /**
   * @brief Pore storage pressure and its tangent @f$ (p_P, p_P') @f$ at
   * porosity @f$ \Phi @f$.
   *
   * - Linear: @f$ p_P = K_\Phi(\Phi - \phi_0) @f$.
   * - Bruinsma: @f$ p_P = -p_0 \ln(1 - (\Phi-\phi_0)/\Delta\Phi_{max}) @f$,
   *   @f$ p_P' = p_0 / (\Delta\Phi_{max} - (\Phi - \phi_0)) @f$. Past
   *   @f$ (1-\varepsilon)\Delta\Phi_{max} @f$ the law is continued linearly
   *   (C1) so Newton iterates never evaluate the logarithm at or beyond the
   *   vascular capacity.
   */
  template <class Input, class Scalar>
  std::pair<Scalar, Scalar> poreStoragePressure(const Input& input, const Scalar phi)
  {
    const Scalar x = phi - input.phi0;
    if (input.poreStorage != Input::PoreStorageLaw::Bruinsma)
      return {input.KPhi * x, input.KPhi};

    const Scalar D = input.dPhiMax;
    const Scalar p0 = input.p0;
    const Scalar xmax = (Scalar(1) - input.poreStorageCutoff) * D;
    if (x < xmax)
      return {-p0 * std::log(Scalar(1) - x / D), p0 / (D - x)};

    const Scalar slope = p0 / (D - xmax);
    return {-p0 * std::log(Scalar(1) - xmax / D) + slope * (x - xmax), slope};
  }
}

#endif
