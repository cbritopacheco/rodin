/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CauchyStress.h
 * @brief Cauchy stress helper from pressure and constitutive deviatoric stress.
 */
#ifndef RODIN_FLUID_FIELDS_CAUCHYSTRESS_H
#define RODIN_FLUID_FIELDS_CAUCHYSTRESS_H

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Fluid/Local/FlowPoint.h"

namespace Rodin::Fluid
{
  /**
   * @brief Combines pressure and deviatoric stress into the Cauchy stress.
   *
   * @f[
   *   \boldsymbol{\sigma} = -p\mathbf{I} + \boldsymbol{\tau}
   * @f]
   * where @f$\boldsymbol{\tau}@f$ is provided by a rheology law.
   *
   * This class intentionally does not compute constitutive response itself.
   */
  template <class LawDerived>
  class CauchyStress
  {
    public:
      using LawType = LawDerived;
      using CacheType = typename LawType::Cache;

      explicit CauchyStress(const LawType& law)
        : m_lawRef(law)
      {}

      void getCauchyStress(
          Math::SpatialMatrix<Real>& sigma,
          const CacheType& cache,
          const FlowPoint& fp) const
      {
        Math::SpatialMatrix<Real> tau;
        m_lawRef.getDeviatoricStress(tau, cache, fp);

        const size_t d = fp.getDimension();
        Math::SpatialMatrix<Real> I;
        I.resize(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        I.setIdentity();

        const Real p = fp.hasPressure() ? fp.getPressure().value() : Real(0);
        sigma = tau - p * I;
      }

      const LawType& getLaw() const
      {
        return m_lawRef;
      }

    private:
      const LawType& m_lawRef;
  };
}

#endif
