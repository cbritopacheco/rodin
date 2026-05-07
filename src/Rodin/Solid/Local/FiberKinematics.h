/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FiberKinematics.h
 * @brief Preferred-direction data at a constitutive point.
 */
#ifndef RODIN_SOLID_LOCAL_FIBERKINEMATICS_H
#define RODIN_SOLID_LOCAL_FIBERKINEMATICS_H

#include <cmath>

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"

namespace Rodin::Solid
{
  /**
   * @brief Preferred-direction kinematics at a constitutive point.
   */
  class FiberKinematics
  {
    public:
      FiberKinematics()
        : m_I4(1.0),
          m_strain(0.0)
      {}

      FiberKinematics(
          const KinematicState& state,
          const Math::SpatialVector<Real>& direction)
        : m_direction(normalize(direction)),
          m_tensor(dyad(m_direction)),
          m_current(state.getDeformationGradient() * m_direction),
          m_I4(m_direction.dot(state.getRightCauchyGreenTensor() * m_direction)),
          m_strain(0.5 * (m_I4 - 1.0))
      {}

      explicit FiberKinematics(const ConstitutivePoint& cp)
        : FiberKinematics(cp.getKinematicState(), direction(cp))
      {}

      const Math::SpatialVector<Real>& direction() const
      {
        return m_direction;
      }

      const Math::SpatialMatrix<Real>& tensor() const
      {
        return m_tensor;
      }

      const Math::SpatialVector<Real>& current() const
      {
        return m_current;
      }

      Real I4() const
      {
        return m_I4;
      }

      Real strain() const
      {
        return m_strain;
      }

      Real dStrain(const Math::SpatialMatrix<Real>& dF) const
      {
        return m_current.dot(dF * m_direction);
      }

      static Math::SpatialVector<Real> direction(const ConstitutivePoint& cp)
      {
        const auto d =
          static_cast<std::uint8_t>(cp.getKinematicState().getDimension());
        Math::SpatialVector<Real> a(d);
        if (cp.has<Tags::FiberDirection>())
          a = cp.get<Tags::FiberDirection>();
        else
        {
          a.setZero();
          a[0] = 1.0;
        }
        return normalize(a);
      }

      static Math::SpatialMatrix<Real> dyad(const Math::SpatialVector<Real>& a)
      {
        const auto d = static_cast<std::uint8_t>(a.size());
        Math::SpatialMatrix<Real> A(d, d);
        A.setZero();
        for (std::uint8_t i = 0; i < d; ++i)
          for (std::uint8_t j = 0; j < d; ++j)
            A(i, j) = a[i] * a[j];
        return A;
      }

      static Math::SpatialVector<Real> normalize(Math::SpatialVector<Real> a)
      {
        const Real n = std::sqrt(a.dot(a));
        if (n > 0.0)
          a = (1.0 / n) * a;
        return a;
      }

    private:
      Math::SpatialVector<Real> m_direction;
      Math::SpatialMatrix<Real> m_tensor;
      Math::SpatialVector<Real> m_current;
      Real m_I4;
      Real m_strain;
  };
}

#endif
