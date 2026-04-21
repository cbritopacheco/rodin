/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FlowPoint.h
 * @brief Central local constitutive input abstraction for fluid closures.
 */
#ifndef RODIN_FLUID_LOCAL_FLOWPOINT_H
#define RODIN_FLUID_LOCAL_FLOWPOINT_H

#include <any>
#include <cassert>
#include <functional>
#include <typeindex>
#include <unordered_map>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialMatrix.h"

namespace Rodin::Fluid
{
  /**
   * @brief Local fluid state bundle for quadrature-point constitutive evaluation.
   *
   * FlowPoint is the fluid analogue of Solid::ConstitutivePoint. It bundles:
   * - Optional geometric context (@c Geometry::Point)
   * - Velocity value @f$ \mathbf{u} @f$
   * - Velocity gradient @f$ \nabla \mathbf{u} @f$
   * - Optional pressure @f$ p @f$
   * - Optional pressure gradient @f$ \nabla p @f$
   * - Extensible typed auxiliary data via tags
   *
   * This lets rheology/closure code query all local information from one place,
   * both in real integration loops and in geometry-free unit tests.
   */
  class FlowPoint
  {
    public:
      using VelocityType = Math::SpatialVector<Real>;
      using VelocityGradientType = Math::SpatialMatrix<Real>;
      using PressureType = Real;
      using PressureGradientType = Math::SpatialVector<Real>;

      FlowPoint(
          const Geometry::Point& point,
          const VelocityType& velocity,
          const VelocityGradientType& velocityGradient)
        : m_point(std::cref(point)),
          m_velocity(velocity),
          m_velocityGradient(velocityGradient)
      {
        assert(m_velocityGradient.rows() == m_velocity.size());
        assert(m_velocityGradient.cols() == m_velocity.size());
      }

      FlowPoint(
          const VelocityType& velocity,
          const VelocityGradientType& velocityGradient)
        : m_velocity(velocity),
          m_velocityGradient(velocityGradient)
      {
        assert(m_velocityGradient.rows() == m_velocity.size());
        assert(m_velocityGradient.cols() == m_velocity.size());
      }

      FlowPoint(const FlowPoint&) = default;
      FlowPoint(FlowPoint&&) = default;
      FlowPoint& operator=(const FlowPoint&) = default;
      FlowPoint& operator=(FlowPoint&&) = default;

      const Optional<std::reference_wrapper<const Geometry::Point>>& getPoint() const
      {
        return m_point;
      }

      const VelocityType& getVelocity() const
      {
        return m_velocity;
      }

      const VelocityGradientType& getVelocityGradient() const
      {
        return m_velocityGradient;
      }

      FlowPoint& setPressure(const PressureType& pressure)
      {
        m_pressure = pressure;
        return *this;
      }

      const Optional<PressureType>& getPressure() const
      {
        return m_pressure;
      }

      bool hasPressure() const
      {
        return m_pressure.has_value();
      }

      FlowPoint& setPressureGradient(const PressureGradientType& pressureGradient)
      {
        assert(pressureGradient.size() == m_velocity.size());
        m_pressureGradient = pressureGradient;
        return *this;
      }

      const Optional<PressureGradientType>& getPressureGradient() const
      {
        return m_pressureGradient;
      }

      bool hasPressureGradient() const
      {
        return m_pressureGradient.has_value();
      }

      size_t getDimension() const
      {
        return m_velocity.size();
      }

      template <class Tag>
      FlowPoint& set(const typename Tag::Type& value)
      {
        m_aux[std::type_index(typeid(Tag))] = value;
        return *this;
      }

      template <class Tag>
      const typename Tag::Type& get() const
      {
        const auto it = m_aux.find(std::type_index(typeid(Tag)));
        assert(it != m_aux.end());
        return std::any_cast<const typename Tag::Type&>(it->second);
      }

      template <class Tag>
      bool has() const
      {
        return m_aux.find(std::type_index(typeid(Tag))) != m_aux.end();
      }

    private:
      Optional<std::reference_wrapper<const Geometry::Point>> m_point;
      VelocityType m_velocity;
      VelocityGradientType m_velocityGradient;
      Optional<PressureType> m_pressure;
      Optional<PressureGradientType> m_pressureGradient;
      std::unordered_map<std::type_index, std::any> m_aux;
  };
}

#endif
