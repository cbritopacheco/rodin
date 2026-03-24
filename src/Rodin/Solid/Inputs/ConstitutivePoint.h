/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ConstitutivePoint.h
 * @brief Central constitutive input abstraction for hyperelastic formulations.
 *
 * The ConstitutivePoint bundles all data available at a quadrature point for
 * constitutive evaluation:
 * - KinematicState (deformation gradient and derived quantities)
 * - Reference coordinates @f$ \boldsymbol{\xi} @f$ (in the reference element)
 * - Physical coordinates @f$ \mathbf{x} @f$ (in the reference configuration)
 * - Geometry attribute / region id (for heterogeneous materials)
 * - Extensible typed auxiliary data (fiber directions, activation, etc.)
 *
 * This abstraction decouples constitutive laws from the integration scheme
 * and finite element space, enabling:
 * - Arbitrary quadrature rules (not just centroid)
 * - Arbitrary FE spaces (not just P1)
 * - Heterogeneous materials (region-dependent parameters)
 * - Anisotropic materials (fiber direction fields)
 * - Active materials (activation parameters)
 */
#ifndef RODIN_SOLID_INPUTS_CONSTITUTIVEPOINT_H
#define RODIN_SOLID_INPUTS_CONSTITUTIVEPOINT_H

#include <any>
#include <typeindex>
#include <unordered_map>
#include <cassert>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Types.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"

namespace Rodin::Solid
{
  /**
   * @brief Central data bundle for constitutive evaluation at a quadrature point.
   *
   * A ConstitutivePoint gathers all information available at a single
   * integration point and passes it to constitutive laws.  Laws may inspect
   * any subset of the data they need (kinematics, coordinates, region id,
   * or user-attached auxiliary data).
   *
   * ## Usage
   *
   * @code
   * ConstitutivePoint cp(state);
   * cp.setReferenceCoordinates(xi);
   * cp.setPhysicalCoordinates(x);
   * cp.setRegionId(3);
   * cp.setAuxiliaryData<FiberDirection>(fiberDir);
   *
   * law.setCache(cache, cp);
   * law.getFirstPiolaKirchhoffStress(P, cache, cp);
   * @endcode
   */
  class ConstitutivePoint
  {
    public:
      /**
       * @brief Constructs a constitutive point from a kinematic state.
       * @param state The kinematic state at this quadrature point
       */
      explicit ConstitutivePoint(const KinematicState& state)
        : m_state(state),
          m_regionId(0)
      {}

      ConstitutivePoint(const ConstitutivePoint&) = default;
      ConstitutivePoint(ConstitutivePoint&&) = default;
      ConstitutivePoint& operator=(const ConstitutivePoint&) = default;
      ConstitutivePoint& operator=(ConstitutivePoint&&) = default;

      /// @brief Gets the kinematic state.
      const KinematicState& getKinematicState() const { return m_state; }

      /**
       * @brief Sets the reference (parent element) coordinates @f$ \boldsymbol{\xi} @f$.
       * @param xi Reference coordinates
       * @returns Reference to this for chaining
       */
      ConstitutivePoint& setReferenceCoordinates(const Math::SpatialVector<Real>& xi)
      {
        m_xi = xi;
        return *this;
      }

      /// @brief Gets the reference coordinates @f$ \boldsymbol{\xi} @f$.
      const Math::SpatialVector<Real>& getReferenceCoordinates() const { return m_xi; }

      /**
       * @brief Sets the physical (material) coordinates @f$ \mathbf{x} @f$.
       * @param x Physical coordinates in the reference configuration
       * @returns Reference to this for chaining
       */
      ConstitutivePoint& setPhysicalCoordinates(const Math::SpatialVector<Real>& x)
      {
        m_x = x;
        return *this;
      }

      /// @brief Gets the physical coordinates @f$ \mathbf{x} @f$.
      const Math::SpatialVector<Real>& getPhysicalCoordinates() const { return m_x; }

      /**
       * @brief Sets the geometry attribute / region id.
       * @param id Region identifier (e.g., material region)
       * @returns Reference to this for chaining
       */
      ConstitutivePoint& setRegionId(Geometry::Attribute id)
      {
        m_regionId = id;
        return *this;
      }

      /// @brief Gets the geometry attribute / region id.
      Geometry::Attribute getRegionId() const { return m_regionId; }

      /**
       * @brief Stores typed auxiliary data (e.g., fiber direction, activation).
       *
       * @tparam Tag A type tag used as a key (e.g., a struct FiberDirection{};)
       * @param value The auxiliary data value
       * @returns Reference to this for chaining
       *
       * @code
       * struct FiberDirection {};
       * cp.setAuxiliaryData<FiberDirection>(fiberVec);
       * @endcode
       */
      template <class Tag, class T>
      ConstitutivePoint& setAuxiliaryData(T&& value)
      {
        m_aux[std::type_index(typeid(Tag))] = std::forward<T>(value);
        return *this;
      }

      /**
       * @brief Retrieves typed auxiliary data by tag.
       *
       * @tparam Tag The type tag used when storing the data
       * @tparam T The expected type of the stored data
       * @returns Const reference to the stored value
       *
       * @code
       * const auto& fiber = cp.getAuxiliaryData<FiberDirection, Math::SpatialVector<Real>>();
       * @endcode
       */
      template <class Tag, class T>
      const T& getAuxiliaryData() const
      {
        auto it = m_aux.find(std::type_index(typeid(Tag)));
        assert(it != m_aux.end());
        return std::any_cast<const T&>(it->second);
      }

      /**
       * @brief Checks whether auxiliary data with the given tag exists.
       * @tparam Tag The type tag to check
       * @returns True if auxiliary data with this tag has been set
       */
      template <class Tag>
      bool hasAuxiliaryData() const
      {
        return m_aux.find(std::type_index(typeid(Tag))) != m_aux.end();
      }

    private:
      const KinematicState& m_state;
      Math::SpatialVector<Real> m_xi;
      Math::SpatialVector<Real> m_x;
      Geometry::Attribute m_regionId;
      std::unordered_map<std::type_index, std::any> m_aux;
  };
}

#endif
