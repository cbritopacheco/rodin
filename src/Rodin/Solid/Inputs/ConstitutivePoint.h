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
 * constitutive evaluation.  It is built by composition over Geometry::Point,
 * which provides the geometric evaluation context (reference/physical
 * coordinates, polytope, Jacobian), while ConstitutivePoint extends it with:
 *
 * - KinematicState (deformation gradient and derived quantities)
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
#include <functional>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Types.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"

namespace Rodin::Solid
{
  /**
   * @brief Central data bundle for constitutive evaluation at a quadrature point.
   *
   * A ConstitutivePoint composes a Geometry::Point (geometric context) with a
   * KinematicState (deformation measures) and extensible typed auxiliary data.
   * Constitutive laws receive a ConstitutivePoint and may inspect any subset
   * of the data they need.
   *
   * ## Usage (with geometric context, typical in integrators)
   *
   * @code
   * Geometry::Point pt(polytope, rc);
   * KinematicState state(d);
   * state.setDisplacementGradient(H);
   *
   * ConstitutivePoint cp(pt, state);
   * cp.setAuxiliaryData<FiberDirection>(fiberDir);
   *
   * law.setCache(cache, cp);
   * law.getFirstPiolaKirchhoffStress(P, cache, cp);
   * @endcode
   *
   * ## Usage (without geometry, for unit testing)
   *
   * @code
   * KinematicState state(2);
   * state.setDisplacementGradient(H);
   * ConstitutivePoint cp(state);
   * @endcode
   */
  class ConstitutivePoint
  {
    public:
      /**
       * @brief Constructs a constitutive point from a geometric point and kinematic state.
       *
       * This is the primary constructor used by integrators.  The Geometry::Point
       * provides reference/physical coordinates and the polytope (from which
       * the region id can be queried).
       *
       * @param point The geometric evaluation point (reference to a local variable is valid
       *              as long as the ConstitutivePoint does not outlive it)
       * @param state The kinematic state at this quadrature point
       */
      explicit ConstitutivePoint(const Geometry::Point& point, const KinematicState& state)
        : m_point(std::cref(point)),
          m_state(state)
      {}

      /**
       * @brief Constructs a constitutive point from a kinematic state only.
       *
       * Use this constructor for unit tests or contexts where geometric
       * context is not needed.  Coordinate and region queries will assert
       * if called without a Geometry::Point.
       *
       * @param state The kinematic state
       */
      explicit ConstitutivePoint(const KinematicState& state)
        : m_state(state)
      {}

      ConstitutivePoint(const ConstitutivePoint&) = default;
      ConstitutivePoint(ConstitutivePoint&&) = default;
      ConstitutivePoint& operator=(const ConstitutivePoint&) = default;
      ConstitutivePoint& operator=(ConstitutivePoint&&) = default;

      /// @brief Gets the kinematic state.
      const KinematicState& getKinematicState() const { return m_state; }

      /// @brief Checks whether a Geometry::Point is available.
      bool hasPoint() const { return m_point.has_value(); }

      /**
       * @brief Gets the underlying Geometry::Point.
       * @pre hasPoint() must be true.
       */
      const Geometry::Point& getPoint() const
      {
        assert(m_point);
        return m_point->get();
      }

      /**
       * @brief Gets the reference (parent element) coordinates @f$ \boldsymbol{\xi} @f$.
       *
       * Delegates to the underlying Geometry::Point.
       * @pre hasPoint() must be true.
       */
      const auto& getReferenceCoordinates() const
      {
        return getPoint().getReferenceCoordinates();
      }

      /**
       * @brief Gets the physical (material) coordinates @f$ \mathbf{x} @f$.
       *
       * Delegates to the underlying Geometry::Point.
       * @pre hasPoint() must be true.
       */
      const auto& getPhysicalCoordinates() const
      {
        return getPoint().getPhysicalCoordinates();
      }

      /**
       * @brief Gets the geometry attribute / region id from the polytope.
       *
       * Delegates to the underlying Geometry::Point's polytope.
       * @pre hasPoint() must be true.
       * @returns The polytope attribute, or empty if no attribute is set.
       */
      Optional<Geometry::Attribute> getRegionId() const
      {
        return getPoint().getPolytope().getAttribute();
      }

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
      Optional<std::reference_wrapper<const Geometry::Point>> m_point;
      const KinematicState& m_state;
      std::unordered_map<std::type_index, std::any> m_aux;
  };
}

#endif
