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
 * Auxiliary data tags define their associated value type via a nested
 * @c Type alias, so @c set<Tag>(value) and @c get<Tag>() do not require
 * the caller to repeat the value type.
 *
 * This abstraction decouples constitutive laws from the integration scheme
 * and finite element space, enabling:
 * - Arbitrary quadrature rules (not just centroid)
 * - Arbitrary FE spaces (not just P1)
 * - Heterogeneous materials (region-dependent parameters)
 * - Anisotropic materials (fiber direction fields)
 * - Active materials (activation parameters)
 */
#ifndef RODIN_SOLID_LOCAL_CONSTITUTIVEPOINT_H
#define RODIN_SOLID_LOCAL_CONSTITUTIVEPOINT_H

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
   * @brief Standard tag types for auxiliary constitutive data.
   *
   * Each tag defines a nested @c Type alias indicating the value type
   * stored and retrieved via @c ConstitutivePoint::set / @c get.
   *
   * @code
   * cp.set<Tags::FiberDirection>(fiberVec);
   * const auto& fiber = cp.get<Tags::FiberDirection>();
   * @endcode
   */
  namespace Tags
  {
    /// @brief Tag for the fiber direction vector.
    struct FiberDirection
    {
      using Type = Math::SpatialVector<Real>;
    };

    /// @brief Tag for the sheet direction vector.
    struct SheetDirection
    {
      using Type = Math::SpatialVector<Real>;
    };

    /// @brief Tag for the sheet-normal direction vector.
    struct SheetNormalDirection
    {
      using Type = Math::SpatialVector<Real>;
    };

    /// @brief Tag for the activation parameter.
    struct Activation
    {
      using Type = Real;
    };
  }

  /**
   * @brief Central data bundle for constitutive evaluation at a quadrature point.
   *
   * A ConstitutivePoint composes a Geometry::Point (geometric context) with a
   * KinematicState (deformation measures) and extensible typed auxiliary data.
   * Constitutive laws receive a ConstitutivePoint and may inspect any subset
   * of the data they need.
   *
   * Geometric context (reference/physical coordinates, region id) is accessed
   * through @c getPoint(), which returns an optional reference to the
   * underlying Geometry::Point.
   *
   * ## Usage (with geometric context, typical in integrators)
   *
   * @code
   * Geometry::Point pt(polytope, rc);
   * KinematicState state(d);
   * state.setDisplacementGradient(H);
   *
   * ConstitutivePoint cp(pt, state);
   * cp.set<Tags::FiberDirection>(fiberDir);
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
       * @warning The caller must ensure the Geometry::Point and KinematicState
       *          outlive this ConstitutivePoint.  Typically all are local
       *          variables in the same quadrature-point loop.
       *
       * @param point The geometric evaluation point
       * @param state The kinematic state at this quadrature point
       */
      explicit ConstitutivePoint(const Geometry::Point& point, const KinematicState& state)
        : m_point(std::cref(point)),
          m_state(std::cref(state))
      {}

      /**
       * @brief Constructs a constitutive point from a kinematic state only.
       *
       * Use this constructor for unit tests or contexts where geometric
       * context is not needed.  getPoint() will return an empty optional.
       *
       * @param state The kinematic state
       */
      explicit ConstitutivePoint(const KinematicState& state)
        : m_state(std::cref(state))
      {}

      ConstitutivePoint(const ConstitutivePoint&) = default;
      ConstitutivePoint(ConstitutivePoint&&) = default;
      ConstitutivePoint& operator=(const ConstitutivePoint&) = default;
      ConstitutivePoint& operator=(ConstitutivePoint&&) = default;

      /// @brief Gets the kinematic state.
      const KinematicState& getKinematicState() const { return m_state.get(); }

      /**
       * @brief Gets the underlying Geometry::Point, if available.
       * @returns An optional reference to the Geometry::Point, or empty if
       *          constructed without geometric context.
       */
      const Optional<std::reference_wrapper<const Geometry::Point>>& getPoint() const
      {
        return m_point;
      }

      /**
       * @brief Stores typed auxiliary data (e.g., fiber direction, activation).
       *
       * The value type is deduced from the tag's @c Type alias.
       *
       * @tparam Tag A type tag with a nested @c Type alias
       * @param value The auxiliary data value
       * @returns Reference to this for chaining
       *
       * @code
       * cp.set<Tags::FiberDirection>(fiberVec);
       * cp.set<Tags::Activation>(0.5);
       * @endcode
       */
      template <class Tag>
      ConstitutivePoint& set(const typename Tag::Type& value)
      {
        m_aux[std::type_index(typeid(Tag))] = value;
        return *this;
      }

      /**
       * @brief Retrieves typed auxiliary data by tag.
       *
       * The return type is deduced from the tag's @c Type alias.
       *
       * @tparam Tag The type tag used when storing the data
       * @returns Const reference to the stored value
       *
       * @code
       * const auto& fiber = cp.get<Tags::FiberDirection>();
       * Real activation = cp.get<Tags::Activation>();
       * @endcode
       */
      template <class Tag>
      const typename Tag::Type& get() const
      {
        auto it = m_aux.find(std::type_index(typeid(Tag)));
        assert(it != m_aux.end());
        return std::any_cast<const typename Tag::Type&>(it->second);
      }

      /**
       * @brief Checks whether auxiliary data with the given tag exists.
       * @tparam Tag The type tag to check
       * @returns True if auxiliary data with this tag has been set
       */
      template <class Tag>
      bool has() const
      {
        return m_aux.find(std::type_index(typeid(Tag))) != m_aux.end();
      }

    private:
      Optional<std::reference_wrapper<const Geometry::Point>> m_point;
      std::reference_wrapper<const KinematicState> m_state;
      std::unordered_map<std::type_index, std::any> m_aux;
  };
}

#endif
