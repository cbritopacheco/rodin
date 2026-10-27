/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Reference-space degree bounds for Lagrange finite elements.
 */
#ifndef RODIN_VARIATIONAL_LAGRANGEFINITEELEMENT_H
#define RODIN_VARIATIONAL_LAGRANGEFINITEELEMENT_H

#include <limits>

#include "FiniteElement.h"

namespace Rodin::Variational
{
  /**
   * @ingroup FiniteElements
   * @brief Common reference-space properties of Lagrange element families.
   *
   * The class describes the reference function spaces. It does not prescribe
   * a basis representation or a degree-of-freedom layout.
   * @tparam Derived Concrete finite element implementing the CRTP interface.
   */
  template <class Derived>
  class LagrangeFiniteElementBase : public FiniteElementBase<Derived>
  {
    public:
      /// @brief Generic finite element base.
      using Parent = FiniteElementBase<Derived>;

      /// @brief Inherits the generic finite element constructors.
      using Parent::Parent;

      /**
       * @brief Bounds a native Lagrange reference space in selected directions.
       * The returned bound describes the space with the supplied family index;
       * it does not necessarily describe this object's own family index.
       * @param order Family index of the Lagrange space.
       * @param directions Reference-coordinate indices.
       * @returns A directional bound, or no bound for invalid directions,
       * unsupported selections, or an unrepresentable degree.
       */
      Optional<size_t> getLagrangeDegree(
        size_t order, std::initializer_list<size_t> directions) const
      {
        using G = Geometry::Polytope::Type;
        const auto geometry = this->getGeometry();
        unsigned valid = 0;
        switch (geometry)
        {
          case G::Point:
            valid = 0;
            break;
          case G::Segment:
            valid = 1;
            break;
          case G::Triangle:
          case G::Quadrilateral:
            valid = 3;
            break;
          case G::Tetrahedron:
          case G::Hexahedron:
          case G::Wedge:
          case G::Pyramid:
            valid = 7;
            break;
        }
        unsigned mask = 0;
        for (const size_t direction : directions)
        {
          if (direction >= 3 || !(valid & (1u << direction)))
            return std::nullopt;
          mask |= 1u << direction;
        }
        if (mask == 0 || order == 0)
          return 0;
        size_t factors = 0;
        switch (geometry)
        {
          case G::Segment:
          case G::Triangle:
          case G::Tetrahedron:
            return order;
          case G::Quadrilateral:
          case G::Hexahedron:
            factors = (mask & 1u ? 1 : 0) + (mask & 2u ? 1 : 0) + (mask & 4u ? 1 : 0);
            break;
          case G::Wedge:
            factors = (mask & 3u ? 1 : 0) + (mask & 4u ? 1 : 0);
            break;
          case G::Pyramid:
            // Pyramid basis components are polynomial in x and y when z is
            // held fixed; their coefficients can be rational functions of z.
            if (mask & 4u)
              return std::nullopt;
            factors = (mask & 1u ? 1 : 0) + (mask & 2u ? 1 : 0);
            break;
          case G::Point:
            return 0;
        }
        if (factors == 0 || order > std::numeric_limits<size_t>::max() / factors)
          return std::nullopt;
        return order * factors;
      }
  };
}

#endif
