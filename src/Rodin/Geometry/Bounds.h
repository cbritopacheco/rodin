/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_BOUNDS_H
#define RODIN_GEOMETRY_BOUNDS_H

/**
 * @file
 * @brief Conservative support bounds for mapped reference cells.
 */

#include <array>
#include <cmath>
#include <vector>
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Geometry/Polytope.h"
#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief A conservative convex enclosure described by support bounds.
   *
   * For an enclosed set @f$C@f$, getSupport() bounds
   * @f$h_C(d)=\sup_{x\in C}d\cdot x@f$. Directions need not be normalized.
   *
   * ## Architecture
   *
   * Native first-order maps use their vertex hull without allocation. Higher
   * orders convert a certified reference-space representation into a positive
   * partition of unity, including the rational collapsed pyramid basis.
   * Reference conversions are shared; physical coefficients are owned.
   * Conversion roundoff is included in every support query. An unreliable
   * conversion gives infinite support, never a sampled membership filter.
   */
  class Bounds
  {
    public:
      /**
       * @brief Constructs the enclosure of a certified reference-space map.
       *
       * The caller must establish containment in the degree-p Bernstein
       * reference space, whether simplex, tensor, triangle-segment product or rational
       * collapsed pyramid. This is a representation bound, not a family index
       * inferred from an arbitrary transformation's samples.
       * @param transformation Coordinate map to enclose.
       * @param geometry Reference geometry of that map.
       * @param degree Certified degree of the containing control space.
       */
      Bounds(const PolytopeTransformation& transformation, Polytope::Type geometry,
        size_t degree);

      /**
       * @brief Returns an outward-rounded support bound.
       * @param direction Physical direction, not necessarily normalized.
       * @returns An upper support bound, or infinity for an uncertified conversion.
       */
      Real getSupport(const Math::SpatialVector<Real>& direction) const
      {
        if (direction.size() != m_dimension)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Direction dimension does not match the enclosure dimension."
            << Alert::Raise;
        }
        size_t axis = m_dimension;
        Real length = 0;
        for (size_t i = 0; i < m_dimension; ++i)
        {
          const Real component = direction[i];
          length += std::abs(component);
          if (component != 0)
          {
            if (axis != m_dimension || std::abs(component) != 1)
              axis = m_dimension + 1;
            else
              axis = i;
          }
        }
        if (axis < m_dimension)
          return direction[axis] > 0 ? m_upper[axis] : -m_lower[axis];
        if (length == 0)
          return 0;
        return getDirectionalSupport(direction, length);
      }

    private:
      /**
       * @brief Returns an outward-rounded support bound in a non-axis direction.
       * @param direction Physical direction, not necessarily normalized.
       * @param length Sum of absolute direction components.
       * @returns An upper support bound, or infinity if it cannot be certified.
       */
      Real getDirectionalSupport(
        const Math::SpatialVector<Real>& direction, Real length) const;

      /// Conversion allowance in machine epsilons, scaled by conditioning.
      static constexpr Real ControlRoundoffFactor = Real(32);

      /// Cubic tensor maps use separable line conversions.
      static constexpr size_t MinSeparableTensorDegree = 3;

      struct ControlBasis
      {
        std::vector<Math::SpatialPoint> samples;
        Math::Matrix<Real> conversion;
        Real roundoffAmplification = 0;
        size_t tensorFactors = 0;
        bool valid = false;
      };

      static Real binomial(size_t n, size_t i);
      static Real bernstein(size_t n, size_t i, Real x);
      static Real multinomial(size_t k, const std::array<size_t, 4>& alpha, size_t count);
      static Real simplexBernstein(size_t k, const std::array<size_t, 4>& alpha,
        size_t rdim, const Math::SpatialPoint& r);
      static Real controlBasisValue(Polytope::Type g, size_t k,
        const std::array<size_t, 4>& mode, const Math::SpatialPoint& r);
      static void makeControlLattice(Polytope::Type g, size_t k,
        std::vector<std::array<size_t, 4>>& modes,
        std::vector<Math::SpatialPoint>& samples);
      static const ControlBasis& getControlBasis(Polytope::Type g, size_t k);

      std::array<Math::SpatialPoint, RODIN_MAXIMUM_POLYTOPE_VERTICES> m_vertices;
      size_t m_count;
      Math::Matrix<Real> m_control;
      std::array<Real, Math::SpatialVector<Real>::MaxSize> m_lower;
      std::array<Real, Math::SpatialVector<Real>::MaxSize> m_upper;
      Real m_error;
      Real m_scale;
      size_t m_dimension;
  };
}
#endif
