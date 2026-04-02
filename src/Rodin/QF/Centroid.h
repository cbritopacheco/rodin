/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_CENTROID_H
#define RODIN_VARIATIONAL_QF_CENTROID_H

/**
 * @file
 * @brief Defines the Centroid single-point quadrature formula.
 */

#include "QuadratureFormula.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Single-point centroid quadrature formula.
   *
   * This class implements a quadrature formula that uses a single integration
   * point located at the centroid (barycenter) of the reference polytope.
   * The quadrature is exact for polynomials up to degree 1 (linear functions).
   *
   * ## Mathematical Foundation
   *
   * The centroid rule approximates integrals as:
   * @f[
   *   \int_K f(x) \, dx \approx w \cdot f(c)
   * @f]
   * where @f$ c @f$ is the centroid of the polytope @f$ K @f$ and @f$ w @f$
   * is the weight equal to the volume (measure) of the reference element.
   *
   * ## Reference Element Centroids
   *
   * The quadrature points are located at:
   * - Point: @f$ (0) @f$
   * - Segment: @f$ (0.5) @f$
   * - Triangle: @f$ (1/3, 1/3) @f$
   * - Quadrilateral: @f$ (0.5, 0.5) @f$
   * - Tetrahedron: @f$ (0.25, 0.25, 0.25) @f$
   * - Wedge: @f$ (1/3, 1/3, 0.5) @f$
   * - Hexahedron: @f$ (0.5, 0.5, 0.5) @f$
   *
   * ## Reference Element Weights
   *
   * The weights correspond to the reference element measures:
   * - Point: @f$ 1 @f$
   * - Segment: @f$ 1 @f$
   * - Triangle: @f$ 1/2 @f$
   * - Quadrilateral: @f$ 1 @f$
   * - Tetrahedron: @f$ 1/6 @f$
   * - Wedge: @f$ 1/2 @f$
   * - Hexahedron: @f$ 1 @f$
   *
   * This is the simplest possible quadrature rule for each geometry type
   * and is computationally efficient when low accuracy is acceptable.
   *
   * @see QuadratureFormulaBase
   * @see GaussLegendre
   */
  class Centroid final : public QuadratureFormulaBase
  {
    public:
      /// Parent class type
      using Parent = QuadratureFormulaBase;

      /**
       * @brief Constructs a centroid quadrature formula for the given geometry.
       * @param g Geometry type (e.g., Triangle, Quadrilateral, Tetrahedron, etc.)
       */
      constexpr
      Centroid(Geometry::Polytope::Type g)
        : m_geometry(g)
      {}

      constexpr
      Geometry::Polytope::Type getGeometry() const
      {
        return m_geometry;
      }

      /**
       * @brief Gets the number of quadrature points (always 1).
       * @return Always returns 1 (single point quadrature)
       */
      size_t getSize() const override
      {
        return 1;
      }

      /**
       * @brief Gets the single quadrature point coordinates (the centroid).
       * @param i Index of the quadrature point (must be 0)
       * @return Reference to the centroid coordinates in reference space
       */
      const Math::SpatialVector<Real>& getPoint(size_t i) const override;

      /**
       * @brief Gets the weight of the single quadrature point.
       * @param i Index of the quadrature point (must be 0)
       * @return Weight equal to the reference element measure
       */
      Real getWeight(size_t i) const override;

      /**
       * @brief Creates a copy of this quadrature formula.
       * @return Pointer to a new Centroid instance (caller takes ownership)
       */
      Centroid* copy() const noexcept override
      {
        return new Centroid(*this);
      }

    private:
      Geometry::Polytope::Type m_geometry;
  };
}

#endif
