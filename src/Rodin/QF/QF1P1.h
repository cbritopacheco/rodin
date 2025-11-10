/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_QF1P1_H
#define RODIN_VARIATIONAL_QF_QF1P1_H

#include "Rodin/Geometry/GeometryIndexed.h"

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Single-point quadrature formula with degree of exactness 1.
   *
   * This class implements quadrature formulas that use a single integration 
   * point (typically the centroid) and are exact for polynomials up to degree 1.
   * This is the simplest possible quadrature rule for each geometry type.
   */
  class QF1P1 final : public QuadratureFormulaBase
  {
    public:
      /// Parent class type
      using Parent = QuadratureFormulaBase;

      /**
       * @brief Constructs a QF1P1 quadrature formula for the given geometry.
       * @param g Geometry type (e.g., triangle, quadrilateral, etc.)
       */
      constexpr
      QF1P1(Geometry::Polytope::Type g)
        : Parent(g)
      {}

      /**
       * @brief Gets the number of quadrature points (always 1).
       * @return Always returns 1 (single point quadrature)
       */
      size_t getSize() const override
      {
        return 1;
      }

      /**
       * @brief Gets the single quadrature point coordinates.
       * @param i Index of the quadrature point (must be 0)
       * @return Reference to the spatial coordinates of the quadrature point
       */
      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        assert(i == 0);
        return s_points[getGeometry()];
      }

      /**
       * @brief Gets the weight of the single quadrature point.
       * @param i Index of the quadrature point (must be 0)
       * @return Weight associated with the quadrature point
       */
      Real getWeight(size_t i) const override
      {
        assert(i == 0);
        return s_weights[getGeometry()];
      }

      /**
       * @brief Creates a copy of this quadrature formula.
       * @return Pointer to a new QF1P1 instance (ownership transferred to caller)
       */
      QF1P1* copy() const noexcept override
      {
        return new QF1P1(*this);
      }

    private:
      /// Static data: single quadrature point coordinates for each geometry
      static const Geometry::GeometryIndexed<Math::SpatialVector<Real>> s_points;
      
      /// Static data: single quadrature point weight for each geometry
      static const Geometry::GeometryIndexed<Real> s_weights;
  };
}

#endif
