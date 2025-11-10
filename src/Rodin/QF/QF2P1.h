/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_QF2P1_H
#define RODIN_VARIATIONAL_QF_QF2P1_H

#include "Rodin/Geometry/GeometryIndexed.h"

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Quadrature formula with degree of exactness 1 and 2 points.
   *
   * This class implements quadrature formulas that are exact for polynomials 
   * up to degree 1. The number and location of quadrature points depends on 
   * the geometry of the element.
   */
  class QF2P1 final : public QuadratureFormulaBase
  {
    public:
      /// Parent class type
      using Parent = QuadratureFormulaBase;

      /**
       * @brief Constructs a QF2P1 quadrature formula for the given geometry.
       * @param g Geometry type (e.g., triangle, quadrilateral, etc.)
       */
      constexpr
      QF2P1(Geometry::Polytope::Type g)
        : Parent(g)
      {}

      /**
       * @brief Gets the number of quadrature points.
       * @return Number of quadrature points for this geometry
       */
      size_t getSize() const override
      {
        return s_size[getGeometry()];
      }

      /**
       * @brief Gets the i-th quadrature point coordinates.
       * @param i Index of the quadrature point
       * @return Reference to the spatial coordinates of the quadrature point
       */
      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        return s_points[getGeometry()][i];
      }

      /**
       * @brief Gets the weight of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Weight associated with the quadrature point
       */
      Real getWeight(size_t i) const override
      {
        return s_weights[getGeometry()][i];
      }

      /**
       * @brief Creates a copy of this quadrature formula.
       * @return Pointer to a new QF2P1 instance (ownership transferred to caller)
       */
      QF2P1* copy() const noexcept override
      {
        return new QF2P1(*this);
      }

    private:
      /// Static data: number of quadrature points for each geometry
      static const Geometry::GeometryIndexed<size_t> s_size;
      
      /// Static data: quadrature point coordinates for each geometry
      static const Geometry::GeometryIndexed<std::vector<Math::SpatialVector<Real>>> s_points;
      
      /// Static data: quadrature weights for each geometry
      static const Geometry::GeometryIndexed<Math::Vector<Real>> s_weights;
  };
}

#endif

