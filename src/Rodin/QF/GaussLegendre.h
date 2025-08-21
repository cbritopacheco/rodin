#ifndef RODIN_VARIATIONAL_QF_GAUSSLEGENDRE_H
#define RODIN_VARIATIONAL_QF_GAUSSLEGENDRE_H

#include "Rodin/Geometry/GeometryIndexed.h"

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Gauss-Legendre quadrature formula.
   *
   * This class implements Gauss-Legendre quadrature rules, which are optimal
   * for integrating polynomials over intervals. The quadrature points are 
   * the roots of Legendre polynomials and the weights are chosen to achieve
   * maximum algebraic degree of exactness.
   *
   * @note Current implementation appears to have inconsistencies between
   *       getSize() returning 2 and getPoint() only accepting i=0.
   */
  class GaussLegendre final : public QuadratureFormulaBase
  {
    public:
      /// Parent class type
      using Parent = QuadratureFormulaBase;

      /**
       * @brief Constructs a Gauss-Legendre quadrature formula for the given geometry.
       * @param g Geometry type (should be appropriate for Gauss-Legendre rule)
       */
      constexpr
      GaussLegendre(Geometry::Polytope::Type g)
        : Parent(g)
      {}

      /**
       * @brief Gets the number of quadrature points.
       * @return Number of Gauss-Legendre quadrature points
       */
      inline
      size_t getSize() const override
      {
        return 2;
      }

      /**
       * @brief Gets the i-th quadrature point coordinates.
       * @param i Index of the quadrature point
       * @return Reference to the spatial coordinates of the quadrature point
       * @note Current implementation only supports i=0
       */
      inline
      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        assert(i == 0);
        return s_points[getGeometry()];
      }

      /**
       * @brief Gets the weight of the i-th quadrature point.
       * @param i Index of the quadrature point  
       * @return Weight associated with the quadrature point
       */
      inline
      Real getWeight(size_t i) const override
      {
        assert(i == 0);
        return s_weights[getGeometry()].coeff(i);
      }

      /**
       * @brief Creates a copy of this quadrature formula.
       * @return Pointer to a new GaussLegendre instance (ownership transferred to caller)
       */
      inline
      GaussLegendre* copy() const noexcept override
      {
        return new GaussLegendre(*this);
      }

    private:
      /// Static data: Gauss-Legendre quadrature point coordinates for each geometry
      static const Geometry::GeometryIndexed<Math::SpatialVector<Real>> s_points;
      
      /// Static data: Gauss-Legendre quadrature weights for each geometry
      static const Geometry::GeometryIndexed<Math::Vector<Real>> s_weights;
  };
}

#endif

