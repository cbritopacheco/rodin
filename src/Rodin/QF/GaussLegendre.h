#ifndef RODIN_VARIATIONAL_QF_GAUSSLEGENDRE_H
#define RODIN_VARIATIONAL_QF_GAUSSLEGENDRE_H

#include "Rodin/Geometry/GeometryIndexed.h"

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   */
  class GaussLegendre final : public QuadratureFormulaBase
  {
    public:
      using Parent = QuadratureFormulaBase;

      constexpr
      GaussLegendre(Geometry::Polytope::Type g)
        : Parent(g)
      {}

      inline
      size_t getSize() const override
      {
        return 2;
      }

      inline
      const Math::SpatialVector& getPoint(size_t i) const override
      {
        assert(i == 0);
        return s_points[getGeometry()];
      }

      inline
      Scalar getWeight(size_t i) const override
      {
        assert(i == 0);
        return s_weights[getGeometry()].coeff(i);
      }

      inline
      GaussLegendre* copy() const noexcept override
      {
        return new GaussLegendre(*this);
      }

    private:
      static const Geometry::GeometryIndexed<Math::SpatialVector> s_points;
      static const Geometry::GeometryIndexed<Math::Vector> s_weights;
  };
}

#endif

