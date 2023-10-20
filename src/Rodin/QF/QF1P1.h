#ifndef RODIN_VARIATIONAL_QF_QF1P1_H
#define RODIN_VARIATIONAL_QF_QF1P1_H

#include "Rodin/Geometry/GeometryIndexed.h"
#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   */
  class QF1P1 final : public QuadratureFormulaBase
  {
    public:
      using Parent = QuadratureFormulaBase;

      constexpr
      QF1P1(Geometry::Polytope::Type g)
        : Parent(g)
      {}

      inline
      size_t getSize() const override
      {
        return s_weights[getGeometry()].size();
      }

      inline
      const Math::SpatialVector& getPoint(size_t i) const override
      {
        return s_points[getGeometry()][i];
      }

      inline
      Scalar getWeight(size_t i) const override
      {
        return s_weights[getGeometry()].coeff(i);
      }

      inline
      QF1P1* copy() const noexcept override
      {
        return new QF1P1(*this);
      }

    private:
      static const Geometry::GeometryIndexed<std::vector<Math::SpatialVector>> s_points;
      static const Geometry::GeometryIndexed<Math::Vector> s_weights;
  };
}

#endif
