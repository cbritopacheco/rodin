#ifndef RODIN_VARIATIONAL_QF_QF1P1_H
#define RODIN_VARIATIONAL_QF_QF1P1_H

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  class QF1P1 final : public QuadratureFormulaBase
  {
    public:
      using Parent = QuadratureFormulaBase;

      constexpr
      QF1P1(Geometry::Polytope::Geometry g)
        : Parent(g)
      {}

      const Math::Matrix& getPoints() const override
      {
        return s_points[getGeometry()];
      }

      const Math::Vector& getWeights() const override
      {
        return s_weights[getGeometry()];
      }

    private:
      static const Geometry::GeometryIndexed<Math::Matrix> s_points;
      static const Geometry::GeometryIndexed<Math::Vector> s_weights;
  };
}

#endif
