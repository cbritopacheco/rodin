#ifndef RODIN_VARIATIONAL_QF_QF1P1_H
#define RODIN_VARIATIONAL_QF_QF1P1_H

#include "Rodin/Geometry/GeometryIndexed.h"
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

      inline
      size_t getSize() const override
      {
        return s_weights[getGeometry()].size();
      }

      inline
      const Math::Vector& getPoint(size_t i) const override
      {
        return s_points[getGeometry()][i];
      }

      inline
      Scalar getWeight(size_t i) const override
      {
        return s_weights[getGeometry()].coeff(i);
      }

    private:
      static const Geometry::GeometryIndexed<std::vector<Math::Vector>> s_points;
      static const Geometry::GeometryIndexed<Math::Vector> s_weights;
  };
}

#endif
