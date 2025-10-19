#ifndef RODIN_VARIATIONAL_QF_QFVertexP1_H
#define RODIN_VARIATIONAL_QF_QFVertexP1_H

#include "Rodin/Geometry/GeometryIndexed.h"
#include "QuadratureFormula.h"

namespace Rodin::QF
{
  class QFVertexP1 final : public QuadratureFormulaBase
  {
    public:
      using Parent = QuadratureFormulaBase;

      constexpr QFVertexP1(Geometry::Polytope::Type g)
        : Parent(g)
      {}

      size_t getSize() const override
      {
        return s_size[getGeometry()];
      }

      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        return s_points[getGeometry()][i];
      }

      Real getWeight(size_t i) const override
      {
        return s_weights[getGeometry()][i];
      }

      QFVertexP1* copy() const noexcept override
      {
        return new QFVertexP1(*this);
      }

    private:
      static const Geometry::GeometryIndexed<size_t> s_size;
      static const Geometry::GeometryIndexed<std::vector<Math::SpatialVector<Real>>> s_points;
      static const Geometry::GeometryIndexed<Math::Vector<Real>> s_weights;
  };
}

#endif
