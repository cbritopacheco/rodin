#ifndef RODIN_VARIATIONAL_QF_QF2P1_H
#define RODIN_VARIATIONAL_QF_QF2P1_H

#include "Rodin/Geometry/GeometryIndexed.h"

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   */
  class QF2P1 final : public QuadratureFormulaBase
  {
    public:
      using Parent = QuadratureFormulaBase;

      constexpr
      QF2P1(Geometry::Polytope::Type g)
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

      QF2P1* copy() const noexcept override
      {
        return new QF2P1(*this);
      }

    private:
      static const Geometry::GeometryIndexed<size_t> s_size;
      static const Geometry::GeometryIndexed<std::vector<Math::SpatialVector<Real>>> s_points;
      static const Geometry::GeometryIndexed<Math::Vector<Real>> s_weights;
  };
}

#endif

