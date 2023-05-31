#ifndef RODIN_VARIATIONAL_QFGG_H
#define RODIN_VARIATIONAL_QFGG_H

#define RODIN_VARIATIONAL_QFGG_MAX_ORDER 256

#include "QuadratureFormula.h"

namespace Rodin::Variational
{
  class QFGG final : public QuadratureFormulaBase
  {
    public:
      constexpr
      QFGG(Geometry::Polytope::Geometry g, size_t order)
        : m_geometry(g), m_order(order)
      {}

      inline
      const Math::Matrix& getPoints() const override
      {
        return s_points[m_order][m_geometry];
      }

      inline
      const Math::Vector& getWeights() const override
      {
        return s_weights[m_order][m_geometry];
      }

      inline
      Geometry::Polytope::Geometry getGeometry() const override
      {
        return m_geometry;
      }

    private:
      static std::array<Geometry::GeometryIndexed<Math::Matrix>, RODIN_VARIATIONAL_QFGG_MAX_ORDER> s_points;
      static std::array<Geometry::GeometryIndexed<Math::Vector>, RODIN_VARIATIONAL_QFGG_MAX_ORDER> s_weights;

      Geometry::Polytope::Geometry m_geometry;
      const size_t m_order;
  };
}

#endif

