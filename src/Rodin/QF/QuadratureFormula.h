#ifndef RODIN_VARIATIONAL_QUADRATUREFORMULA_H
#define RODIN_VARIATIONAL_QUADRATUREFORMULA_H

#include <vector>
#include <utility>

#include "Rodin/Math.h"
#include "Rodin/Geometry.h"

namespace Rodin::QF
{
  /**
   * @brief Abstract base class for quadrature formulas.
   */
  class QuadratureFormulaBase
  {
    public:
      constexpr
      QuadratureFormulaBase(Geometry::Polytope::Geometry g)
        : m_geometry(g)
      {}

      inline
      constexpr
      Geometry::Polytope::Geometry getGeometry() const
      {
        return m_geometry;
      }

      virtual size_t getSize() const = 0;

      virtual Scalar getWeight(size_t i) const = 0;

      virtual const Math::Vector& getPoint(size_t i) const = 0;

    private:
      Geometry::Polytope::Geometry m_geometry;
  };
}

#endif

