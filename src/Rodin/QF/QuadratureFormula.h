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
      QuadratureFormulaBase(Geometry::Polytope::Type g)
        : m_geometry(g)
      {}

      inline
      constexpr
      Geometry::Polytope::Type getGeometry() const
      {
        return m_geometry;
      }

      virtual size_t getSize() const = 0;

      virtual Scalar getWeight(size_t i) const = 0;

      /**
       * @brief Returns a reference to a vector containing the coordinates in
       * reference space.
       *
       * @note The reference must be valid throughout the whole lifetime of the
       * program.
       */
      virtual const Math::SpatialVector& getPoint(size_t i) const = 0;

    private:
      Geometry::Polytope::Type m_geometry;
  };
}

#endif
