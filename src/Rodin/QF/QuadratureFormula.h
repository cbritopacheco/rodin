#ifndef RODIN_VARIATIONAL_QUADRATUREFORMULA_H
#define RODIN_VARIATIONAL_QUADRATUREFORMULA_H

#include <vector>
#include <utility>

#include "Rodin/Math.h"
#include "Rodin/Copyable.h"
#include "Rodin/Geometry.h"

namespace Rodin::QF
{
  /**
   * @defgroup RodinQuadrature Quadrature formulae
   * @brief Quadrature formulae utilized by Rodin
   * @see QuadratureFormulaBase
   */

  /**
   * @brief Abstract base class for quadrature formulas.
   */
  class QuadratureFormulaBase : public Copyable
  {
    public:
      constexpr
      QuadratureFormulaBase(Geometry::Polytope::Type g)
        : m_geometry(g)
      {}

      constexpr
      QuadratureFormulaBase(const QuadratureFormulaBase& other)
        : m_geometry(other.m_geometry)
      {}

      virtual ~QuadratureFormulaBase() = default;

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
      virtual const Math::SpatialVector<Scalar>& getPoint(size_t i) const = 0;

      virtual QuadratureFormulaBase* copy() const noexcept override = 0;

    private:
      Geometry::Polytope::Type m_geometry;
  };
}

#endif

