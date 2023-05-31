#ifndef RODIN_VARIATIONAL_QUADRATUREFORMULA_H
#define RODIN_VARIATIONAL_QUADRATUREFORMULA_H

#include <vector>
#include <utility>

#include "Rodin/Math.h"
#include "Rodin/Geometry.h"

namespace Rodin::Variational
{
  class QuadratureFormulaBase
  {
    public:
      inline
      size_t getSize() const
      {
        assert(getPoints().cols() == getWeights().size());
        return getWeights().size();
      }

      inline
      auto getPoint(size_t i) const
      {
        return getPoints().col(i);
      }

      inline
      Scalar getWeight(size_t i) const
      {
        return getWeights().coeff(i);
      }

      virtual const Math::Matrix& getPoints() const = 0;

      virtual const Math::Vector& getWeights() const = 0;

      virtual Geometry::Polytope::Geometry getGeometry() const = 0;
  };
}

#endif

