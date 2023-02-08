#include "Rodin/Variational/FiniteElement.h"

#include "SimplexTransformation.h"

namespace Rodin::Geometry
{
  const Simplex& IsoparametricTransformation::getSimplex() const
  {
    return m_fe.get().getSimplex();
  }
}
