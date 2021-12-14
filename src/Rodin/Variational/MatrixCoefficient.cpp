#include "Transpose.h"

#include "MatrixCoefficient.h"

namespace Rodin::Variational
{
   Transpose MatrixCoefficientBase::T() const
   {
      return Transpose(*this);
   }
}

