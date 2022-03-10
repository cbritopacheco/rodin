#include "Rodin/Alert.h"

#include "Transpose.h"

#include "MatrixCoefficient.h"

namespace Rodin::Variational
{
   std::unique_ptr<mfem::MatrixCoefficient> MatrixCoefficientBase::build() const
   {
      return std::unique_ptr<mfem::MatrixCoefficient>(new Internal::ProxyMatrixCoefficient(*this));
   }

   Transpose<MatrixCoefficientBase> MatrixCoefficientBase::T() const
   {
      return Transpose(*this);
   }
}

