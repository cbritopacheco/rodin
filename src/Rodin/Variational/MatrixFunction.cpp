#include "Rodin/Alert.h"

#include "Transpose.h"

#include "MatrixFunction.h"

namespace Rodin::Variational
{
   std::unique_ptr<mfem::MatrixCoefficient> MatrixFunctionBase::build() const
   {
      return std::unique_ptr<mfem::MatrixCoefficient>(new Internal::ProxyMatrixFunction(*this));
   }

   Transpose<MatrixFunctionBase> MatrixFunctionBase::T() const
   {
      return Transpose(*this);
   }
}

