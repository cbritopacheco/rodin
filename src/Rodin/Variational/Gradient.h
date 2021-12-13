#ifndef RODIN_VARIATIONAL_GRADIENT_H
#define RODIN_VARIATIONAL_GRADIENT_H

#include "VectorCoefficient.h"
#include "MatrixCoefficient.h"

namespace Rodin::Variational
{
   class GradientBase : public MatrixCoefficientBase
   {
      public:
         virtual void buildMFEMMatrixCoefficient() override = 0;
         virtual GradientBase* copy() const noexcept = 0;
   };

   class Gradient : public GradientBase
   {
      public:
         Gradient(GridFunctionBase& u)
            : m_u(u)
         {}

      private:
         GridFunctionBase& m_u;
   };
}

#endif
