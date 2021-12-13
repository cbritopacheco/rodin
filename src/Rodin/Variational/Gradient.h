/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
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
