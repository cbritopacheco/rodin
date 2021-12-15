/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRADIENT_H
#define RODIN_VARIATIONAL_GRADIENT_H

#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   class Gradient : public VectorCoefficientBase
   {
      public:
         Gradient(GridFunction<H1>& u);

         size_t getDimension() const override;

         void buildMFEMVectorCoefficient() override;

         mfem::VectorCoefficient& getMFEMVectorCoefficient() override;

         VectorCoefficientBase* copy() const noexcept override
         {
            return new Gradient(*this);
         }

      private:
         GridFunction<H1>& m_u;
         std::optional<mfem::GradientGridFunctionCoefficient> m_mfemVectorCoefficient;
   };
}

#endif
