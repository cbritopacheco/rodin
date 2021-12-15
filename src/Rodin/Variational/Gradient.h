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
   /**
    * @brief Represents the gradient @f$ \nabla u @f$ of the scalar function
    * @f$ u @f$.
    *
    * For @f$ u : \mathbb{R}^n \rightarrow \mathbb{R} @f$, the gradient
    * @f$ \nabla u : \mathbb{R}^n \rightarrow \mathbb{R} @f$ at the point
    * @f$ x = (x_1, \ldots, x_n) @f$ is defined by
    * @f[
    *    \nabla u (x) =
    *    \left[
    *       \dfrac{\partial u}{\partial x_1}(x), \ldots,
    *       \dfrac{\partial u}{\partial x_n}(x)
    *    \right]^T
    * @f]
    */
   class Gradient : public VectorCoefficientBase
   {
      public:
         /**
          * @brief Constructs the Gradient of an @f$ H^1 @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
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
