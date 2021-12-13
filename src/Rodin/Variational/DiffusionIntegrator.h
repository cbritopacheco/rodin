/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIFFUSIONINTEGRATOR_H
#define RODIN_VARIATIONAL_DIFFUSIONINTEGRATOR_H

#include <optional>
#include <functional>

#include "ForwardDecls.h"
#include "ScalarCoefficient.h"

#include "FormLanguage/BilinearFormExpr.h"

namespace Rodin::Variational
{

   /**
    * @brief Object used to integrate the typical diffusion operator.
    *
    * Represents the integration of the bilinear form
    * @f[
    *    a(u, v) = \int_\Omega \lambda \nabla u \cdot \nabla v \ dx
    * @f]
    * where @f$ \lambda @f$ is a scalar coefficient.
    *
    * ----
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | H1                                           |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  Continuous operator  | @f$ - \nabla \cdot (\lambda \nabla u) @f$    |
    * |  @f$ \lambda @f$      | ScalarCoefficient                            |
    *
    * ----
    *
    */
   class DiffusionIntegrator
      :  public FormLanguage::BilinearFormExpr<DiffusionIntegrator>
   {
      public:
         /**
          * @brief Constructs a DiffusionIntegrator with a scalar coefficient
          * @f$ \lambda = 1@f$.
          */
         DiffusionIntegrator()
            : m_lambda(new ScalarCoefficient(1.0))
         {}

         /**
          * @brief Constructs a DiffusionIntegrator with a scalar coefficient
          * @f$ \lambda @f$.
          *
          * @param[in] lambda Diffusion coefficient
          */
         template <class T>
         DiffusionIntegrator(const ScalarCoefficient<T>& lambda)
            : m_lambda(lambda)
         {}

         DiffusionIntegrator(const DiffusionIntegrator& other);

         virtual DiffusionIntegrator& setBilinearForm(BilinearFormBase& bf) override;

         void eval() override;

         DiffusionIntegrator& toggleSign() override;

         virtual DiffusionIntegrator* copy() const noexcept override
         {
            return new DiffusionIntegrator(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_lambda;
         std::optional<std::reference_wrapper<BilinearFormBase>> m_bf;
   };
}

#endif
