/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MASSINTEGRATOR_H
#define RODIN_VARIATIONAL_MASSINTEGRATOR_H

#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    * Represents the integration of the linear form
    * @f[
    *    L(v) = \int_{\Omega} \lambda u v \ dx
    * @f]
    * where @f$ \lambda @f$ is a scalar coefficient.
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | L2, H1                                       |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  Continuous operator  | @f$ f @f$                                    |
    * |  @f$ \lambda @f$      | ScalarCoefficient                            |
    *
    */
   class MassIntegrator : public BilinearFormDomainIntegrator
   {
      public:
         /**
          * @brief Constructs a MassIntegrator with a scalar coefficient
          * @f$ \lambda = 1 @f$.
          */
         MassIntegrator();

         /**
          * @brief Constructs a MassIntegrator with scalar coefficient
          * @f$ \lambda @f$.
          *
          * @param[in] lambda Coefficient to integrate.
          */
         MassIntegrator(const ScalarCoefficientBase& lambda);

         MassIntegrator(const MassIntegrator& other);

         void buildMFEMBilinearFormIntegrator() override;

         mfem::BilinearFormIntegrator& getMFEMBilinearFormIntegrator() override;
         mfem::BilinearFormIntegrator* releaseMFEMBilinearFormIntegrator() override;

         MassIntegrator* copy() const noexcept override
         {
            return new MassIntegrator(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_lambda;
         std::unique_ptr<mfem::MassIntegrator> m_mfemBFI;
   };

}

#endif
