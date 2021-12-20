/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DOMAINLFINTEGRATOR_H
#define RODIN_VARIATIONAL_DOMAINLFINTEGRATOR_H

#include "ForwardDecls.h"
#include "ScalarCoefficient.h"
#include "LinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    *
    * @brief Object used to integrate a coefficient against a fixed coefficient.
    *
    * Represents the integration of the linear form
    * @f[
    *    L(v) = \int_{\Omega} f v \ dx
    * @f]
    * where @f$ f @f$ is a scalar coefficient.
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | L2, H1                                       |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  Continuous operator  | @f$ f @f$                                    |
    * |  @f$ f @f$            | ScalarCoefficient                            |
    *
    */
   class DomainLFIntegrator : public LinearFormIntegratorBase
   {
      public:
         /**
          * @brief Constructs a DomainLFIntegrator with scalar coefficient
          * @f$ f @f$.
          *
          * @param[in] f Coefficient to integrate.
          */
         DomainLFIntegrator(const ScalarCoefficientBase& f);

         DomainLFIntegrator(const DomainLFIntegrator& other);

         void buildMFEMLinearFormIntegrator() override;

         mfem::LinearFormIntegrator& getMFEMLinearFormIntegrator() override;
         mfem::LinearFormIntegrator* releaseMFEMLinearFormIntegrator() override;

         DomainLFIntegrator* copy() const noexcept override;

      private:
         std::unique_ptr<ScalarCoefficientBase> m_f;
         std::unique_ptr<mfem::DomainLFIntegrator> m_mfemLFI;
   };
}

#endif

