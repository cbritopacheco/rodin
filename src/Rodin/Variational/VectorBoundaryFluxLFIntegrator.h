/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_VECTORBOUNDARYFLUXLFINTEGRATOR_H
#define RODIN_VARIATIONAL_VECTORBOUNDARYFLUXLFINTEGRATOR_H

#include "LinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    * Represents the integration of the linear form
    * @f[
    *    L(v) = \int_{\Gamma} \lambda v \cdot n \ dx
    * @f]
    * where @f$ \lambda @f$ is a scalar coefficient and @f$ \Gamma @f$ is a
    * segment of the boundary with normal vector @f$ n @f$.
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | L2, H1                                       |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  Continuous operator  | @f$ f @f$                                    |
    * |  @f$ \lambda @f$            | ScalarCoefficient                      |
    *
    */
   class VectorBoundaryFluxLFIntegrator : public LinearFormBoundaryIntegrator
   {
      public:
         /**
          * @brief Constructs a VectorBoundaryFluxLFIntegrator with scalar coefficient
          * @f$ \lambda @f$.
          *
          * @param[in] lambda Coefficient to integrate.
          */
         VectorBoundaryFluxLFIntegrator(const ScalarCoefficientBase& lambda);

         VectorBoundaryFluxLFIntegrator(
               int bdrAttr,
               const ScalarCoefficientBase& lambda);

         VectorBoundaryFluxLFIntegrator(
               const std::vector<int>& bdrAttrs,
               const ScalarCoefficientBase& lambda);

         VectorBoundaryFluxLFIntegrator(const VectorBoundaryFluxLFIntegrator& other);

         void buildMFEMLinearFormIntegrator() override;

         mfem::LinearFormIntegrator& getMFEMLinearFormIntegrator() override;
         mfem::LinearFormIntegrator* releaseMFEMLinearFormIntegrator() override;

         VectorBoundaryFluxLFIntegrator* copy() const noexcept override
         {
            return new VectorBoundaryFluxLFIntegrator(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_f;
         std::unique_ptr<mfem::VectorBoundaryFluxLFIntegrator> m_mfemLFI;
   };

}

#endif