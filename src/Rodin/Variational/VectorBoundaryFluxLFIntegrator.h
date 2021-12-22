#ifndef RODIN_VARIATIONAL_VECTORBOUNDARYFLUXLFINTEGRATOR_H
#define RODIN_VARIATIONAL_VECTORBOUNDARYFLUXLFINTEGRATOR_H

#include "LinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    *
    * @brief 
    *
    * Represents the integration of the linear form
    * @f[
    *    L(v) = \int_{\Gamma} f v \cdot n \ dx
    * @f]
    * where @f$ f @f$ is a scalar coefficient and @f$ \Gamma @f$ is a segment
    * of the boundary.
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | L2, H1                                       |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  Continuous operator  | @f$ f @f$                                    |
    * |  @f$ f @f$            | ScalarCoefficient                            |
    *
    */
   class VectorBoundaryFluxLFIntegrator : public LinearFormIntegratorBase
   {
      public:
         /**
          * @brief Constructs a VectorBoundaryFluxLFIntegrator with scalar coefficient
          * @f$ f @f$.
          *
          * @param[in] f Coefficient to integrate.
          */
         VectorBoundaryFluxLFIntegrator(const ScalarCoefficientBase& f);

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
