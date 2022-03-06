/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_VECTORDOMAINLFINTEGRATOR_H
#define RODIN_VARIATIONAL_VECTORDOMAINLFINTEGRATOR_H

#include "ForwardDecls.h"
#include "VectorCoefficient.h"
#include "LinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    *
    * @brief Object used to integrate a coefficient against a fixed coefficient.
    *
    * Represents the integration of the linear form
    * @f[
    *    L(v) = \int_{\Omega} f \cdot v \ dx
    * @f]
    * where @f$ f @f$ is a vector coefficient.
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | L2, H1                                       |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  Continuous operator  | @f$ f @f$                                    |
    * |  @f$ f @f$            | VectorCoefficient                            |
    *
    */
   class VectorDomainLFIntegrator : public LinearFormDomainIntegrator
   {
      public:
         /**
          * @brief Constructs a VectorDomainLFIntegrator with scalar coefficient
          * @f$ f @f$.
          *
          * @param[in] f Coefficient to integrate.
          */
         VectorDomainLFIntegrator(const VectorCoefficientBase& f)
            : m_f(f.copy()),
              m_mfemVector(f.build()),
              m_mfemLFI(*m_mfemVector)
         {}

         VectorDomainLFIntegrator(const VectorDomainLFIntegrator& other)
            :  LinearFormDomainIntegrator(other),
               m_f(other.m_f->copy()),
               m_mfemVector(other.m_f->build()),
               m_mfemLFI(*m_mfemVector)
         {}

         void getElementVector(
               const mfem::FiniteElement& fe, mfem::ElementTransformation&
               trans, mfem::Vector& vec) override
         {
            m_mfemLFI.AssembleRHSElementVect(fe, trans, vec);
         }

         VectorDomainLFIntegrator* copy() const noexcept override
         {
            return new VectorDomainLFIntegrator(*this);
         }

      private:
         std::unique_ptr<VectorCoefficientBase> m_f;
         std::unique_ptr<mfem::VectorCoefficient> m_mfemVector;
         mfem::VectorDomainLFIntegrator m_mfemLFI;
   };
}

#endif


