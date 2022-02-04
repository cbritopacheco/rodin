/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_VECTORMASSINTEGRATOR_H
#define RODIN_VARIATIONAL_VECTORMASSINTEGRATOR_H

#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    * Represents the integration of the linear form
    * @f[
    *    L(v) = \int_{\Omega} \lambda u \cdot v \ dx
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
   class VectorMassIntegrator : public BilinearFormDomainIntegrator
   {
      public:
         /**
          * @brief Constructs a VectorMassIntegrator with a scalar coefficient
          * @f$ \lambda = 1 @f$.
          */
         VectorMassIntegrator();

         /**
          * @brief Constructs a VectorMassIntegrator with scalar coefficient
          * @f$ \lambda @f$.
          *
          * @param[in] lambda Coefficient to integrate.
          */
         VectorMassIntegrator(const ScalarCoefficientBase& lambda);

         VectorMassIntegrator(const VectorMassIntegrator& other);

         const std::set<int>& getAttributes() const override
         {
            return m_attr;
         }

         VectorMassIntegrator& over(int attr) override
         {
            return over(std::set<int>{attr});
         }

         VectorMassIntegrator& over(const std::set<int>& attrs) override
         {
            m_attr = attrs;
            return *this;
         }

         void buildMFEMBilinearFormIntegrator() override;

         mfem::BilinearFormIntegrator& getMFEMBilinearFormIntegrator() override;
         mfem::BilinearFormIntegrator* releaseMFEMBilinearFormIntegrator() override;

         VectorMassIntegrator* copy() const noexcept override
         {
            return new VectorMassIntegrator(*this);
         }

      private:
         std::set<int> m_attr;
         std::unique_ptr<ScalarCoefficientBase> m_lambda;
         std::unique_ptr<mfem::VectorMassIntegrator> m_mfemBFI;
   };

}

#endif


