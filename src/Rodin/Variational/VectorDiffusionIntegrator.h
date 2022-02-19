/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_VECTORDIFFUSIONINTEGRATOR_H
#define RODIN_VARIATIONAL_VECTORDIFFUSIONINTEGRATOR_H

#include <optional>
#include <functional>

#include "ForwardDecls.h"
#include "ScalarCoefficient.h"

#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{

   /**
    * @brief Object used to integrate the typical diffusion operator.
    *
    * Represents the integration of the bilinear form
    * @f[
    *    a(u, v) = \int_\Omega \lambda \nabla u \colon \nabla v \ dx
    * @f]
    * where @f$ \lambda @f$ is a scalar coefficient and @f$ \nabla u @f$ is the
    * Jacobian matrix of @f$ u @f$.
    *
    * ----
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | H1                                           |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  @f$ \lambda @f$      | ScalarCoefficient                            |
    *
    * ----
    *
    */
   class VectorDiffusionIntegrator : public BilinearFormDomainIntegrator
   {
      public:
         /**
          * @brief Constructs a VectorDiffusionIntegrator with a scalar coefficient
          * @f$ \lambda = 1@f$.
          */
         VectorDiffusionIntegrator();

         /**
          * @brief Constructs a VectorDiffusionIntegrator with scalar coefficient
          * @f$ \lambda @f$.
          *
          * @param[in] lambda VectorDiffusion coefficient
          */
         VectorDiffusionIntegrator(const ScalarCoefficientBase& lambda);

         VectorDiffusionIntegrator(const VectorDiffusionIntegrator& other);

         const std::set<int>& getAttributes() const override
         {
            return m_attr;
         }

         VectorDiffusionIntegrator& over(int attr) override
         {
            return over(std::set<int>{attr});
         }

         VectorDiffusionIntegrator& over(const std::set<int>& attrs) override
         {
            m_attr = attrs;
            return *this;
         }

         void getElementMatrix(
               const mfem::FiniteElement& trial,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override
         {
            m_bfi.AssembleElementMatrix(trial, trans, mat);
         }

         VectorDiffusionIntegrator* copy() const noexcept override
         {
            return new VectorDiffusionIntegrator(*this);
         }

      private:
         std::set<int> m_attr;
         std::unique_ptr<ScalarCoefficientBase> m_lambda;

         std::unique_ptr<Internal::ScalarCoefficient> m_mfemLambda;
         mfem::VectorDiffusionIntegrator m_bfi;
   };
}

#endif

