/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ELASTICITYINTEGRATOR_H
#define RODIN_VARIATIONAL_ELASTICITYINTEGRATOR_H

#include <functional>

#include "ForwardDecls.h"
#include "ScalarCoefficient.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    *
    * @brief Object used to integrate the linear elasticity form.
    *
    * Represents the integration of the bilinear form
    *
    * @f[
    *    a(u, v) = \int_\Omega A e(u) \colon e(v) \ dx
    *    .
    * @f]
    * The stress tensor @f$ \sigma(u) @f$ is related to the strain tensor
    * @f$ e(u) := \frac{1}{2} \left( \nabla u + \nabla u^T \right) @f$ via
    * Hooke's law:
    * @f[
    *    \sigma(u) = A e(u),
    * @f]
    * where for any @f$ e @f$ in the set of real symmetric @f$ d \times d @f$
    * matrices,
    * @f[
    *   Ae = 2\mu e + \lambda \mathrm{tr}(e)I,
    * @f]
    * @f$ I @f$ is the identity @f$ d \times d @f$ matrix, @f$ \lambda @f$ and
    * @f$ \mu @f$ are the Lamé parameters of the constituent material,
    * satisfying @f$ \mu > 0 @f$ and @f$ \lambda + 2 \mu / d > 0 @f$.
    *
    * ----
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  Continuous operator  | @f$ - \nabla \cdot \sigma (u) @f$            |
    * |  @f$ \mu @f$          | ScalarCoefficient                            |
    * |  @f$ \lambda @f$      | ScalarCoefficient                            |
    *
    * ----
    *
    */
   class ElasticityIntegrator : public BilinearFormDomainIntegrator
   {
      public:
         /**
          * @brief Creates an ElasticityIntegrator with the given Lamé
          * coefficients @f$ \lambda @f$ and @f$ \mu @f$.
          */
         ElasticityIntegrator(
            const ScalarCoefficientBase& lambda, const ScalarCoefficientBase& mu);

         ElasticityIntegrator(const ElasticityIntegrator& other);

         void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override
         {
            assert(&trial == &test);
            m_bfi.AssembleElementMatrix(trial, trans, mat);
         }

         ElasticityIntegrator* copy() const noexcept override
         {
            return new ElasticityIntegrator(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_lambda;
         std::unique_ptr<ScalarCoefficientBase> m_mu;

         std::unique_ptr<Internal::ScalarCoefficient> m_mfemLambda;
         std::unique_ptr<Internal::ScalarCoefficient> m_mfemMu;
         mfem::ElasticityIntegrator m_bfi;
   };
}

#endif
