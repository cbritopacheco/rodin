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

namespace Rodin::Variational
{
   /**
    *
    * @brief Object used to integrate the linear elasticity form.
    *
    * Represents the integration of the bilinear form
    *
    * @f[
    *    a(\vec{u}, \vec{v}) = \int_\Omega A e(\vec{u}) \colon e(\vec{v}) \ dx
    *    .
    * @f]
    * The stress tensor @f$ \sigma(\vec{u}) @f$ is related to the strain tensor
    * @f$ e(\vec{u}) := \frac{1}{2} \left( \nabla \vec{u} + \nabla \vec{u}^T
    * \right) @f$ via Hooke's law:
    * @f[
    *    \sigma(\vec{u}) = A e(\vec{u}),
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
    * |  Continuous operator  | @f$ - \nabla \cdot \sigma (\vec{u}) @f$      |
    * |  @f$ \mu @f$          | Scalar                                       |
    * |  @f$ \lambda @f$      | Scalar                                       |
    *
    * ----
    *
    */
   template <class L, class M>
   class ElasticityIntegrator
      :  public FormLanguage::BilinearFormExpr<ElasticityIntegrator<L, M>>
   {
      public:
         /**
          * @brief Creates an ElasticityIntegrator with the given Lamé
          * coefficients @f$ \lambda @f$ and @f$ \mu @f$.
          */
         ElasticityIntegrator(const L& lambda, const M& mu);

         ElasticityIntegrator(const ElasticityIntegrator& other);

         virtual ElasticityIntegrator& setBilinearForm(BilinearFormBase& bf) override;

         void eval() override;

         ElasticityIntegrator& toggleSign() override;

         template <class ... Args>
         static ElasticityIntegrator* create(Args&&... args) noexcept;

         virtual ElasticityIntegrator* copy() const noexcept override;

      private:
         ScalarCoefficient<L> m_lambda;
         ScalarCoefficient<M> m_mu;

         std::optional<std::reference_wrapper<BilinearFormBase>> m_bf;
   };
}

#include "ElasticityIntegrator.hpp"

#endif
