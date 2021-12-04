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
    * @brief Object used to integrate the typical elasticity operator.
    *
    * Represents the integration of the bilinear forms
    *
    * @f$ \int_{\Omega} c_{ijkl} \nabla u_j \cdot \nabla v_i \ dx @f$
    *
    * where:
    * @f$ c_{ijkl} := \lambda \delta_{ik} \delta_{jl} + \mu (\delta_{ij}
    * \delta_{kl} + \delta_{il} \delta_{jk}) @f$
    *
    * for @f$ 1 \leq i, j, k, l \leq d @f$.
    *
    *
    * ----
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  Continuous operator  | @f$ - \nabla \cdot \sigma (u) @f$            |
    * |  @f$ \mu @f$          | Scalar                                       |
    * |  @f$ \lambda @f$      | Scalar                                       |
    *
    * ----
    *
    * where:
    *
    * @f[
    *   \sigma(u) := \lambda (\nabla \cdot u) I + \mu (\nabla u + {\nabla u}^T)
    * @f]
    *
    */
   template <class L, class M>
   class ElasticityIntegrator
      :  public FormLanguage::BilinearFormExpr<ElasticityIntegrator<L, M>>
   {
      public:
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
