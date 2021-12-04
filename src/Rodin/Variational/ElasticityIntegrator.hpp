/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ELASTICITYINTEGRATOR_HPP
#define RODIN_VARIATIONAL_ELASTICITYINTEGRATOR_HPP

#include "Problem.h"

#include "ElasticityIntegrator.h"

namespace Rodin::Variational
{
   template <class L, class M>
   ElasticityIntegrator<L, M>::ElasticityIntegrator(
         const L& lambda, const M& mu)
      : m_lambda(ScalarCoefficient<L>(lambda)), m_mu(ScalarCoefficient<M>(mu))
   {}

   template <class L, class M>
   ElasticityIntegrator<L, M>::ElasticityIntegrator(const ElasticityIntegrator& other)
      : m_lambda(other.m_lambda), m_mu(other.m_mu), m_bf(other.m_bf)
   {}

   template <class L, class M>
   ElasticityIntegrator<L, M>&
   ElasticityIntegrator<L, M>::setBilinearForm(BilinearFormBase& bf)
   {
      m_bf.emplace(bf);
      return *this;
   }

   template <class L, class M>
   void ElasticityIntegrator<L, M>::eval()
   {
      m_lambda.eval();
      m_mu.eval();

      m_bf->get()
         .getHandle()
         .AddDomainIntegrator(
               new mfem::ElasticityIntegrator(m_lambda.coeff(), m_mu.coeff()));
   }

   template <class L, class M>
   ElasticityIntegrator<L, M>& ElasticityIntegrator<L, M>::toggleSign()
   {
      m_lambda.toggleSign();
      m_mu.toggleSign();
      return *this;
   }

   template <class L, class M>
   template <class ... Args>
   ElasticityIntegrator<L, M>*
   ElasticityIntegrator<L, M>::create(Args&&... args) noexcept
   {
      return new ElasticityIntegrator(std::forward<Args>(args)...);
   }

   template <class L, class M>
   ElasticityIntegrator<L, M>*
   ElasticityIntegrator<L, M>::copy() const noexcept
   {
      return new ElasticityIntegrator(*this);
   }
}

#endif
