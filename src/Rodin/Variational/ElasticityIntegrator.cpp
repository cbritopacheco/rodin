/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Problem.h"

#include "ElasticityIntegrator.h"

namespace Rodin::Variational
{
   ElasticityIntegrator::ElasticityIntegrator(const ElasticityIntegrator& other)
      : m_lambda(other.m_lambda->copy()), m_mu(other.m_mu->copy()), m_bf(other.m_bf)
   {}

   ElasticityIntegrator&
   ElasticityIntegrator::setBilinearForm(BilinearFormBase& bf)
   {
      m_bf.emplace(bf);
      return *this;
   }

   void ElasticityIntegrator::eval()
   {
      m_lambda->buildMFEMCoefficient();
      m_mu->buildMFEMCoefficient();

      m_bf->get()
         .getHandle()
         .AddDomainIntegrator(
               new mfem::ElasticityIntegrator(
                  m_lambda->getMFEMCoefficient(),
                  m_mu->getMFEMCoefficient()));
   }

   ElasticityIntegrator& ElasticityIntegrator::toggleSign()
   {
      m_lambda.reset(new ScalarCoefficient(-(*m_lambda)));
      m_mu.reset(new ScalarCoefficient(-(*m_mu)));
      return *this;
   }

   ElasticityIntegrator* ElasticityIntegrator::copy() const noexcept
   {
      return new ElasticityIntegrator(*this);
   }
}

