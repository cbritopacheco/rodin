/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "FormLanguage/ScalarUnaryMinus.h"

#include "Problem.h"

#include "ElasticityIntegrator.h"

namespace Rodin::Variational
{

   ElasticityIntegrator::ElasticityIntegrator(
         const ScalarCoefficientBase& lambda, const ScalarCoefficientBase& mu)
      : m_lambda(lambda.copy()), m_mu(mu.copy())
   {}

   ElasticityIntegrator::ElasticityIntegrator(const ElasticityIntegrator& other)
      :  m_attr(other.m_attr),
         m_lambda(other.m_lambda->copy()), m_mu(other.m_mu->copy())
   {}

   void ElasticityIntegrator::build()
   {
      m_lambda->build();
      m_mu->build();
      m_bfi = std::make_unique<mfem::ElasticityIntegrator>(m_lambda->get(), m_mu->get());
   }

   mfem::BilinearFormIntegrator&
   ElasticityIntegrator::get()
   {
      assert(m_bfi);
      return *m_bfi;
   }

   mfem::BilinearFormIntegrator*
   ElasticityIntegrator::release()
   {
      return m_bfi.release();
   }
}

