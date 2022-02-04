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

   void ElasticityIntegrator::buildMFEMBilinearFormIntegrator()
   {
      m_lambda->buildMFEMCoefficient();
      m_mu->buildMFEMCoefficient();
      m_bfi = std::make_unique<mfem::ElasticityIntegrator>(
            m_lambda->getMFEMCoefficient(), m_mu->getMFEMCoefficient());
   }

   mfem::BilinearFormIntegrator&
   ElasticityIntegrator::getMFEMBilinearFormIntegrator()
   {
      assert(m_bfi);
      return *m_bfi;
   }

   mfem::BilinearFormIntegrator*
   ElasticityIntegrator::releaseMFEMBilinearFormIntegrator()
   {
      return m_bfi.release();
   }
}

