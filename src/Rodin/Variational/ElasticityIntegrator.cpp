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
      : m_lambda(lambda.copy()), m_mu(mu.copy()),
        m_mfemLambda(m_lambda->build()), m_mfemMu(m_mu->build()),
        m_bfi(*m_mfemLambda, *m_mfemMu)
   {}

   ElasticityIntegrator::ElasticityIntegrator(const ElasticityIntegrator& other)
      :  BilinearFormDomainIntegrator(other),
         m_lambda(other.m_lambda->copy()), m_mu(other.m_mu->copy()),
         m_mfemLambda(m_lambda->build()), m_mfemMu(m_mu->build()),
         m_bfi(*m_mfemLambda, *m_mfemMu)
   {}
}

