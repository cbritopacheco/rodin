/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "FormLanguage/ScalarUnaryMinus.h"

#include "Problem.h"

#include "VectorDiffusionIntegrator.h"

namespace Rodin::Variational
{
   VectorDiffusionIntegrator::VectorDiffusionIntegrator()
      : VectorDiffusionIntegrator(ScalarCoefficient(1.0))
   {}

   VectorDiffusionIntegrator::VectorDiffusionIntegrator(const ScalarCoefficientBase& lambda)
      : m_lambda(lambda.copy()),
        m_mfemLambda(m_lambda->build()),
        m_bfi(*m_mfemLambda)
   {}

   VectorDiffusionIntegrator::VectorDiffusionIntegrator(const VectorDiffusionIntegrator& other)
      :  m_attr(other.m_attr),
         m_lambda(other.m_lambda->copy()),
         m_mfemLambda(m_lambda->build()),
         m_bfi(*m_mfemLambda)
   {}
}


