#include "ScalarCoefficient.h"

#include "MassIntegrator.h"

namespace Rodin::Variational
{
   MassIntegrator::MassIntegrator()
      : MassIntegrator(ScalarCoefficient(1.0))
   {}

   MassIntegrator::MassIntegrator(const ScalarCoefficientBase& lambda)
      : m_lambda(lambda.copy()),
        m_mfemLambda(m_lambda->build()),
        m_mfemBFI(*m_mfemLambda)
   {}

   MassIntegrator::MassIntegrator(const MassIntegrator& other)
      :  m_attr(other.m_attr),
         m_lambda(other.m_lambda->copy()),
         m_mfemLambda(m_lambda->build()),
         m_mfemBFI(*m_mfemLambda)
   {}
}

