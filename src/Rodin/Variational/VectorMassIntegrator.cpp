#include "ScalarCoefficient.h"

#include "VectorMassIntegrator.h"

namespace Rodin::Variational
{
   VectorMassIntegrator::VectorMassIntegrator()
      : VectorMassIntegrator(ScalarCoefficient(1.0))
   {}

   VectorMassIntegrator::VectorMassIntegrator(const ScalarCoefficientBase& lambda)
      : m_lambda(lambda.copy()),
        m_mfemLambda(lambda.build()),
        m_mfemBFI(*m_mfemLambda)
   {}

   VectorMassIntegrator::VectorMassIntegrator(const VectorMassIntegrator& other)
      :  m_attr(other.m_attr),
         m_lambda(other.m_lambda->copy()),
         m_mfemLambda(m_lambda->build()),
         m_mfemBFI(*m_mfemLambda)
   {}
}


