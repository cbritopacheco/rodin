#include "ScalarCoefficient.h"

#include "VectorDomainLFDivIntegrator.h"

namespace Rodin::Variational
{
   VectorDomainLFDivIntegrator
   ::VectorDomainLFDivIntegrator(const ScalarCoefficientBase& f)
      : m_f(f.copy()),
        m_mfemScalar(m_f->build()),
        m_mfemLFI(*m_mfemScalar)
   {}

   VectorDomainLFDivIntegrator
   ::VectorDomainLFDivIntegrator(const VectorDomainLFDivIntegrator& other)
      :  LinearFormDomainIntegrator(other),
         m_f(other.m_f->copy()),
         m_mfemScalar(m_f->build()),
         m_mfemLFI(*m_mfemScalar)
   {}
}
