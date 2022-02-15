#include "ScalarCoefficient.h"

#include "VectorDomainLFDivIntegrator.h"

namespace Rodin::Variational
{
   VectorDomainLFDivIntegrator
   ::VectorDomainLFDivIntegrator(const ScalarCoefficientBase& f)
      : m_f(f.copy())
   {}

   VectorDomainLFDivIntegrator
   ::VectorDomainLFDivIntegrator(const VectorDomainLFDivIntegrator& other)
      :  m_attr(other.m_attr),
         m_f(other.m_f->copy())
   {}

   void VectorDomainLFDivIntegrator::buildMFEMLinearFormIntegrator()
   {
      m_f->buildMFEMCoefficient();
      m_mfemLFI
         = std::make_unique<Internal::VectorDomainLFDivIntegrator>(
               m_f->getMFEMCoefficient());
   }

   mfem::LinearFormIntegrator&
   VectorDomainLFDivIntegrator::getMFEMLinearFormIntegrator()
   {
      assert(m_mfemLFI);
      return *m_mfemLFI;
   }

   mfem::LinearFormIntegrator*
   VectorDomainLFDivIntegrator::releaseMFEMLinearFormIntegrator()
   {
      return m_mfemLFI.release();
   }
}
