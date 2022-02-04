#include "ScalarCoefficient.h"

#include "DivergenceIntegrator.h"

namespace Rodin::Variational
{
   DivergenceIntegrator<Domain, Linear>
   ::DivergenceIntegrator(const ScalarCoefficientBase& f)
      : m_f(f.copy())
   {}

   DivergenceIntegrator<Domain, Linear>
   ::DivergenceIntegrator(const DivergenceIntegrator& other)
      :  m_attr(other.m_attr),
         m_f(other.m_f->copy())
   {}

   void DivergenceIntegrator<Domain, Linear>::buildMFEMLinearFormIntegrator()
   {
      m_f->buildMFEMCoefficient();
      m_mfemLFI
         = std::make_unique<Internal::VectorDomainLFDivIntegrator>(
               m_f->getMFEMCoefficient());
   }

   mfem::LinearFormIntegrator&
   DivergenceIntegrator<Domain, Linear>::getMFEMLinearFormIntegrator()
   {
      assert(m_mfemLFI);
      return *m_mfemLFI;
   }

   mfem::LinearFormIntegrator*
   DivergenceIntegrator<Domain, Linear>::releaseMFEMLinearFormIntegrator()
   {
      return m_mfemLFI.release();
   }
}
