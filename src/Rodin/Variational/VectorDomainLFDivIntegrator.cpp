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

   void VectorDomainLFDivIntegrator::build()
   {
      m_f->build();
      m_mfemLFI
         = std::make_unique<Internal::VectorDomainLFDivIntegrator>(
               m_f->get());
   }

   mfem::LinearFormIntegrator&
   VectorDomainLFDivIntegrator::get()
   {
      assert(m_mfemLFI);
      return *m_mfemLFI;
   }

   mfem::LinearFormIntegrator*
   VectorDomainLFDivIntegrator::release()
   {
      return m_mfemLFI.release();
   }
}
