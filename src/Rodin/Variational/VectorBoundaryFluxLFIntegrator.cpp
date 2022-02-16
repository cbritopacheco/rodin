#include "ScalarCoefficient.h"

#include "VectorBoundaryFluxLFIntegrator.h"

namespace Rodin::Variational
{
   VectorBoundaryFluxLFIntegrator::VectorBoundaryFluxLFIntegrator(
         const ScalarCoefficientBase& f)
      :  m_f(f.copy())
   {}

   VectorBoundaryFluxLFIntegrator
   ::VectorBoundaryFluxLFIntegrator(const VectorBoundaryFluxLFIntegrator& other)
      : m_attr(other.m_attr),
        m_f(other.m_f->copy())
   {}

   void VectorBoundaryFluxLFIntegrator
   ::build()
   {
      m_f->build();
      m_mfemLFI
         = std::make_unique<mfem::VectorBoundaryFluxLFIntegrator>(
               m_f->get());
   }

   mfem::LinearFormIntegrator&
   VectorBoundaryFluxLFIntegrator::get()
   {
      assert(m_mfemLFI);
      return *m_mfemLFI;
   }

   mfem::LinearFormIntegrator*
   VectorBoundaryFluxLFIntegrator::release()
   {
      return m_mfemLFI.release();
   }
}
