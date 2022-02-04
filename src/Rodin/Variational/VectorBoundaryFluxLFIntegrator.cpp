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
   ::buildMFEMLinearFormIntegrator()
   {
      m_f->buildMFEMCoefficient();
      m_mfemLFI
         = std::make_unique<mfem::VectorBoundaryFluxLFIntegrator>(
               m_f->getMFEMCoefficient());
   }

   mfem::LinearFormIntegrator&
   VectorBoundaryFluxLFIntegrator::getMFEMLinearFormIntegrator()
   {
      assert(m_mfemLFI);
      return *m_mfemLFI;
   }

   mfem::LinearFormIntegrator*
   VectorBoundaryFluxLFIntegrator::releaseMFEMLinearFormIntegrator()
   {
      return m_mfemLFI.release();
   }
}
