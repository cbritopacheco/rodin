#include "ScalarCoefficient.h"

#include "VectorBoundaryFluxLFIntegrator.h"

namespace Rodin::Variational
{
   VectorBoundaryFluxLFIntegrator
   ::VectorBoundaryFluxLFIntegrator(
         const std::vector<int>& bdrAttrs, const ScalarCoefficientBase& f)
      :  LinearFormBoundaryIntegrator(bdrAttrs),
         m_f(f.copy())
   {}

   VectorBoundaryFluxLFIntegrator
   ::VectorBoundaryFluxLFIntegrator(const ScalarCoefficientBase& f)
      : VectorBoundaryFluxLFIntegrator({}, f)
   {}

   VectorBoundaryFluxLFIntegrator
   ::VectorBoundaryFluxLFIntegrator(const VectorBoundaryFluxLFIntegrator& other)
      :  LinearFormBoundaryIntegrator(other),
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
