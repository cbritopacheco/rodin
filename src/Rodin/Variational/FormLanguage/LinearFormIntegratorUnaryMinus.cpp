#include "LinearFormIntegratorUnaryMinus.h"

namespace Rodin::Variational::FormLanguage
{
   LinearFormIntegratorUnaryMinus
   ::LinearFormIntegratorUnaryMinus(const LinearFormIntegratorBase& lfi)
      : m_lfi(lfi.copy())
   {}

   LinearFormIntegratorUnaryMinus
   ::LinearFormIntegratorUnaryMinus(const LinearFormIntegratorUnaryMinus& other)
      : m_lfi(other.m_lfi->copy())
   {}

   LinearFormIntegratorBase& LinearFormIntegratorUnaryMinus::getLFI()
   {
      return *m_lfi;
   }

   void LinearFormIntegratorUnaryMinus::buildMFEMLinearFormIntegrator()
   {
      m_lfi->buildMFEMLinearFormIntegrator();
      m_mfemLFI = std::make_unique<Internal::LFIUnaryMinus>(
            m_lfi->getMFEMLinearFormIntegrator());
   }

   mfem::LinearFormIntegrator&
   LinearFormIntegratorUnaryMinus::getMFEMLinearFormIntegrator()
   {
      assert(m_mfemLFI);
      return *m_mfemLFI;
   }

   mfem::LinearFormIntegrator*
   LinearFormIntegratorUnaryMinus::releaseMFEMLinearFormIntegrator()
   {
      return m_mfemLFI.release();
   }

   LinearFormIntegratorUnaryMinus operator-(const LinearFormIntegratorBase& lfi)
   {
      return LinearFormIntegratorUnaryMinus(lfi);
   }
}
