#include "LinearFormIntegratorUnaryMinus.h"

namespace Rodin::Variational::FormLanguage
{
   LinearFormIntegratorUnaryMinus<LinearFormIntegratorSum>
   ::LinearFormIntegratorUnaryMinus(const LinearFormIntegratorSum& lfi)
      : LinearFormIntegratorSum(lfi)
   {
      for (auto& p : getLinearFormDomainIntegratorList())
         p.reset(new LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>(*p));
      for (auto& p : getLinearFormBoundaryIntegratorList())
         p.reset(new LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>(*p));
   }

   LinearFormIntegratorUnaryMinus<LinearFormIntegratorSum>
   ::LinearFormIntegratorUnaryMinus(const LinearFormIntegratorUnaryMinus& other)
      : LinearFormIntegratorSum(other)
   {}

   LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>
   ::LinearFormIntegratorUnaryMinus(const LinearFormIntegratorBase& lfi)
      : m_lfi(lfi.copy())
   {}

   LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>
   ::LinearFormIntegratorUnaryMinus(const LinearFormIntegratorUnaryMinus& other)
      : m_lfi(other.m_lfi->copy())
   {}

   LinearFormIntegratorBase&
   LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>
   ::getLFI()
   {
      return *m_lfi;
   }

   void
   LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>
   ::buildMFEMLinearFormIntegrator()
   {
      m_lfi->buildMFEMLinearFormIntegrator();
      m_mfemLFI = std::make_unique<Internal::LFIUnaryMinus>(
            m_lfi->getMFEMLinearFormIntegrator());
   }

   mfem::LinearFormIntegrator&
   LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>
   ::getMFEMLinearFormIntegrator()
   {
      assert(m_mfemLFI);
      return *m_mfemLFI;
   }

   mfem::LinearFormIntegrator*
   LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>
   ::releaseMFEMLinearFormIntegrator()
   {
      return m_mfemLFI.release();
   }

   LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>
   operator-(const LinearFormIntegratorBase& lfi)
   {
      return LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>(lfi);
   }

   LinearFormIntegratorUnaryMinus<LinearFormIntegratorSum>
   operator-(const LinearFormIntegratorSum& lfi)
   {
      return LinearFormIntegratorUnaryMinus<LinearFormIntegratorSum>(lfi);
   }
}
