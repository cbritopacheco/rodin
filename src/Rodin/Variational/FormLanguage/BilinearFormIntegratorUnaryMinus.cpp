#include "BilinearFormIntegratorUnaryMinus.h"

namespace Rodin::Variational::FormLanguage
{
   // ---- BilinearFormIntegratorSum -----------------------------------------
   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorSum>
   ::BilinearFormIntegratorUnaryMinus(const BilinearFormIntegratorSum& lfi)
      : BilinearFormIntegratorSum(lfi)
   {
      for (auto& p : getBilinearFormDomainIntegratorList())
         p.reset(new BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase>(*p));
   }

   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorSum>
   ::BilinearFormIntegratorUnaryMinus(const BilinearFormIntegratorUnaryMinus& other)
      : BilinearFormIntegratorSum(other)
   {}

   // ---- BilinearFormIntegratorBase ----------------------------------------
   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase>
   ::BilinearFormIntegratorUnaryMinus(const BilinearFormIntegratorBase& bfi)
      : m_bfi(bfi.copy())
   {}

   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase>
   ::BilinearFormIntegratorUnaryMinus(
         const BilinearFormIntegratorUnaryMinus& other)
      : m_bfi(other.copy())
   {}

   BilinearFormIntegratorBase& 
   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase>
   ::getBFI()
   {
      return *m_bfi;
   }

   void
   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase>
   ::buildMFEMBilinearFormIntegrator()
   {
      m_bfi->buildMFEMBilinearFormIntegrator();
      m_mfemBFI = std::make_unique<Internal::BFIUnaryMinus>(
            m_bfi->getMFEMBilinearFormIntegrator());
   }

   mfem::BilinearFormIntegrator&
   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase>
   ::getMFEMBilinearFormIntegrator()
   {
      assert(m_mfemBFI);
      return *m_mfemBFI;
   }

   mfem::BilinearFormIntegrator*
   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase>
   ::releaseMFEMBilinearFormIntegrator()
   {
      return m_mfemBFI.release();
   }

   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase>
   operator-(const BilinearFormIntegratorBase& bfi)
   {
      return BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase>(bfi);
   }

   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorSum>
   operator-(const BilinearFormIntegratorSum& lfi)
   {
      return BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorSum>(lfi);
   }
}
