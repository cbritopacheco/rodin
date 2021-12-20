#include "BilinearFormIntegratorUnaryMinus.h"

namespace Rodin::Variational::FormLanguage
{
   BilinearFormIntegratorUnaryMinus
   ::BilinearFormIntegratorUnaryMinus(const BilinearFormIntegratorBase& bfi)
      : m_bfi(bfi.copy())
   {}

   BilinearFormIntegratorUnaryMinus
   ::BilinearFormIntegratorUnaryMinus(
         const BilinearFormIntegratorUnaryMinus& other)
      : m_bfi(other.copy())
   {}

   BilinearFormIntegratorBase& BilinearFormIntegratorUnaryMinus::getBFI()
   {
      return *m_bfi;
   }

   void BilinearFormIntegratorUnaryMinus::buildMFEMBilinearFormIntegrator()
   {
      m_bfi->buildMFEMBilinearFormIntegrator();
      m_mfemBFI = std::make_unique<Internal::BFIUnaryMinus>(
            m_bfi->getMFEMBilinearFormIntegrator());
   }

   mfem::BilinearFormIntegrator&
   BilinearFormIntegratorUnaryMinus::getMFEMBilinearFormIntegrator()
   {
      assert(m_mfemBFI);
      return *m_mfemBFI;
   }

   mfem::BilinearFormIntegrator*
   BilinearFormIntegratorUnaryMinus::releaseMFEMBilinearFormIntegrator()
   {
      return m_mfemBFI.release();
   }

   BilinearFormIntegratorUnaryMinus operator-(
         const BilinearFormIntegratorBase& bfi)
   {
      return BilinearFormIntegratorUnaryMinus(bfi);
   }
}
