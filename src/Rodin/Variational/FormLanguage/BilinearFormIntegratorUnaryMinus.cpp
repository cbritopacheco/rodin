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
