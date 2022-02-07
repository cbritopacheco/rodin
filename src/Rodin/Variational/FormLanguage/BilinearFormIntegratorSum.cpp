#include "BilinearFormIntegratorUnaryMinus.h"

#include "BilinearFormIntegratorSum.h"

namespace Rodin::Variational::FormLanguage
{
   BilinearFormIntegratorSum
   ::BilinearFormIntegratorSum(
         const BilinearFormIntegratorBase& lhs,
         const BilinearFormIntegratorBase& rhs)
      : m_lhs(lhs.copy()), m_rhs(rhs.copy())
   {}

   BilinearFormIntegratorSum
   ::BilinearFormIntegratorSum(const BilinearFormIntegratorSum& other)
      : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
   {
      assert(std::equal(
               lhs.getAttributes().begin(), lhs.getAttributes().end(),
               rhs.getAttributes().begin()
               ));
   }

   BilinearFormIntegratorBase& BilinearFormIntegratorSum::getLHS()
   {
      return *m_lhs;
   }

   BilinearFormIntegratorBase& BilinearFormIntegratorSum::getRHS()
   {
      return *m_rhs;
   }

   void BilinearFormIntegratorSum::buildMFEMBilinearFormIntegrator()
   {
      m_lhs->buildMFEMBilinearFormIntegrator();
      m_rhs->buildMFEMBilinearFormIntegrator();
      m_mfemBFI = std::make_unique<Internal::BilinearFormIntegratorSum>(
            m_lhs->getMFEMBilinearFormIntegrator(),
            m_rhs->getMFEMBilinearFormIntegrator());
   }

   mfem::BilinearFormIntegrator&
   BilinearFormIntegratorSum::getMFEMBilinearFormIntegrator()
   {
      assert(m_mfemBFI);
      return *m_mfemBFI;
   }

   mfem::BilinearFormIntegrator*
   BilinearFormIntegratorSum::releaseMFEMBilinearFormIntegrator()
   {
      return m_mfemBFI.release();
   }

   BilinearFormIntegratorSum operator+(
         const BilinearFormIntegratorBase& lhs,
         const BilinearFormIntegratorBase& rhs)
   {
      return BilinearFormIntegratorSum(lhs, rhs);
   }

   BilinearFormIntegratorSum operator-(
         const BilinearFormIntegratorBase& lhs,
         const BilinearFormIntegratorBase& rhs)
   {
      return BilinearFormIntegratorSum(lhs, BilinearFormIntegratorUnaryMinus(rhs));
   }
}
