#include "ScalarProduct.h"

namespace Rodin::Variational::FormLanguage
{
   ScalarProduct::ScalarProduct(
         const ScalarCoefficientBase& a, const ScalarCoefficientBase& b)
      : m_a(a.copy()), m_b(b.copy())
   {}

   ScalarProduct::ScalarProduct(const ScalarProduct& other)
      :  m_a(other.m_a->copy()), m_b(other.m_b->copy()),
         m_mfemCoefficient(other.m_mfemCoefficient)
   {}

   void ScalarProduct::buildMFEMCoefficient()
   {
      m_a->buildMFEMCoefficient();
      m_b->buildMFEMCoefficient();
      m_mfemCoefficient.emplace(m_a->getMFEMCoefficient(), m_b->getMFEMCoefficient());
   }

   mfem::Coefficient& ScalarProduct::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }

   ScalarProduct
   operator*(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return ScalarProduct(lhs, rhs);
   }
}
