#include "ScalarProduct.h"

namespace Rodin::Variational::FormLanguage
{
   ScalarProduct::ScalarProduct(
         const ScalarCoefficientBase& a, const ScalarCoefficientBase& b)
      : m_a(a.copy()), m_b(b.copy())
   {}

   ScalarProduct::ScalarProduct(const ScalarProduct& other)
      :  m_a(other.m_a->copy()), m_b(other.m_b->copy())
   {}

   void ScalarProduct::build()
   {
      m_a->build();
      m_b->build();
      m_mfemCoefficient.emplace(m_a->get(), m_b->get());
   }

   mfem::Coefficient& ScalarProduct::get()
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
