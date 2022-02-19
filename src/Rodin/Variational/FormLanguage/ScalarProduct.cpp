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

   double ScalarProduct::getValue(
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   {
      return m_a->getValue(trans, ip) * m_b->getValue(trans, ip);
   }

   ScalarProduct
   operator*(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return ScalarProduct(lhs, rhs);
   }
}
