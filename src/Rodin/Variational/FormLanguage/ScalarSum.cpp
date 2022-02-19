#include "Rodin/Variational/ScalarCoefficient.h"

#include "ScalarSum.h"

namespace Rodin::Variational::FormLanguage
{
   ScalarSum::ScalarSum(
         const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
      : m_lhs(lhs.copy()), m_rhs(rhs.copy())
   {}

   ScalarSum::ScalarSum(const ScalarSum& other)
      : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
   {}

   ScalarCoefficientBase& ScalarSum::getLHS()
   {
      assert(m_lhs);
      return *m_lhs;
   }

   ScalarCoefficientBase& ScalarSum::getRHS()
   {
      assert(m_rhs);
      return *m_rhs;
   }

   double ScalarSum::getValue(
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) 
   {
      return m_lhs->getValue(trans, ip) + m_rhs->getValue(trans, ip);
   }

   ScalarSum
      operator+(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return ScalarSum(lhs, rhs);
   }
}

