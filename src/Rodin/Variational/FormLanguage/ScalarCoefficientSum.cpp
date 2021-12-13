#include "Rodin/Variational/ScalarCoefficient.h"

#include "ScalarCoefficientSum.h"

namespace Rodin::Variational::FormLanguage
{
   ScalarCoefficientSum::ScalarCoefficientSum(
         const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
      : m_lhs(lhs.copy()), m_rhs(rhs.copy())
   {}

   ScalarCoefficientSum::ScalarCoefficientSum(const ScalarCoefficientSum& other)
      : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
   {}

   ScalarCoefficientBase& ScalarCoefficientSum::getLHS()
   {
      assert(m_lhs);
      return *m_lhs;
   }

   ScalarCoefficientBase& ScalarCoefficientSum::getRHS()
   {
      assert(m_rhs);
      return *m_rhs;
   }

   ScalarCoefficientSum
      operator+(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return ScalarCoefficientSum(lhs, rhs);
   }
}

