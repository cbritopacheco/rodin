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

   void ScalarCoefficientSum::buildMFEMCoefficient()
   {
      m_lhs->buildMFEMCoefficient();
      m_rhs->buildMFEMCoefficient();
      m_mfemCoefficient.emplace(
            m_lhs->getMFEMCoefficient(),
            m_rhs->getMFEMCoefficient());
   }

   mfem::Coefficient& ScalarCoefficientSum::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }

   ScalarCoefficientSum
      operator+(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return ScalarCoefficientSum(lhs, rhs);
   }
}

