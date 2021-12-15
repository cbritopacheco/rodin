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

   void ScalarSum::buildMFEMCoefficient()
   {
      m_lhs->buildMFEMCoefficient();
      m_rhs->buildMFEMCoefficient();
      m_mfemCoefficient.emplace(
            m_lhs->getMFEMCoefficient(),
            m_rhs->getMFEMCoefficient());
   }

   mfem::Coefficient& ScalarSum::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }

   ScalarSum
      operator+(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return ScalarSum(lhs, rhs);
   }
}

