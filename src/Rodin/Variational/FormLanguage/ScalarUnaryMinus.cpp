#include "Rodin/Variational/ScalarCoefficient.h"

#include "ScalarSum.h"

#include "ScalarUnaryMinus.h"

namespace Rodin::Variational::FormLanguage
{
   ScalarCoefficientUnaryMinus
   ::ScalarCoefficientUnaryMinus(const ScalarCoefficientBase& s)
      : m_s(s.copy())
   {}

   ScalarCoefficientUnaryMinus
   ::ScalarCoefficientUnaryMinus(const ScalarCoefficientUnaryMinus& other)
      : m_s(other.m_s->copy())
   {}

   ScalarCoefficientBase&
   ScalarCoefficientUnaryMinus::getScalarCoefficient()
   {
      return *m_s;
   }

   void ScalarCoefficientUnaryMinus::buildMFEMCoefficient()
   {
      m_s->buildMFEMCoefficient();
      m_mfemCoefficient.emplace(0, m_s->getMFEMCoefficient(), 0, -1.0);
   }

   mfem::Coefficient& ScalarCoefficientUnaryMinus::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }

   ScalarCoefficientUnaryMinus operator-(const ScalarCoefficientBase& s)
   {
      return ScalarCoefficientUnaryMinus(s);
   }

   ScalarSum
      operator-(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return ScalarSum(
            lhs, ScalarCoefficientUnaryMinus(rhs));
   }
}
