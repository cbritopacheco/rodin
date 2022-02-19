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

   double ScalarCoefficientUnaryMinus::getValue(
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   {
      return -1.0 * m_s->getValue(trans, ip);
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
