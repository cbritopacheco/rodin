#include "Rodin/Variational/ScalarCoefficient.h"

#include "ScalarCoefficientUnaryMinus.h"

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

   ScalarCoefficientUnaryMinus operator-(const ScalarCoefficientBase& s)
   {
      return ScalarCoefficientUnaryMinus(s);
   }
}
