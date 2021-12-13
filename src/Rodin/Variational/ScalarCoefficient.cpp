#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   // ---- CoeffUnaryMinus ---------------------------------------------------
   // ------------------------------------------------------------------------
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus>
   ::ScalarCoefficient(const FormLanguage::ScalarCoefficientUnaryMinus& expr)
      : m_expr(expr.copy())
   {}

   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus>
   ::ScalarCoefficient(const ScalarCoefficient& other)
      :  m_expr(other.m_expr->copy())
   {}

   void
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus>::buildMFEMCoefficient()
   {
      m_expr->getScalarCoefficient().buildMFEMCoefficient();
      m_mfemCoefficient.emplace(
            0, m_expr->getScalarCoefficient().getMFEMCoefficient(),
            0, -1.0);
   }

   mfem::Coefficient&
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus>::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }

   // ---- FormLanguage::ScalarCoefficientSum --------------------------------
   // ------------------------------------------------------------------------
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum>
   ::ScalarCoefficient(const FormLanguage::ScalarCoefficientSum& expr)
      : m_expr(expr.copy())
   {}

   ScalarCoefficient<FormLanguage::ScalarCoefficientSum>
   ::ScalarCoefficient(const ScalarCoefficient& other)
      :  m_expr(other.m_expr->copy()),
         m_mfemCoefficient(other.m_mfemCoefficient)
   {}

   void
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum>::buildMFEMCoefficient()
   {
      m_expr->getLHS().buildMFEMCoefficient();
      m_expr->getRHS().buildMFEMCoefficient();

      m_mfemCoefficient.emplace(
            m_expr->getLHS().getMFEMCoefficient(),
            m_expr->getRHS().getMFEMCoefficient());
   }

   mfem::Coefficient&
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum>::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }
}
