#include "ScalarMatrixProduct.h"

namespace Rodin::Variational::FormLanguage
{
   ScalarMatrixProduct::ScalarMatrixProduct(
         const ScalarCoefficientBase& s, const MatrixCoefficientBase& m)
      :  m_scalar(s.copy()),
         m_matrix(m.copy())
   {}

   ScalarMatrixProduct::ScalarMatrixProduct(const ScalarMatrixProduct& other)
      :  m_scalar(other.m_scalar->copy()),
         m_matrix(other.m_matrix->copy()),
         m_mfemMatrixCoefficient(other.m_mfemMatrixCoefficient)
   {}

   int ScalarMatrixProduct::getRows() const
   {
      return m_matrix->getRows();
   }

   int ScalarMatrixProduct::getColumns() const
   {
      return m_matrix->getColumns();
   }

   void ScalarMatrixProduct::buildMFEMMatrixCoefficient()
   {
      m_scalar->buildMFEMCoefficient();
      m_matrix->buildMFEMMatrixCoefficient();
      m_mfemMatrixCoefficient.emplace(
            m_scalar->getMFEMCoefficient(), m_matrix->getMFEMMatrixCoefficient());
   }

   mfem::MatrixCoefficient& ScalarMatrixProduct::getMFEMMatrixCoefficient()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }

   ScalarMatrixProduct
   operator*(const ScalarCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
   {
      return ScalarMatrixProduct(lhs, rhs);
   }

   ScalarMatrixProduct
   operator*(const MatrixCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return ScalarMatrixProduct(rhs, lhs);
   }
}
