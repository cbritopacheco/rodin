#include "../ScalarCoefficient.h"

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
         m_matrix(other.m_matrix->copy())
   {}

   int ScalarMatrixProduct::getRows() const
   {
      return m_matrix->getRows();
   }

   int ScalarMatrixProduct::getColumns() const
   {
      return m_matrix->getColumns();
   }

   void ScalarMatrixProduct::getValue(
         mfem::DenseMatrix& value,
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   {
      m_matrix->getValue(value, trans, ip);
      value *= m_scalar->getValue(trans, ip);
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
