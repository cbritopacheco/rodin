#include "Rodin/Variational/MatrixCoefficient.h"

#include "MatrixSum.h"

namespace Rodin::Variational::FormLanguage
{
   MatrixSum::MatrixSum(
         const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
      : m_lhs(lhs.copy()), m_rhs(rhs.copy())
   {
      assert(lhs.getRows() == rhs.getRows());
      assert(lhs.getColumns() == rhs.getColumns());
   }

   MatrixSum::MatrixSum(const MatrixSum& other)
      : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
   {}

   void MatrixSum::getValue(
         mfem::DenseMatrix& value,
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   {
      mfem::DenseMatrix m;
      m_lhs->getValue(m, trans, ip);
      m_rhs->getValue(value, trans, ip);
      value += m;
   }

   MatrixCoefficientBase& MatrixSum::getLHS()
   {
      return *m_lhs;
   }

   MatrixCoefficientBase& MatrixSum::getRHS()
   {
      return *m_rhs;
   }

   int MatrixSum::getRows() const
   {
      return m_lhs->getRows();
   }

   int MatrixSum::getColumns() const
   {
      return m_lhs->getColumns();
   }

   MatrixSum
   operator+(const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
   {
      return MatrixSum(lhs, rhs);
   }
}
