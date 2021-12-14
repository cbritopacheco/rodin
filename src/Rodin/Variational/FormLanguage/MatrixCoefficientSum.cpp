#include "Rodin/Variational/MatrixCoefficient.h"

#include "MatrixCoefficientSum.h"

namespace Rodin::Variational::FormLanguage
{
   MatrixCoefficientSum::MatrixCoefficientSum(
         const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
      : m_lhs(lhs.copy()), m_rhs(rhs.copy())
   {
      assert(lhs.getRows() == rhs.getRows());
      assert(lhs.getColumns() == rhs.getColumns());
   }

   MatrixCoefficientSum::MatrixCoefficientSum(const MatrixCoefficientSum& other)
      : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
   {}

   MatrixCoefficientBase& MatrixCoefficientSum::getLHS()
   {
      return *m_lhs;
   }

   MatrixCoefficientBase& MatrixCoefficientSum::getRHS()
   {
      return *m_rhs;
   }

   int MatrixCoefficientSum::getRows() const
   {
      return m_lhs->getRows();
   }

   int MatrixCoefficientSum::getColumns() const
   {
      return m_lhs->getColumns();
   }

   void MatrixCoefficientSum::buildMFEMMatrixCoefficient()
   {
      m_lhs->buildMFEMMatrixCoefficient();
      m_rhs->buildMFEMMatrixCoefficient();

      m_mfemMatrixCoefficient.emplace(
            m_lhs->getMFEMMatrixCoefficient(), m_rhs->getMFEMMatrixCoefficient());
   }

   mfem::MatrixCoefficient& MatrixCoefficientSum::getMFEMMatrixCoefficient()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }

   MatrixCoefficientSum
   operator+(const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
   {
      return MatrixCoefficientSum(lhs, rhs);
   }
}
