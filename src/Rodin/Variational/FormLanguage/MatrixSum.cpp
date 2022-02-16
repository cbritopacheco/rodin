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

   void MatrixSum::build()
   {
      m_lhs->build();
      m_rhs->build();

      m_mfemMatrixCoefficient.emplace(
            m_lhs->get(), m_rhs->get());
   }

   mfem::MatrixCoefficient& MatrixSum::get()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }

   MatrixSum
   operator+(const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
   {
      return MatrixSum(lhs, rhs);
   }
}
