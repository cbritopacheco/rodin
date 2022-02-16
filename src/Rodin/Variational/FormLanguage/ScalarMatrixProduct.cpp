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

   void ScalarMatrixProduct::build()
   {
      m_scalar->build();
      m_matrix->build();
      m_mfemMatrixCoefficient.emplace(
            m_scalar->get(), m_matrix->get());
   }

   mfem::MatrixCoefficient& ScalarMatrixProduct::get()
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
