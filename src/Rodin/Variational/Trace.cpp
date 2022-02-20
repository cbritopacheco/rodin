#include "MatrixCoefficient.h"

#include "Trace.h"

namespace Rodin::Variational
{
   Trace::Trace(const MatrixCoefficientBase& m)
      : m_matrix(m.copy())
   {
      assert(m.getColumns() == m.getRows());
   }

   Trace::Trace(const Trace& other)
      :  m_matrix(other.m_matrix->copy())
   {}

   double Trace::getValue(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   {
      mfem::DenseMatrix mat;
      m_matrix->getValue(mat, trans, ip);
      return mat.Trace();
   }
}
