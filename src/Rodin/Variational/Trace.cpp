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
      :  m_matrix(other.m_matrix->copy()),
         m_mfemCoefficient(other.m_mfemCoefficient)
   {}

   void Trace::buildMFEMCoefficient()
   {
      m_matrix->buildMFEMMatrixCoefficient();
      m_mfemCoefficient.emplace(m_matrix->getMFEMMatrixCoefficient());
   }

   mfem::Coefficient& Trace::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }
}
