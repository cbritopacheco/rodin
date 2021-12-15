#include "MatrixCoefficient.h"

#include "Trace.h"

namespace Rodin::Variational
{
   Trace::Trace(const MatrixCoefficientBase& m)
      : m_matrix(m.copy())
   {}

   Trace::Trace(const Trace& other)
      : m_matrix(other.m_matrix->copy())
   {}

   void Trace::buildMFEMCoefficient()
   {
      m_matrix->getMFEMMatrixCoefficient();
      m_mfemCoefficient.emplace(m_matrix->getMFEMMatrixCoefficient());
   }

   mfem::Coefficient& Trace::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }
}
