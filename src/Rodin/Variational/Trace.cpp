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

   void Trace::build()
   {
      m_matrix->build();
      m_mfemCoefficient.emplace(m_matrix->get());
   }

   mfem::Coefficient& Trace::get()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }
}
