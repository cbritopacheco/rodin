#include "IdentityMatrix.h"

namespace Rodin::Variational
{
   int IdentityMatrix::getRows() const
   {
      return m_n;
   }

   int IdentityMatrix::getColumns() const
   {
      return m_n;
   }

   void IdentityMatrix::build()
   {
      m_mfemMatrixCoefficient.emplace(m_n);
   }

   mfem::MatrixCoefficient& IdentityMatrix::get()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }
}
