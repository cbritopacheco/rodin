#include "IdentityMatrix.h"

namespace Rodin::Variational
{
   IdentityMatrix::IdentityMatrix(int n)
      : m_n(n)
   {}

   IdentityMatrix::IdentityMatrix(const IdentityMatrix& other)
      : m_n(other.m_n), m_mfemMatrixCoefficient(other.m_mfemMatrixCoefficient)
   {}

   int IdentityMatrix::getRows() const
   {
      return m_n;
   }

   int IdentityMatrix::getColumns() const
   {
      return m_n;
   }

   void IdentityMatrix::buildMFEMMatrixCoefficient()
   {
      m_mfemMatrixCoefficient.emplace(m_n);
   }

   mfem::MatrixCoefficient& IdentityMatrix::getMFEMMatrixCoefficient()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }
}
