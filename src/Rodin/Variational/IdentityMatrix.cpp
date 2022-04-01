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
}
