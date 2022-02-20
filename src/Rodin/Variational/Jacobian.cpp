#include "H1.h"

#include "Jacobian.h"

namespace Rodin::Variational
{
   int Jacobian::getRows() const
   {
      return m_u.getFiniteElementSpace().getRangeDimension();
   }

   int Jacobian::getColumns() const
   {
      return m_u.getFiniteElementSpace().getMesh().getDimension();
   }
}
