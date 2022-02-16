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

   void Jacobian::build()
   {
      m_mfemMatrixCoefficient.emplace(m_u.getHandle());
   }

   mfem::MatrixCoefficient& Jacobian::get()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }
}
