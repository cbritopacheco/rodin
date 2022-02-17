#include "H1.h"

#include "Jacobian.h"

namespace Rodin::Variational
{
   Jacobian& Jacobian::setTraceDomain(int traceDomain)
   {
      m_traceDomain = traceDomain;
      return *this;
   }

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
      if (m_traceDomain)
         m_mfemMatrixCoefficient.emplace(m_u.getHandle(), *m_traceDomain);
      else
         m_mfemMatrixCoefficient.emplace(m_u.getHandle());
   }

   mfem::MatrixCoefficient& Jacobian::get()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }
}
