#include "H1.h"

#include "Jacobian.h"

namespace Rodin::Variational
{
   Jacobian::Jacobian(GridFunction<H1>& u)
      :  m_u(u)
   {}

   Jacobian::Jacobian(const Jacobian& other)
      : m_u(other.m_u)
   {}

   int Jacobian::getRows() const
   {
      return m_u.getFiniteElementSpace().getDimension();
   }

   int Jacobian::getColumns() const
   {
      return m_u.getFiniteElementSpace().getMesh().getDimension();
   }

   void Jacobian::buildMFEMMatrixCoefficient()
   {
      m_mfemMatrixCoefficient.emplace(m_u.getHandle());
   }

   mfem::MatrixCoefficient& Jacobian::getMFEMMatrixCoefficient()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }
}
