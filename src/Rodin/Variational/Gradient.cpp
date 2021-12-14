#include "H1.h"

#include "Gradient.h"

namespace Rodin::Variational
{
   Gradient::Gradient(GridFunction<H1>& u)
      :  m_u(u)
   {}

   Gradient::Gradient(const Gradient& other)
      : m_u(other.m_u)
   {}

   int Gradient::getRows() const
   {
      return m_u.getFiniteElementSpace().getDimension();
   }

   int Gradient::getColumns() const
   {
      return m_u.getFiniteElementSpace().getMesh().getDimension();
   }

   void Gradient::buildMFEMMatrixCoefficient()
   {
      m_mfemMatrixCoefficient.emplace(m_u.getHandle());
   }

   mfem::MatrixCoefficient& Gradient::getMFEMMatrixCoefficient()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }
}
