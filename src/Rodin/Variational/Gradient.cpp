#include "Gradient.h"

#include "GridFunction.h"
#include "H1.h"

namespace Rodin::Variational
{
   constexpr
   Gradient::Gradient(GridFunction<H1>& u)
      : m_u(u)
   {}

   size_t Gradient::getDimension() const
   {
      return m_u.getFiniteElementSpace().getMesh().getDimension();
   }

   void Gradient::buildMFEMVectorCoefficient()
   {
      m_mfemVectorCoefficient.emplace(&m_u.getHandle());
   }

   mfem::VectorCoefficient& Gradient::getMFEMVectorCoefficient()
   {
      assert(m_mfemVectorCoefficient);
      return *m_mfemVectorCoefficient;
   }
}

