#include "Gradient.h"

#include "GridFunction.h"
#include "H1.h"

namespace Rodin::Variational
{
   size_t Gradient::getDimension() const
   {
      return m_u.getFiniteElementSpace().getMesh().getDimension();
   }

   void Gradient::build()
   {
      m_mfemVectorCoefficient.emplace(&m_u.getHandle());
   }

   mfem::VectorCoefficient& Gradient::get()
   {
      assert(m_mfemVectorCoefficient);
      return *m_mfemVectorCoefficient;
   }
}

