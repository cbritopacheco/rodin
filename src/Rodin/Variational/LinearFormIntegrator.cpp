#include "LinearFormIntegrator.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      LinearFormIntegrator::LinearFormIntegrator(const LinearFormIntegratorBase& lfi)
         : m_lfi(lfi.copy())
      {}

      void LinearFormIntegrator::AssembleRHSElementVect(
            const mfem::FiniteElement& fe,
            mfem::ElementTransformation& trans,
            mfem::Vector& vec)
      {
         m_lfi->getElementVector(fe, trans, vec);
      }
   }
}
